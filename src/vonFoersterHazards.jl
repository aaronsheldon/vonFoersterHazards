"""
    vonFoersterHazards

Engine and utilities to forward propagate the von Foerster evolution
from the hazard rate, birth rate, and initial demographics. Remember
that our transition matrix convention is that SOURCES are columns 
and TARGETS are rows.
"""
module vonFoersterHazards

import Base.iterate,
       Base.length,
       Base.eltype,
       Base.firstindex,
       Base.lastindex,
       Base.getindex,
       Base.setindex!,
       Base.size

export randomtruncate,
       conservesum!,
       covariance,
       abstractcohort,
       abstractpopulation,
       abstractevolve,
       cohort,
       population,
       evolve,
       scatterrate,
       hazardrate,
       birthrate

"""
    abtractcohort

Parametric container type for concrete cohort types.
"""
abstract type abstractcohort{
    Q<:Real,
    R<:Real,
    S<:AbstractVector{U} where U<:Real,
    T<:AbstractMatrix{V} where V<:Real
}
end

"""
    cohort(elapsed, age, stratum, covariance, conserving)

Return type for indexing into the population. Container for the state occupancies
of a specific age group, and the covariances between the occupancies. The occupancies
are conserved in groups in boundaries indicated by the true values of conserving.
"""
struct cohort{Q, R, S, T} <: abstractcohort{Q, R, S, T}
    elapsed::Q
    age::R
    stratum::S
    covariance::T
    conserving::BitArray{1}
end
cohort(a, s) = cohort(
    0.0,
    a,
    s,
    zeros(Float64, size(s, 1), size(s, 1)),
    BitArray(zeros(size(s, 1)))
)

"""
    abstractpopulation

Parametric container type for concrete population types.
"""
abstract type abstractpopulation{
    Q<:Real, 
    R<:AbstractVector{U} where U<:Real,
    S<:AbstractMatrix{V} where V<:Real,
    T<:AbstractArray{W, 3} where W<:Real
}
end

"""
    population(elapsed, ages, strata, covariances, conserving)

Indexable container for the demographics of a population at a single time step.
Conceptually the ages labels the rows of the population matrix, and the columns
are the states. Cohorts are stored in reverse order, so that births can be pushed
to the vector. The occupancies are conserved in groups in boundaries indicated by
the true values of conserving.
"""
struct population{Q, R, S, T} <: abstractpopulation{Q, R, S, T}
    elapsed::Q
    ages::R
    strata::S
    covariances::T
    conserving::BitArray{1}
end
population(a, s) = population(
    0.0,
    a,
    s,
    zeros(Float64, size(s, 1), size(s, 2), size(s, 2)),
    BitArray(zeros(size(s, 2)))
)
Base.eltype(::Type{A}) where {
    Q<:Real, 
    U<:Real,
    R<:AbstractVector{U},
    V<:Real,
    S<:AbstractMatrix{V},
    W<:Real,
    T<:AbstractArray{W, 3},
    A<:abstractpopulation{Q, R, S, T}
} = cohort{
    Q, 
    U, 
    SubArray{V, 1, S, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true},
    SubArray{W, 2, T, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}}, true}
}
Base.length(P::abstractpopulation) = length(P.ages)
Base.size(P::abstractpopulation, d=1) = ((d==1) ? length(P) : 1)
Base.firstindex(P::abstractpopulation) = 1
Base.lastindex(P::abstractpopulation) = length(P)
Base.getindex(P::abstractpopulation, i) = cohort(
    P.elapsed,
    P.ages[i],
    view(P.strata, i, :),
    view(P.covariances, i, :, :),
    P.conserving
)
function Base.setindex!(P::abstractpopulation, C::abstractcohort, i)
    P.ages[i] = C.age
    P.strata[i, :] .= C.stratum[:]
    P.covariances[i, :, :] .= C.covariance[:, :]
    C
end

"""
    abstractevolve

Parametric container type for concrete evolve types.
"""
abstract type abstractevolve{
    Q<:abstractpopulation{T, U, V, W} where {
        T<:Real, 
        X<:Real,
        U<:AbstractVector{X},
        Y<:Real,
        V<:AbstractMatrix{Y},
        Z<:Real,
        W<:AbstractArray{Z, 3}
    },
    R<:Real,
    S<:Real
}
end

"""
    evolve(initial, size, count, gestation)

Iterable container for the population evolution engine. Computes count steps of
duration size from starting population initial, generating a new cohort after
every gestation has elapsed
"""
struct evolve{R, S, T} <: abstractevolve{R, S, T}
    initial::R
    size::S
    gestation::T
    count::Int64
end
Base.eltype(::Type{A}) where {
    T<:Real, 
    X<:Real,
    U<:AbstractVector{X},
    Y<:Real,
    V<:AbstractMatrix{Y},
    Z<:Real,
    W<:AbstractArray{Z, 3},
    Q<:abstractpopulation{T, U, V, W},
    R<:Real,
    S<:Real,
    A<:abstractevolve{Q, R, S}
} = Q
Base.length(E::abstractevolve) = E.count
Base.size(E::abstractevolve, d=1) = ((d==1) ? length(E) : 1)
function Base.iterate(E::abstractevolve)

    # Sanity checks
    (E.size > 0) ||
        throw(DomainError(E.size, "time increment must be positive."))

    (E.gestation > 0) ||
        throw(DomainError(E.gestation, "gestation must be positive."))

    (E.count > 0) ||
        throw(DomainError(E.gestation, "count must be positive."))

    (issorted(E.initial.ages, rev=true)) ||
        throw(DomainError(E.initial.ages, "ages must be sorted in descending order."))

    (0 <= E.initial.ages[end]) ||
        throw(DomainError(E.initial.ages, "ages must be non-negative."))

    (0 <= minimum(E.initial.strata)) ||
        throw(DomainError(E.initial.strata, "cohorts must be non-negative."))

    (minimum(E.initial.covariances) == 0 == maximum(E.initial.covariances)) ||
        throw(DomainError(E.initial.covariances, "covariances must be zero."))

    (al,) = size(E.initial.ages)
    (bl,) = size(birthrate(E.initial))
    (cl, cw) = size(E.initial.strata)
    (dl,) = size(scatterrate(E.initial.ages[1]))
    (el, ew, eh) = size(hazardrate(E.initial))
    (fl, fw, fh) = size(E.initial.covariances)
    (gl,) = size(E.initial.conserving)

    (al == cl == fl) ||
        throw(DimensionMismatch("number of cohorts and covariances must equal number of ages."))

    (dl == el) ||
        throw(DimensionMismatch("number of summands in spectral decomposition returned by hazard rate and scatter rate must be equal."))

    (bl  == cw == ew == eh == fw == fh == gl) ||
        throw(DimensionMismatch("strata in covariances, birth rate, hazard rate, conserving, and cohorts must be equal."))

    # Send
    (E.initial, (0, E.initial))
end
function Base.iterate(E::abstractevolve, S)
    (S[1] <= E.count) || return nothing
    
    # One time computation of the exogenous birth rate
    b = birthrate(S[2])

    # One time computation of the exogenous hazard rates
    H = hazardrate(S[2])

# # # Curried transition function to update the cohort. # # # # # # # # # # # # # # #
    function transitioncohort(c)                                                    #
        P = exp(-E.size * sum(scatterrate(c.age) .* H, dims=1)[1, :, :])            #
        cohort(                                                                     #
            c.elapsed + E.size,                                                     #
            c.age + E.size,                                                         #
            conservesum!(randomtruncate.(P * c.stratum), c.stratum, c.conserving),  #
            covariance(c.covariance, P, c.stratum),                                 #
            c.conserving                                                            #
        )                                                                           #
    end                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                                                                                    #
# # # SIMD main loop, compute the transitions within each cohort. # # # # # # # # # #
    S[2] .= transitioncohort.(S[2])                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Pass the pointers, youngest cohort is less than gestation time old, keep the cohorts
    if S[2].ages[end] < E.gestation
        R = population(
            S[2] + E.size,
            S[2].ages,
            S[2].strata,
            S[2].covariances,
            S[2].conserving
        )

    # Pass the pointers, youngest cohort is more than gestation time old, add a cohort
    else
        R = population(
            S[2] + E.size,
            [S[2].ages; zero(eltype(S[2].ages))],
            [S[2].strata; zeros(eltype(S[2].strata), 1, size(S[2].strata, 2))],
            [S[2].covariances; zeros(eltype(S[2].covariances), 1, size(S[2].covariances, 2), size(S[2].covariances, 3))],
            S[2].conserving
        )
    end
    
    # Add births to the youngest cohort
    R.strata[end, :] .= R.strata[end, :] .+ b'
    
    # Send
    (R, (1 + S[1], R))
end

"""
    birthrate(P)

Stub function to be overloaded in implementation. Compute the exogenous
birth rate vector from the population P.
"""
function birthrate(P::abstractpopulation) end

"""
    scatterrate(a)

Stub function to be overloaded in implementation. For the age a return
a vector where the elements correspond to the spectral decomposition of
the endogeneous scattering rate.
"""
function scatterrate(a::Real) end

"""
    hazardrate(P)

Stub function to be overloaded in implementation. For the population P return
a three dimensional array where the length corresponds to the spectral
decomposition and the height and width are square corresponding to the
exogenous hazard rate.
"""
function hazardrate(P::abstractpopulation) end

"""
    randomtruncate(x)

Randomly return the floor or the ceiling by comparing the fractional part of the
number to a sample from the uniform distribution on [0,1). If the fraction is
greater than the random sample return the ceiling otherwise return the floor.
"""
function randomtruncate(x)
    y = convert(Int64, trunc(x))
    y + convert(Int64, rand() < (x-y))
end

"""
    conservesum!(A, a, c)

In place enforcement that the subset sums of A equals the subset sums of a, where
the boundary of the subsets are set by c. Note this assumes the inputs are well formed,
there are no bounds or sanity checks.
"""
function conservesum!(A, a, c)
    d = a - sum(A)
    I = sortperm(A, rev=(d<0))
    A[I] .= A[I] .+ [sign(d) .* ones(Int64, abs(d)) ; zeros(Int64, length(A)-abs(d))]
    A
end

"""
    covariance(C, P, n)

Update of the covariances C from the transition probabilities P and the previous
state n. Each column of P is assumed to represent the exit probabilities from a
single state and thus should add to 1 with all entries non-negative. Note this
assumes the inputs are well formed, there are no bounds or sanity checks.
"""
function covariance(C, P, n)
    U = (P.*n')*P'
    I = [1:1+stride(U, 2):length(U)...]
    U[I] .= (P*n) .- U[I]
    P*C*P' .+ U
end

end
