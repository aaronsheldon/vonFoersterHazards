"""
    vonFoersterHazards

Engine and utilities to forward propagate the von Foerster evolution
from the hazard rate, birth rate, and initial demographics.
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
       cohort,
       population,
       evolve,
       hazardrate,
       birthrate

"""
    cohort(elapsed, age, stratum, covariance)

Return type for indexing into the population. Container for the state occupancies
of a specific age group, and the covariances between the occupancies.
"""
struct cohort{
        Q<:Real,
        R<:Real,
        S<:AbstractVector{U} where U<:Real,
        T<:AbstractMatrix{V} where V<:Real
    }
    elapsed::Q
    age::R
    stratum::S
    covariance::T
end
cohort(a, o) = cohort(0.0, a, o, zeros(Float64, size(o, 1), size(o, 1)))

"""
    population(elapsed, ages, strata, covariances)

Indexable container for the demographics of a population at a single time step.
Conceptually the ages labels the rows of the population matrix, and the columns
are the states. Cohorts are stored in reverse order, so that births can be pushed
to the vector.
"""
struct population{
        Q<:Real, 
        R<:AbstractVector{U} where U<:Real,
        S<:AbstractMatrix{V} where V<:Real,
        T<:AbstractArray{W, 3} where W<:Real
    }
    elapsed::Q
    ages::R
    strata::S
    covariances::T
end
population(a, c) = population(0.0, a, c, zeros(Float64, size(c, 1), size(c, 2), size(c, 2)))
Base.eltype(::Type{population{Q, R{U}, S{V}, T{W}}}) = cohort{
    Q, 
    U, 
    SubArray{V, 1, S, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true},
    SubArray{W, 2, T, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}}, true}
}
Base.length(P::population) = length(P.ages)
Base.size(P::population, d=1) = ((d==1) ? length(P) : 1)
Base.firstindex(P::population) = 1
Base.lastindex(P::population) = length(P)
Base.getindex(P::population, i) = cohort(P.elapsed, P.ages[i], view(P.strata, i, :), view(P.covariances, i, :, :))
function Base.setindex!(P::population, C::cohort, i)
    P.ages[i] = C.age
    P.strata[i, :] .= C.stratum[:]
    P.covariances[i, :, :] .= C.covariance[:, :]
    C
end

"""
    evolve(initial, size, count)

Iterable container for the population evolution engine. Computes count steps of
duration size from starting population initial.
"""
struct evolve{S<:population, T<:Real}
    initial::S
    size::T
    count::Int64
end
Base.eltype(::Type{evolve{S, T}}) = S
Base.length(E::evolve) = E.count
Base.size(E::evolve, d=1) = ((d==1) ? length(E) : 1)
function Base.iterate(E::evolve)
    
    # Sanity checks
    (E.size > 0) ||
        throw(DomainError(size, "time increment must be positive."))
    
    (issorted(E.initial.ages, rev=true)) ||
        throw(DomainError(E.initial.ages, "ages must be sorted in descending order."))
    
    (convert(eltype(E.initial.ages), 0) <= E.initial.ages[end]) ||
        throw(DomainError(E.initial.ages, "ages must be non-negative."))
    
    (convert(eltype(E.initial.strata), 0) <= minimum(E.initial.strata)) ||
        throw(DomainError(E.initial.strata, "cohorts must be non-negative."))
    
    (minimum(E.initial.covariances) == convert(eltype(E.initial.covariances), 0) == maximum(E.initial.covariances)) ||
        throw(DomainError(E.initial.covariances, "covariances must be zero."))
    
    H = hazardrate(E.initial)
    (al,) = size(E.initial.ages)
    (bl,) = size(birthrate(E.initial))
    (cl, cw) = size(E.initial.strata)
    (dl, dw) = size(H)
    (el, ew) = size(hazardrate(E.initial.ages[1], H))
    (fl, fw, fh) = size(E.initial.covariances)
    
    (al == cl == fl) ||
        throw(DimensionMismatch("number of cohorts and covariances must equal number of ages."))
    
    (bl  == cw == dl == dw == el == ew == fw == fh) ||
        throw(DimensionMismatch("strata in covariances, birth rate, extensize hazard rate, intensize hazard rate, and cohorts must be equal."))
    
    # Send
    (E.initial, (0, E.initial))
end
function Base.iterate(E::evolve, S)
    (S[1] <= E.count) || return nothing
        
    # One time computation of the extensize birth rate
    b = birthrate(S[2])

    # One time computation of the extensive hazard rate
    H = hazardrate(S[2])

    # Curried transition function
    function t(c)
        P = exp(-E.size * (hazardrate(a) + scatterrate(a) * H))
        cohort(
            c.elapsed + E.size,
            c.age + E.size,
            conservesum!(randomtruncate.(P * c.stratum), sum(c.stratum)),
            covariance(c.covariance, P, c.stratum)
        )
    end

    # Compute the transitions within each cohort
    S[2] .= t.(S[2])

    # Youngest cohort is less than 1 year old, add births to youngest cohort
    if S[1].ages[end] < convert(eltype(S[1].ages), 1) then
        R = population(
            S[2] + E.size,
            S[2].ages,
            S[2].strata,
            S[2].covariances
        )
        R.cohorts[end, :] .= R.cohorts[end, :] .+ b'

    # Youngest cohort is more than 1 year old, generate a new youngest cohort
    else
        R = population(
            S[2] + E.size,
            S[2].ages,
            [S[2].strata; b'],
            [S[2].covariances; zeros(eltype(S[2].covariances), 1, size(S[2].covariances, 2), size(S[2].covariances, 3))]
        )
        push!(R.ages, convert(eltype(R.ages), 0))
    end

    # Send
    (R, (1 + S[1], R))
end

"""
    birthrate(P)

Stub function to be overloaded in implementation. Compute the exogenous
birth rate vector from the population P.
"""
function birthrate(P::population) end

"""
    hazardrate(P)

Stub function to be overloaded in implementation. Compute the exogenous
hazard rate matrix from the population P.
"""
function hazardrate(P::population) end

"""
    hazardrate(a)

Stub function to be overloaded in implementation. Compute the endogenous
hazard rate matrix a given age a.
"""
function hazardrate(a::Real) end

"""
    scatterrate(a)

Stub function to be overloaded in implementation. Compute the scattering
rate Eigen function.
"""
function scatterrate(a::Real) end

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
    conservesum!(A, a)

In place enforcement that the sum of the columns of A equals a. Note this 
assumes the inputs are well formed, there are no bounds or sanity checks.
"""
function conservesum!(A, a)
    d = a - sum(A)
    I = sortperm(A, rev=(d<0))
    A[I] .= A[I] .+ [sign(d) .* ones(Int64, abs(d)) ; zeros(Int64, length(A)-abs(d))]
    A
end

"""
    covariance(C, P, n)

Update of the covariances C from the transition probabilities P and the previous
state n. Each column of P is assumed to represent the exit  probabilities from a
single state and thus should add to 1 with all entries non-negative. Note this
assumes the inputs are well formed, there are no bounds or sanity checks.
"""
function covariance(C, P, n)
    U = (P.*n')*P'
    I = [1:stride(U, 2):length(U)...]
    U[I] .= (P*n) .- U[I]
    P*C*P' .+ U
end

end
