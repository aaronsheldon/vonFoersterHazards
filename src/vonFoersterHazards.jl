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
       covariance!,
       cohort,
       population,
       evolve,
       hazardrate,
       birthrate

"""
    strata(elapsed, age, occupancy, covariance)

Return type for indexing into the population. Container for the state occupancies
of a specific age group, and the covariances between the occupancies.
"""
struct strata{
        Q<:Real,
        R<:Real,
        S<:AbstractVector{U} where U<:Real,
        T<:AbstractMatrix{V} where V<:Real
    }
    elapsed::Q
    age::R
    cohort::S
    covariance::T
end
strata(a, o) = strata(0.0, a, o, zeros(Float64, size(o, 1), size(o, 1)))

"""
    population(elapsed, ages, cohorts, covariances)

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
    cohorts::S
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
Base.getindex(P::population, i) = cohort(P.elapsed, P.ages[i], view(P.cohorts, i, :), view(P.covariances, i, :, :))
function Base.setindex!(P::population, C::cohort, i)
    P.ages[i] = C.age
    P.cohorts[i, :] .= C.cohort[:]
    P.covariances[i, :, :] .= C.covariance[:, :]
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
    (E.size > 0) ||
        throw(DomainError(size, "time increment must be positive."))
    
    (issorted(E.initial.ages, rev=true)) ||
        throw(DomainError(E.initial.ages, "ages must be sorted in descending order."))
    
    (convert(eltype(E.initial.ages), 0) <= E.initial.ages[end]) ||
        throw(DomainError(E.initial.ages, "ages must be non-negative."))
    
    (convert(eltype(E.initial.cohorts), 0) <= minimum(E.initial.cohorts)) ||
        throw(DomainError(E.initial.cohorts, "cohorts must be non-negative."))
    
    (minimum(E.initial.covariances) == convert(eltype(E.initial.covariances), 0) == maximum(E.initial.covariances)) ||
        throw(DomainError(E.initial.covariances, "covariances must be zero."))   

    H = hazardrate(E.initial)
    (al,) = size(E.initial.ages)
    (bl,) = size(birthrate(E.initial))
    (cl, cw) = size(E.initial.cohorts)
    (dl, dw) = size(H)
    (el, ew) = size(hazardrate(E.initial.ages[1], H))
    (fl, fw, fh) = size(E.initial.covariances) 
    
    (al == cl == fl) ||
        throw(DimensionMismatch("stratas in ages, cohorts, and covariances must be equal.")) 

    (bl  == cw == dl == dw == el == ew == fw == fh) ||
        throw(DimensionMismatch("states in the covariances, birth rate, extensize hazard rate, intensize hazard rate, and cohorts must be equal."))
    
    (E.initial, (0, E.initial))
end
function Base.iterate(E::evolve, S)
    if S[1] > E.count
        nothing
    else
        
        # One time computation of the extensize birth rate
        b = birthrate(S[2])
        
        # One time computation of the extensive hazard rate
        H = hazardrate(S[2])
        
        # Curried transition function
        function t(c)
            P = exp(-E.size * hazardrate(a, H))
            r = cohort(
                c.elapsed + E.size,
                c.age + E.size,
                conservesum!(randomtruncate.(P * c.cohort), sum(c.cohort)),
                c.covariance
            )
            covariance!(r.covariance, P, r.cohort)
            r
        end
        
        # Compute the transitions within each cohort
        S[2] .= t.(S[2])
        
        # Youngest cohort is less than 1 year old, add births to youngest cohort
        if S[1].ages[end] < 1.0 then
            R = population(
                S[2] + E.size,
                S[2].ages,
                S[2].cohorts,
                S[2].covariances
            )
            R.cohorts[end, :] .= R.cohorts[end, :] .+ b'
            
        # Youngest cohort is more than 1 year old, generate a new youngest cohort
        else
            R = population(
                S[2] + E.size,
                S[2].ages,
                [S[2].cohorts; b'],
                S[2].covariances
            )
            push!(R.ages, convert(eltype(R.ages), 0))
        end
        
        # Send
        (R, (1 + S[1], R))
    end
end

"""
    birthrate(P)

Stub function to be overloaded in implementation. Compute the extensive
birth rate vector from the population P.
"""
function birthrate(P::population) end

"""
    hazardrate(P)

Stub function to be overloaded in implementation. Compute the extensive
hazard rate matrix from the population P.
"""
function hazardrate(P::population) end

"""
    hazardrate(a, H)

Stub function to be overloaded in implementation. Compute the intensive
hazard rate matrix from the extensive hazard rate matrix H and a given age a.
"""
function hazardrate(a, H) end

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
    covariance!(C, P, n)

Update of the covariances C from the transition probabilities P and the previous
state n. Each column of P is assumed to represent the exit  probabilities from a
single state and thus should add to 1 with all entries non-negative. Note this
assumes the inputs are well formed, there are no bounds or sanity checks.
"""
function covariance!(C, P, n)
    U = (P.*n')*P'
    i = [1:(1 + size(U, 1)):length(U)...]
    U[i] .= (P*n) .- U[i]
    C .= P*C*P' .+ U
    C
end

end
