"""
    vonFoersterHazards

Engine and utilities to forward propagate the von Foerster evolution
from the hazard rate, birth rate, and initial demographics.
"""
module vonFoersterHazards

import Base.iterate
export randomtruncate,
       conservesum!,
       covariance,
       evolve,
       hazardrate,
       birthrate

"""
    evolve(ages, population, count, size)

Iterable container for the population model. Conceptually the ages labels
the rows of the population matrix, and the columns are the states. Cohorts
are stored in reverse order, so that births can be pushed to the vector.
"""
struct evolve{R<:AbstractVector{Float64}, S<:AbstractVector{U} where U<:AbstractVector{Int64}, T<:AbstractVector{V} where V<:AbstractMatrix{Float64}}
    ages::R
    population::S
    covariances::T
    elapsed::Float64
    count::Int64
    size::Float64
end

function Base.iterate(E::evolve)
    (E.size > 0) ||
        throw(DomainError(size, "time increment size must be positive"))
    (issorted(E.ages, rev=true)) ||
        throw(DomainError(E.ages, "ages must be sorted in descending order"))
    (0.0 <= E.ages[end]) ||
        throw(DomainError(E.ages, "ages must be non-negative"))
    (0 <= minimum(minimum.(E.population))) ||
        throw(DomainError(E.population, "state occupancies in population must be non-negative"))
    (length(E.ages) == length(E.population)) ||
        throw(DimensionMismatch("cohorts in ages must equal cohorts in population"))
    (l, u) = extrema(length.(E.population))
    (l == u) ||
        throw(DimensionMismatch("states in the cohorts of population must be equal"))
    (length(birthrate(E.ages, E.population)) == u) ||
        throw(DimensionMismatch("states in the birth rate must equal the states in the cohorts of population"))
    H = hazardrate(E.ages, E.population)
    (size(H) == (l, u)) ||
        throw(DimensionMismatch("states in the extensize hazard rate must equal the states in the cohorts of population"))
    (size(hazardrate(E.ages[1], H)) == (l, u)) ||
        throw(DimensionMismatch("states in the intensize hazard rate must equal the states in the cohorts of population"))
    (E, 0)
end
function Base.iterate(E::evolve, step::Int64)
    if step > E.count
        nothing
    else
        
        # One time computation of the extensize birth rate
        b = birthrate(E.ages, E.population)
        
        # One time computation of the extensive hazard rate
        H = hazardrate(E.ages, E.population)
        
        # Curried transition function
        t(a, n) = conservesum!(randomtruncate.(exp(-E.size * hazardrate(a, H)) * n), sum(n))
        
        # Compute the transitions within each cohort from the intensive hazard rate
        E.population .= t.(E.ages, E.population)
        E.ages .= E.ages .+ E.size
        E.elapsed = E.elapsed + E.size
        
        # Youngest cohort is less than 1 year old, add births to youngest cohort
        if E.ages[end] < 1.0 then
            @inbounds E.population[end] = E.population[end] + b
            
        # Youngest cohort is more than 1 year old, generate a new youngest cohort
        else
            push!(E.ages, 0.0)
            push!(E.population, b)
        end
        (E, step + 1)
    end
end

"""
    birthrate(a, n)

Stub function to be overloaded in implementation. Compute the extensive
birth rate vector from the ages a and occupancy n of the population.
"""
function birthrate(a::AbstractVector{Float64}, n::AbstractVector{T})::T where T<:AbstractVector{Int64} end

"""
    hazardrate(a, n)

Stub function to be overloaded in implementation. Compute the extensive
hazard rate matrix from the ages a and occupancy n of the population.
"""
function hazardrate(a::AbstractVector{Float64}, n::AbstractVector{T} where T<:AbstractVector{Int64})::AbstractMatrix{Float64} end

"""
    hazardrate(a, H)

Stub function to be overloaded in implementation. Compute the intensive
hazard rate matrix from the extensive hazard rate matrix and a given age a.
"""
function hazardrate(a::Float64, H::T)::T where T<:AbstractMatrix{Float64} end

"""
    randomtruncate(x)

Randomly return the floor or the ceiling by comparing the fractional part of the
number to a sample from the uniform distribution on [0,1). If the fraction is
greater than the random sample return the ceiling otherwise return the floor.
"""
function randomtruncate(x::Float64)::Float64
    y = convert(Int64, trunc(x))
    y + convert(Int64, rand() < (x-y))
end

"""
    conservesum!(A, a)

In place enforcement that the sum of the columns of A equals the a. Note this 
assumes the inputs are well formed, there are no bounds or sanity checks.
"""
function conservesum!(A::T, a::Int64)::T where T<:AbstractVector{Int64}
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
function covariance(C::T, P::T, n::AbstractVector{Int64})::T where T<:AbstractMatrix{Float64}
    U = P*(P.*n)'
    i = [1:(1 + size(U, 1)):length(U)...]
    U[i] .= (P*n) .- U[i]
    P*C*P' + U
end

end
