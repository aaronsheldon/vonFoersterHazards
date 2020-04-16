module vonFoersterHazards

import Base.iterate
export randomtruncate,
       conservesum!,
       covariance!,
       evolve,
       hazardrate

"""
    evolve(boundary, population, ages, count, size)

Iterable container for the population model. Conceptually the boundary
labels the columns and the ages labels the rows of the population matrix.
"""
struct evolve
    boundary::AbstractVector{Int64}
    ages::AbstractVector{Float64}
    population::AbstractVector{AbstractVector{Int64}}
    count::Int64
    size::Float64
end
Base.iterate(E::evolve) = ((E.ages, E.population), 0)
function Base.iterate(E::evolve, step::Int64)
    if step > E.count
        nothing
    else
        
        # One time computation of the extensive hazard rate
        H = hazardrate(E.ages, E.population)
        
        # Compute the transitions across the cohorts from the intensive hazard rate
        for i in eachindex(E.ages)
            E.population[i] = conservesum!(randomtruncate.(exp(-t * hazardrate(E.ages[i], H)) * E.population[i]), sum(E.population[i]))
        end

        # Youngest cohort is less than 1 year old, add births to youngest cohort
        if E.ages[end] < 1 then
            E.population[end] = E.population[end] + E.boundary
            
        # Youngest cohort is more than 1 year old, generate a new youngest cohort
        else
            push!(E.ages, 0)
            push!(E.population, E.boundary)
        end
        ((E.ages, E.population), step + 1)
    end
end

"""
    hazardrate(a, H)

Stub function to be overloaded in implementation. Compute the intensive
hazard rate matrix from the extensive hazard rate matrix and a given age a.
"""
function hazardrate(a::Float64, H::AbstractMatrix{Float64}) end

"""
    hazardrate(a, n)

Stub function to be overloaded in implementation. Compute the extensive
hazard rate matrix from ages a and population occupancies n.
"""
function hazardrate(a::AbstractVector{Float64}, n::AbstractVector{Float64}) end

"""
    randomtruncate(x)

Randomly return the floor or the ceiling by comparing the fractional part of the
number to a sample from the uniform distribution on [0,1). If the fraction is
greater than the random sample return the ceiling otherwise return the floor.
"""
function randomtruncate(x::Float64)
    y = convert(Int64, trunc(x))
    y + convert(Int64, rand() < (x-y))
end

"""
    conservesum!(A, a)

In place enforcement that the sum of the columns of A equals the corresponding
element of a. Note this assumes the inputs are well formed, there are no bounds
or sanity checks.
"""
function conservesum!(A::AbstractVector{Int64}, a::Int64)
    d = a - sum(A)
    I = sortperm(A, rev=(d<0))
    A[I] .= A[I] .+ [sign(d) .* ones(Int64, abs(d)) ; zeros(Int64, length(A)-abs(d))]
    A
end

"""
    covariance!(C, P, n)

In place update of the covariances C from the transition probabilities P and
the previous state n. Each column of P is assumed to represent the exit 
probabilities from a single state and thus should add to 1 with all entries
non-negative. Note this assumes the inputs are well formed, there are no bounds
or sanity checks.
"""
function covariance!(C::AbstractMatrix{Float64}, P::AbstractMatrix{Float64}, n::AbstractVector{Int64})
    C
end

end # module
