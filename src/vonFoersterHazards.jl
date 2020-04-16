module vonFoersterHazards

import Base.iterate
export randomtruncate,
       conservesum!,
       covariance!,
       evolve,
       hazardrate,
       probability,
       transition!,
       birth!

"""
    evolve(boundary, population, ages, count, size)

Iterable container for the population model.
"""
struct evolve
    boundary::AbstractVector{Int64}
    population::AbstractVector{AbstractVector{Int64}}
    ages::AbstractVector{Float64}
    count::Int64
    size::Float64
end
Base.iterate(E::evolve) = ((E.population, E.ages), 0)
function Base.iterate(E::evolve, step::Int64)
    if step > E.count
        nothing
    else
        ((a, n) -> transition!(a, n, E)).(E.ages, E.population)
        birth!(E)
        ((E.population, E.ages), step + 1)
    end
end

"""
    transition!(a, n, E)

Single cohort transition update from the hazard rate matrix. Meant to be
called through a vectorized broadcast across the whole population.
"""
function transition!(a::Float64, n::AbstractVector{Int64}), E::evolve)
    conservesum!(randomtruncate.(hazardrate(a, E) * n), sum(n))
end

"""
    birth!(E)

If the youngest cohort is less than 1 year old, add births to youngest cohort.
Otherwise youngest cohort is more than 1 year old, generate a new youngest cohort.
"""
function birth!(E::evolve)
        if E.ages[end] < 1 then
            E.population[end] = E.population[end] + E.boundary
        else
            push!(E.ages, 0)
            push!(E.population, E.boundary)
        end
        E
end

"""
    hazardrate(a, n)

Stub function to be overloaded in implementation. Compute the hazard rate
matrix for a given cohort state occupancy n and age a.
"""
function hazardrate(a, n) end

"""
    probability(H, t)

Compute the transition probabilities from the hazard rate matrix H, given a small
time step t, assuming the diagonal is non-negative exit rates, off diagonal are
non-positive entrance rates, and the columns represent the exit rates from a
single state and thus should add to 0. Note this assumes the inputs are well formed,
there are no bounds or sanity checks.
"""
function probability(H::AbstractMatrix{Float64}, t::Float64)
    exp(-t * H)
end

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
function conservesum!(A::AbstractMatrix{Int64}, a::AbstractVector{Int64})
    l = length(a)
    d = a .- sum(A, dims=1)
    for i in eachindex(a)
        I = sortperm(A[:, i], rev=(d[i]<0))
        A[I, i] .= A[I, i] .+ [sign(d[i]) .* ones(Int64, abs(d[i])) ; zeros(Int64, length(l)-abs(d[i]))]
    end
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
