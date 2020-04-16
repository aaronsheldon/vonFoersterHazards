module vonFoersterHazards

export randomtruncate, conservesum!, updatecovariance!

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
function conservesum!(A::Vector{Int64}, a::Int64)
    d = a - sum(A)
    I = sortperm(A, rev=(d<0))
    A[I] .= A[I] .+ [sign(d) .* ones(Int64, abs(d)) ; zeros(Int64, length(A)-abs(d))]
    A
end
function conservesum!(A::Matrix{Int64}, a::Vector{Int64})
    l = length(a)
    d = a .- sum(A, dims=1)
    for i in eachindex(a)
        I = sortperm(A[:, i], rev=(d[i]<0))
        A[I, i] .= A[I, i] .+ [sign(d[i]) .* ones(Int64, abs(d[i])) ; zeros(Int64, length(l)-abs(d[i]))]
    end
    A
end

"""
    updatecovariance!(C, H, n)

In place update of the covariances from the hazard rates and the previous state.
Note this assumes the inputs are well formed, there are no bounds or sanity checks.
"""
function updatecovariance!(C, H, n)
    0
end

end # module
