module utilities

"""
    randomtruncate(x)

Randomly return the floor or the ceiling by comparing the fractional part of the
number to a sample from the uniform distribution on [0,1). If the fraction is
greater than the random sample return the ceiling otherwise return the floor.
"""
function randomtruncate(x)
    y = convert(Int64, trunc(x))
    y + ((x-y) < rand() ? 0 : 1)
end

"""
    conservesum(A, a)

For each dimension in A matching a dimension in a ensure the sum over the unmatched
dimensions equals the element at a.
"""
function conservesum(A, a)
    b = sum(A, (1+ndims(a)):ndims(A))
end

end