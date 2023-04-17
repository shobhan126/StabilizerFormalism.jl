# TODO Autogenerate levi-civita dictionary indexed by tuples and keep it in the namespace.
ε = [
    0 1 -1;
    -1 0 1;
    1 -1 0
] * im

global const paulibits = Dict(
    (false, false) => 0,
    (true, false) => 1,
    (true, true) => 2,
    (false, true) => 3,
)

"""
Levi Civita AntiSymmetric Function 
    $(FUNCTIONNAME)(i,j)
"""
levicivita(i::AbstractArray, j::AbstractArray) =  iszero(i) | iszero(j) | (i==j) ? 1 : ε[paulibits[i...], paulibits[j...]]
