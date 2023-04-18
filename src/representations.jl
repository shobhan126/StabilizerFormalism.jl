import Base.Matrix

"""
    symplectic(s::Pauli)
Creates the symplectic 2N boolean representation of the PauliOperator. [X...|Z...]
`X...` bitvector is 1 if X operator at that site, 0 if Z operator at the site

# Examples
```julia-repl
julia> symplectic(p"IXY")
[0 1 1 0 0 1]
julia> length(symplectic(p"IXYZ))
6
"""
symplectic(p::AbstractPauli) = [xbits(p); zbits(p)]

"""
Given a vector of Pauli Operators with same number of qubits, return a matrix rows corresponding to the 
symplectic representation of each Pauli{N}.
"""
checkmatrix(s::Vector{<:AbstractPauli{N}}) where N= vcat(symplectic.(s)'...)

"""
    Matrix(x::Pauli)

Returns matrix representation of the 1 qubit PauliPrimitive
"""
function Matrix(p::Pauli, sparse=true)
    # using a tuple has less allocations than using a vector!

    # c  sign * im 
    c = (-1^p.signbit  * 1im^p.imagbit)
    if sparse
        f = b -> (
            SparseMatrixCSC([1 0; 0 1]), 
            SparseMatrixCSC([0 1; 1 0]),
            SparseMatrixCSC([0 -im; im 0]), 
            SparseMatrixCSC([1 0; 0 -1])
        )[sum(b)+1] * c
    else
        f = b -> (
            [1 0; 0 1], 
            [0 1; 1 0],
            [0 -im; im 0], 
            [1 0; 0 -1]
        )[sum(b)+1] * c
    end
    mapreduce(f, kron, eachrow(bits(p))) * c
        
end