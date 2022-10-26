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
symplectic(p::Pauli) = vec(p.bits)

"""
Given a vector of Pauli Operators with same number of qubits, return a matrix rows corresponding to the 
symplectic representation of each Pauli{N}.
"""
checkmatrix(s::Vector{<:Pauli}) = vcat(symplectic.(s)'...)

"""
    Matrix(x::PauliPrimitive)

Returns matrix representation of the 1 qubit PauliPrimitive
"""
# function Matrix(x::PauliPrimitive)
#     v = (
#         [1 0; 0 1],
#         [0 1; 1 0],
#         im*[0 -1; 1 0],
#         [1 0; -1 0],
#     )
#     v[Integer(x) + 1]
# end

"""
    Matrix(x::PauliPrimitive)

Returns a tensor product matrix of the Pauli{N} operator.
"""
# Matrix(x::Pauli) = x.coeff * kron(Matrix.(x.val)...)