import Base.Matrix

"""
Creates the symplectic boolean representation of the PauliOperator
"""
function symplectic(s::Pauli)
    v = vec(s)
    r = [
        [e in [σx, σy] ? 1 : 0 for e in v]...,
        [e in [σz, σy] ? 1 : 0 for e in v]...
    ]
    BitVector(r)
end

checkmatrix(s::Vector{<:Pauli}) = vcat(symplectic.(s)'...)

function Matrix(x::PauliPrimitive)
    v = (
        [1 0; 0 1],
        [0 1; 1 0],
        im*[0 -1; 1 0],
        [1 0; -1 0],
    )
    v[Integer(x) + 1]
end

Matrix(x::Pauli) = x.coeff * kron(Matrix.(x.val)...)