module Cliffords

    using ..StabilizerFormalism
    using SparseArrays: SparseVector, dropzeros!
    export Clifford, ⋊, CNOT, H, S, X, Y, Z
    using StabilizerFormalism.Paulis: bits, zbits, xbits
    # I don't think i need this kind of subtyping at the moment? 
    abstract type Gate end
    abstract type Unitary <: Gate end
    abstract type Clifford <: Unitary end

    macro cliffordgenerators(args...)
        for arg in args
            expr = quote 
                    struct ($arg) <: Clifford
                    qubits::Vector{<:Int}
                    ($arg)(qubits...) = new([qubits...])
                    end
            end
            eval(expr)
        end
    end

    @cliffordgenerators CNOT H S X Y Z

    # \to operator

    # Hadamard flips the Z->X,X->Z, Y-> -Y
    function ⋊(p::Pauli, op::H)
        signbit = xor(p.signbit, isequal.(zbits(p)[op.qubits], xbits(p)[op.qubits])...)
        q = Pauli(signbit, p.imagbit, bits(p))
        # flip the particular qubits
        q.bits[op.qubits, :] = q.bits[op.qubits, [2,1]]
        dropzeros!(q.bits)
        q
    end


    function ⋊(x::Pauli, y::S)
        z = Pauli(x.signbit, x.imagbit, deepcopy(x.bits))
        for i in y.qubits
            # X => Y;  Y => X bitflip on Z
            z.bits[i,2] ⊻= z.bits[i,1]
            dropzeros!(z.bits[i,:])
        end
        z
    end


    function ⋊(x::Pauli, y::X)
        z = Pauli(x.signbit, x.imagbit, deepcopy(x.bits))
        for i in y.qubits
            # Z -> -Z # flips the sign if z bit is present
            z.signbit ⊻= x.bits[i,2] 
            dropzeros!(z.bits[i,:])
        end
        z
    end


    function ⋊(x::Pauli, y::Y)
        z = Pauli(x.signbit, x.imagbit, deepcopy(x.bits))
        for i in y.qubits
            # X-> - X# flips the sign if x bit is present
            z.signbit ⊻= (x.bits[i,1] ⊻ x.bits[i,2])
            dropzeros!(z.bits[i,:])
        end
        z
    end

    function ⋊(x::Pauli, y::Z)
        z = Pauli(x.signbit, x.imagbit, deepcopy(x.bits))
        for i in [y.qubits...]
            # X -> -X # flips the sign if z bit is present
            z.signbit ⊻= x.bits[i,1]
            dropzeros!(z.bits[i,:])
        end
        z
    end

    function ⋊(x::Pauli, y::CNOT)
        z = Pauli(x.signbit, x.imagbit, deepcopy(x.bits))
        control = y.qubits[1]
        for i in y.qubits[2:end]
            # X1 -> X1X2; Z2 => Z1Z2
            z.bits[i,1] ⊻= x.bits[control,1]
            z.bits[control, 2] ⊻= x.bits[i, 2]
            dropzeros!(z.bits[i,:])
        end
        dropzeros!(z.bits[control,:])
        z
    end

    # Which one is better?  Composed Function or using Vector?
    ⋊(x::Clifford, y::Clifford) = x ∘ y
    ⋊(x::Clifford, y::ComposedFunction{<:Clifford}) = x ∘ y
    ⋊(x::Pauli, y::ComposedFunction{<:Clifford}) = (x ⋊ y.outer) ⋊ y.inner
    ⋊(f, g, h...) = ⋊(⋊(f, g), h...)

end


# ⋊(g::Clifford, h::Clifford) = [g, h]
# ⋊(g::Clifford, h::Vector{<:Clifford}) = [g, h...]
# ⋊(g::Vector{<:Clifford}, h::Clifford) = [g..., h]
# ⋊(f::Pauli, g, h...) = ⋊(⋊(f,g), h...)
# ⋊(g::Pauli, h::Vector{<:Clifford}) = ⋊(g,h...)
