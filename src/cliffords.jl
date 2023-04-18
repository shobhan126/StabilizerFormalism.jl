module Cliffords

using ..StabilizerFormalism
export Clifford, ⋊, CNOT, H, S, X, Y, Z, H!, S!, Y!, Z!
using ..StabilizerFormalism.Paulis: bits, zbits, xbits, Pauli
# I don't think i need this kind of subtyping at the moment? 
abstract type Gate <: Function end
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

@cliffordgenerators CZ CNOT H S X Y Z

@cliffordgenerators CZ! CNOT! H! S! X! Y! Z!



"""
In-place Hadamard gate conjugation
   
    $(SIGNATURES)

Passing Multiple arguments will be read as concurrent / transveral operation.
"""
function (g::H!)(p::AbstractPauli)
    # Hadamard flips the Z->X,X->Z, Y-> -Y
    @inbounds for i ∈ g.qubits
        p.signbit ⊻= (zbits(p)[i] & xbits(p)[i])
        zbits(p)[i], xbits(p)[i] = xbits(p)[i], zbits(p)[i]
    end
    p
end

"""
In-place S (Phase Gate) gate conjugation
    
    $(SIGNATURES)

Passing Multiple arguments will be read as concurrent / transveral operation.

"""
function (g::S!)(p::Pauli)
    # X => Y;  Y => X bitflip on Z
    @inbounds for i ∈ g.qubits
        zbits(p)[i] ⊻= xbits(p)[i, 1]
    end
    p
end


"""
In-place X gate conjugation
    
    $(SIGNATURES)

Passing Multiple arguments will be read as concurrent / transveral operation.

"""
function (g::X!)(p::AbstractPauli)
    # Z -> -Z # flips the sign if z bit is present
    @inbounds for i ∈ g.qubits
        p.signbit ⊻= (zbits(p))
    end
    p
end


"""
In-place Y gate conjugation
    
    $(SIGNATURES)

Passing Multiple arguments will be read as concurrent / transveral operation.

"""
function (g::Y!)(p::AbstractPauli)
    # X-> - X# flips the sign if x bit is present
    for i ∈ g.qubits
        p.signbit ⊻= xbits(p)[i]
    end
    p
end


"""
In-place Z gate conjugation
    
    $(SIGNATURES)

Passing Multiple arguments will be read as concurrent / transveral operation.

"""
function (g::Z!)(p::AbstractPauli)
    # X -> -X # flips the sign if z bit is present
    @inbounds for i ∈ g.qubits
        p.signbit ⊻= zbits(p)[i]
    end
    p
end


"""
In-place CNOT gate conjugation.

    $(SIGNATURES)

The first qubit is the control qubit, subsequent qubits are target qubits.
Passing Multiple arguments interpretted as C-ⁿNOT.
"""
function (g::CNOT!)(p::AbstractPauli)
    @inbounds for i ∈ g.qubits[2:end]
        # X1 -> X1X2; Z2 => Z1Z2
        p.xbits[i] ⊻= p.xbits[g.qubits[1]]
        p.zbits[g.qubits[1]] ⊻= p.zbits[i]

    end
    p
end

"""
In-place CZ gate conjugation.

    $(SIGNATURES)

The first qubit is the control qubit, subsequent qubits are target qubits.
Passing Multiple arguments interpretted as C-ⁿZ.
"""
(g::CZ!)(p::AbstractPauli) = CNOT!(g.qubits...) |> H!(g.qubits[2:end]...)



(g::H)(p::AbstractPauli) = (g!)(copy(p))
(g::S)(p::AbstractPauli) = (g!)(copy(p))
(g::X)(p::AbstractPauli) = (g!)(copy(p))
(g::Y)(p::AbstractPauli) = (g!)(copy(p))
(g::Z)(p::AbstractPauli) = (g!)(copy(p))
(g::CNOT)(p::AbstractPauli) = (g!)(copy(p))
(g::CZ)(p::AbstractPauli) = (g!)(copy(p))




### Legacy #### 


# Hadamard flips the Z->X,X->Z, Y-> -Y
function ⋊(p::Pauli, op::H)
    signbit = xor(p.signbit, isequal.(zbits(p)[op.qubits], xbits(p)[op.qubits])...)
    q = Pauli(signbit, p.imagbit, bits(p))
    # flip the particular qubits
    q.bits[op.qubits, :] = q.bits[op.qubits, [2, 1]]
    dropzeros!(q.bits)
    q
end

function ⋊(x::Pauli, y::S)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    for i in y.qubits
        # X => Y;  Y => X bitflip on Z
        z.bits[i, 2] ⊻= z.bits[i, 1]
        dropzeros!(z.bits[i, :])
    end
    z
end




function ⋊(x::Pauli, y::X)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    for i in y.qubits
        # Z -> -Z # flips the sign if z bit is present
        z.signbit ⊻= x.bits[i, 2]
        dropzeros!(z.bits[i, :])
    end
    z
end


function ⋊(x::Pauli, y::Y)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    for i in y.qubits
        # X-> - X# flips the sign if x bit is present
        z.signbit ⊻= (x.bits[i, 1] ⊻ x.bits[i, 2])
        dropzeros!(z.bits[i, :])
    end
    z
end

function ⋊(x::Pauli, y::Z)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    for i in [y.qubits...]
        # X -> -X # flips the sign if z bit is present
        z.signbit ⊻= x.bits[i, 1]
        dropzeros!(z.bits[i, :])
    end
    z
end


function ⋊(x::Pauli, y::CNOT)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    control = y.qubits[1]
    for i in y.qubits[2:end]
        # X1 -> X1X2; Z2 => Z1Z2
        z.bits[i, 1] ⊻= x.bits[control, 1]
        z.bits[control, 2] ⊻= x.bits[i, 2]
        dropzeros!(z.bits[i, :])
    end
    dropzeros!(z.bits[control, :])
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