module Cliffords
using DocStringExtensions
using ..StabilizerFormalism
export Clifford, ⋊, CZ, CNOT, H, S, X, Y, Z, CZ!, CNOT!, H!, S!, Y!, Z!
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
(g::CZ!)(p::AbstractPauli) = p |> CNOT!(g.qubits...) |> H!(g.qubits[2:end]...)

# todo macro?
(g::H)(p::AbstractPauli) = H!(g.qubits...)(copy(p))
(g::S)(p::AbstractPauli) = S!(g.qubits...)(copy(p))
(g::X)(p::AbstractPauli) = X!(g.qubits...)(copy(p))
(g::Y)(p::AbstractPauli) = Y!(g.qubits...)(copy(p))
(g::Z)(p::AbstractPauli) = Z!(g.qubits...)(copy(p))
(g::CNOT)(p::AbstractPauli) = CNOT!(g.qubits...)(copy(p))
(g::CZ)(p::AbstractPauli) = CZ!(g.qubits...)(copy(p))


#= 
1Q Gates Remaining: H_XY, H_YZ, √X, √X†, √Y, √Y†
2Q Gates Remaining: CY, XCX, XCY, XCZ, YCX, YCY, YCZ, SWAP, ISWAP, ISWAP_DAG

=# 


### Legacy #### 
# TODO: Deprecate & Remove this insanity

# Hadamard flips the Z->X,X->Z, Y-> -Y
function ⋊(p::Pauli, op::H)
    signbit = xor(p.signbit, isequal.(zbits(p)[op.qubits], xbits(p)[op.qubits])...)
    q = Pauli(signbit, p.imagbit, bits(p))
    # flip the particular qubits
    q.bits[op.qubits, :] = q.bits[op.qubits, [2, 1]]
    q
end

function ⋊(x::Pauli, y::S)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    for i in y.qubits
        # X => Y;  Y => X bitflip on Z
        z.bits[i, 2] ⊻= z.bits[i, 1]
    end
    z
end




function ⋊(x::Pauli, y::X)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    for i in y.qubits
        # Z -> -Z # flips the sign if z bit is present
        z.signbit ⊻= x.bits[i, 2]
    end
    z
end


function ⋊(x::Pauli, y::Y)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    for i in y.qubits
        # X-> - X# flips the sign if x bit is present
        z.signbit ⊻= (x.bits[i, 1] ⊻ x.bits[i, 2])
    end
    z
end

function ⋊(x::Pauli, y::Z)
    z = Pauli(x.signbit, x.imagbit, copy(bits(x)))
    for i in [y.qubits...]
        # X -> -X # flips the sign if z bit is present
        z.signbit ⊻= x.bits[i, 1]
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
    end
    z
end

# Which one is better?  Composed Function or using Vector?
⋊(x::Clifford, y::Clifford) = x ∘ y
⋊(x::Clifford, y::ComposedFunction{<:Clifford}) = x ∘ y
⋊(x::Pauli, y::ComposedFunction{<:Clifford}) = (x ⋊ y.outer) ⋊ y.inner
⋊(f, g, h...) = ⋊(⋊(f, g), h...)
# ⋊(g::Clifford, h::Clifford) = [g, h]
# ⋊(g::Clifford, h::Vector{<:Clifford}) = [g, h...]
# ⋊(g::Vector{<:Clifford}, h::Clifford) = [g..., h]
# ⋊(f::Pauli, g, h...) = ⋊(⋊(f,g), h...)
# ⋊(g::Pauli, h::Vector{<:Clifford}) = ⋊(g,h...)

end
