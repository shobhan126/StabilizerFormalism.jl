"""
Module defining Pauli Operators
    $(EXPORTS)
"""
module Paulis
# Pauli Primitive Enum
import Base: show, *, -, ==, abs, vec
import Base.Iterators.product

using DocStringExtensions
import LinearAlgebra: mul!, lmul!, rmul!, kron
using LoopVectorization: @inbounds

include("paulitable.jl")

#### Common Interface for different Representations #####
abstract type AbstractPauli{N} end

### Accessors 
function xbits(p::AbstractPauli) end
function zbits(p::AbstractPauli) end
function bits(p::AbstractPauli) end
function signbit(p::AbstractPauli) end
function imagbit(p::AbstractPauli) end

### Math Operations ####
function Base.adjoint(p::AbstractPauli) end
function Base.:-(x::AbstractPauli) end
function Base.:*(p::AbstractPauli, q::AbstractPauli) end
function Base.:*(y::AbstractPauli, x::Number) end
function Base.:*(x::Number, y::AbstractPauli) end
function abs(p::AbstractPauli) end
function commuting(x::AbstractPauli, y::AbstractPauli) end
function anticommuting(x::AbstractPauli, y::AbstractPauli) end
function ==(x::AbstractPauli, y::AbstractPauli) end

#### Inplace Operations
function lmul!(a::Number, B::AbstractPauli) end
function rmul!(A::AbstractPauli, b::Number) end
function lmul!(A::AbstractPauli, B::AbstractPauli) end
function rmul!(A::AbstractPauli, B::AbstractPauli) end
function mul!(C::AbstractPauli, A::AbstractPauli, B::AbstractPauli, α::Number, β::Number) end
function mul!(Y::AbstractPauli, A::AbstractPauli, B::AbstractPauli) end


# Common Operations

# REPL representation for the Pauli Operator
function Base.show(io::IO, p::AbstractPauli)
    prefix = signbit(p) ? "-" : ""
    prefix *= imagbit(p) ? "im" : ""
    print(io, prefix)
    s = ["I", "X", "Z", "Y"]
    foreach(x -> print(io, s[x+1]), xbits(p) + 2 * zbits(p))
end

isneg(p::AbstractPauli) = signbit(p)


# """
# Create an N-Qubit Pauli Operator
#     $(TYPEDEF)
# with fields:
#     $(FIELDS)

# # Examples
# ```julia-repl
# julia> Pauli(:XIZI)
# XIZI
# ```
# """
mutable struct Pauli{N} <: AbstractPauli{N}
    """
        true if coefficient has `-` sign
    """
    signbit::Bool
    """
        true if coefficient is imaginary
    """
    imagbit::Bool
    """
        Sparse X-Bit Representation of the Pauli Operator.
    """
    xbits::BitArray{1}
    zbits::BitArray{1}
end


"""
Create a Pauli Operator from its symplectic representation [xbits; zbits]

    $(FUNCTIONNAME)(xibts, zbits, neg=false, imag=false)

"""
Pauli(xbits::AbstractVector, zbits::AbstractVector, signbit=false, imagbit=false) = Pauli{length(xbits)}(signbit, imagbit, BitArray(xbits), BitArray(zbits))


"""
Create a Pauli Operator from its symplectic representation [xbits; zbits]
    
    $(FUNCTIONNAME)(xz, neg=false, imag=false)
"""
Pauli(xz::AbstractVector, neg=false, imag=false) = Pauli(xz[1:Int(length(x) / 2)], xz[Int(length(x) / 2)+1:end], neg, imag)


"""
Macro to instantiate a Pauli{N} type from a string
    @p_str("IXYZ")
    p"IXYZ..."
"""
macro p_str(p)
    conv = Dict(
        'x' => [true, false],
        'y' => [true, true],
        'z' => [false, true],
        'i' => [false, false]
    )
    m = BitArray(hcat([conv[lowercase(v)] for v in p]...)')
    return Pauli(m[:, 1], m[:, 2])
end

"""
@pauli XYZXII...

Convenience macro for constructing a Pauli{N} type.
"""
macro pauli(ex)
    s = string(ex)
    :(@p_str $s)
end

"""
Construct a Puali type by passing a symbol.
    
    $(FUNCTIONNAME)(sym)

"""
function Pauli(sym::Symbol)
    eval(:(@pauli $sym))
end

"""
Construct a Pauli Operator type by passing a string representation
    
    $(FUNCTIONNAME)(s)
    
"""
Pauli(s::AbstractString) = eval(:(@p_str $s))


# Accessors
xbits(p::Pauli) = p.xbits
zbits(p::Pauli) = p.zbits
bits(p::Pauli) = hcat(xbits(p), zbits(p))
signbit(p::Pauli) = p.signbit
imagbit(p::Pauli) = p.imagbit

# In-Place Operations
function mul!(C::Pauli{N}, A::Pauli{N}, B::Pauli{N}) where {N}
    C.signbit = A.signbit ⊻ B.signbit
    C.imagbit = A.imagbit ⊻ B.imagbit
    @inbounds for x in 1:N
        C.signbit ⊻= signlogic(A.xbits[x], A.zbits[x], B.xbits[x], B.zbits[x])
        C.imagbit ⊻= imaglogic(A.xbits[x], A.zbits[x], B.xbits[x], B.zbits[x])
        C.xbits[x] = A.xbits[x] ⊻ B.xbits[x]
        C.zbits[x] = A.zbits[x] ⊻ B.zbits[x]
    end
end

mul!(A::Pauli{N}, B::Pauli{N}) where {N} = mul!(A, A, B)

function Base.:*(p::Pauli{N}, q::Pauli{N}) where {N}
    o = Pauli{N}(false, false, similar(xbits(p)), similar(zbits(p)))
    mul!(o, p, q)
    o
end


"""
Take a adjoint of the pauli operator `iXYZ...`  => `-iXYZ...`.

    $(FUNCTIONNAME)(p::Pauli)
"""
Base.adjoint(p::Pauli) = p.imagbit ? Pauli(p.signbit ⊻ true, p.imagbit) : p



# Scalar Multiplication
*(y::Pauli, x::Number) = Pauli((x < 0) ⊻ y.signbit, y.imagbit, y.bits)
*(y::Pauli, x::Complex) = isequal(x, im) ?
                          Pauli(y.signbit ⊻ y.imagbit, ~y.imagbit, y.bits) : isequal(x, -im) ?
                          Pauli(~y.signbit ⊻ y.imagbit, ~y.imagbit, y.bits) : throw(MethodError)
*(x::Number, y::Pauli) = *(y, x)


#### Negative sign -Pauli(1, [...]) = Pauli(-1, [...])
Base.:-(x::Pauli) = Pauli(x.signbit ⊻ true, x.imagbit, x.bits)

"""
Return a new Pauli with the coefficient of `p` set to 1
    abs(p::Pauli)
"""
abs(p::Pauli) = Pauli(false, false, p.bits)

==(x::Pauli, y::Pauli) = (x.signbit == y.signbit) & isequal(bits(x), bits(y)) & (x.imagbit == y.imagbit)


# todo: is it better to have a an empty struct or just pass around types
# abstract type PauliPrimitive end
# for x in "IXYZ"
#     s = Symbol("Pauli$(x)")
#     constname = Symbol(lowercase(x))
#     consttype = QuoteNode(s)
#     :(struct $s <: PauliPrimitive end) |> eval
#     :($constname = eval(Expr(:call, $consttype))) |> eval
# end
# Also, might help if the internals is binary? 

"""
Alias for kronecker product 
    
    ⊗ === kron
"""
const ⊗ = kron
function Base.kron(p1::Pauli, p2::Pauli)
    Pauli(
        (p1.signbit ⊻ p2.signbit ⊻ (p1.imagbit & p2.imagbit)),
        (p1.imagbit ⊻ p2.imagbit),
        vcat(p1.bits, p2.bits)
    )
end

"""
Return true if two paulis commute
    $(FUNCTIONNAME)(x,y)

"""
commuting(x::AbstractPauli, y::AbstractPauli) = imagbit(x * y) & (imagbit(x) ⊻ imagbit(y))

"""
Return true if two paulis anticommute
    $(FUNCTIONNAME)(x,y)
"""
anticommuting(x::AbstractPauli, y::AbstractPauli) = ~commuting(x, y)

# setting up iteration interface
Base.length(x::Pauli) = 1
Base.iterate(x::Pauli) = (x, nothing)
Base.iterate(x::Pauli, n::Nothing) = nothing

include("representations.jl")

export Pauli, @p_str, @pauli
export bits, xbits, zbits, symplectic, checkmatrix
export ⊗, commuting, anticommuting
end