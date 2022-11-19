# Pauli Primitive Enum
import Base: show, *, -, ==, abs, vec
import Base.Iterators.product
using SparseArrays

include("paulitable.jl")

"""
N-Qubit Pauli Group
    Pauli(:XYZI)
Creates a Pauli{N} object where N is the number of qubits.

# Examples
```julia-repl
julia> Pauli(:XIZI)
XIZI
```
"""
struct Pauli
    signbit::Bool
    imagbit::Bool
    bits::SparseMatrixCSC
end

isneg(p::Pauli) = p.signbit

function Pauli(x::AbstractVector,z::AbstractVector, neg=false, imag=false)
    Pauli(neg, imag, SparseMatrixCSC(BitMatrix(hcat(x,z))))
end

Pauli(x::AbstractVector, neg=false, imag=false) = Pauli(x[1:Int(length(x)/2)], x[Int(length(x)/2)+1:end], neg, imag)

bits(p::Pauli) = p.bits
xbits(p::Pauli) = p.bits[:,1]
zbits(p::Pauli) = p.bits[:,2]
Base.adjoint(p::Pauli) = p.imagbit ? p.signbit ⊻ true : p


# REPL representation for the Pauli Operator
function Base.show(io::IO, p::Pauli)
    prefix = p.signbit ? "-" : ""
    prefix *= p.imagbit ? "im" : ""
    print(io, prefix)
    s = ["I", "X", "Z",  "Y"]
    foreach(x-> print(io, s[x+1]), xbits(p) + 2*zbits(p))    
end

function Base.:*(p::Pauli, q::Pauli)
    newbits = bits(p) .⊻ bits(q)
    coeff = (im)^(p.imagbit+q.imagbit) * prod([ϵ(i,j) for (i, j) in zip(eachrow(p.bits), eachrow(q.bits))])
    Pauli((coeff.re + coeff.im < 0) ⊻ p.signbit ⊻ q.signbit, ~isreal(coeff), newbits)
end


"""
    Pauli(s::Symbol)

Construct a Puali{N} type by passing a symbol
"""
function Pauli(s::Symbol)
    eval(:(@pauli $s))
end

"""
    Pauli(s::AbstractString)
    
Construct a Pauli{N} type by passing a string representation
"""
function Pauli(s::AbstractString)
    eval(:(@p_str $s))
end

# function Pauli(coeff::Number, p::PauliPrimitive)
#     @assert coeff in [1, -1, im, -im]
#     @assert p isa PauliPrimitive
#     Pauli(coeff,(p,))
# end

*(x::Real, y::Pauli) = Pauli((x < 0) ⊻ y.signbit , y.imagbit, y.bits)
*(x::Complex, y::Pauli) = Pauli(Integer(x.im) <0 , y.imagbit ⊻ true, y.bits)
*(x::Pauli, y::Number) = *(y, x)

#### Negative sign -Pauli(1, [...]) = Pauli(-1, [...])
Base.:-(x::Pauli) = Pauli(x.signbit ⊻ true, x.imagbit, x.bits)

"""
    abs(p::Pauli)
Return a new Pauli with the coefficient of `p` set to 1
"""
abs(p::Pauli) = Pauli(false, false, p.bits)

==(x::Pauli, y::Pauli) = (x.signbit == y.signbit) & isequal(bits(x), bits(y)) & (x.imagbit == y.imagbit)

"""
    @p_str("IXYZ")
    p"IXYZ..."
Macro to instantiate a Pauli{N} type from a string
"""
macro p_str(p)
    conv = Dict(
        'x' => [true, false],
        'y' => [true, true],
        'z' => [false, true],
        'i' => [false, false]
    )
    m = SparseMatrixCSC(hcat([conv[lowercase(v)] for v in p]...)')
    return Pauli(false, false, m)
end

"""
   @pauli XYZXII...

Convenience macro for constructing a Pauli{N} type.
"""
macro pauli(ex)
    s = string(ex)
    :(@p_str $s)
end

# todo: is it better to have a an empty struct or just pass around types?
# abstract type PauliPrimitive end
# for x in "IXYZ"
#     s = Symbol("Pauli$(x)")
#     constname = Symbol(lowercase(x))
#     consttype = QuoteNode(s)
#     :(struct $s <: PauliPrimitive end) |> eval
#     :($constname = eval(Expr(:call, $consttype))) |> eval
# end
# Also, might help if the internals is binary? 

⊗ = kron
function Base.kron(p1::Pauli, p2::Pauli) 
    Pauli(
        (p1.signbit ⊻ p2.signbit ⊻ (p1.imagbit & p2.imagbit)),
        (p1.imagbit ⊻ p2.imagbit), 
        vcat(p1.bits, p2.bits)
    )
end


commuting(x::Pauli, y::Pauli) = (x * y).imagbit == (x.imagbit ⊻ y.imagbit)
anticommuting(x::Pauli, y::Pauli) = ~commuting(x,y)
commuting(x::AbstractVector{<:Pauli}, y::AbstractVector{<:Pauli}) = [commuting.(i,y) for i in x]
anticommuting(x::AbstractVector{<:Pauli}, y::AbstractVector{<:Pauli}) = [anticommuting.(i,y) for i in x]

Base.length(x::Pauli) = 1
Base.iterate(x::Pauli) = (x, nothing)
Base.iterate(x::Pauli, n::Nothing) = nothing