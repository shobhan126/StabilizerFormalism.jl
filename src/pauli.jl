# Pauli Primitive Enum
import Base: show, *, -, ==, abs, vec
import Base.Iterators.product

@enum PauliPrimitive σi=0 σx=1 σy=2 σz=3

# should 
# abstract type PauliPrimitive end
# for x in "IXYZ"
#     s = Symbol("Pauli$(x)")
#     constname = Symbol(lowercase(x))
#     consttype = QuoteNode(s)
#     :(struct $s <: PauliPrimitive end) |> eval
#     :($constname = eval(Expr(:call, $consttype))) |> eval
# end

# helpers to create group table
# anti-syemmtric tensor

ε = [
    0 1 -1;
    -1 0 1;
    1 -1 0
] * im

ϵ(i::Int,j::Int) = (0 in [i, j]) | (i==j) ? 1 : ε[i,j]
k(i,j) = (i == j) ? 0 : 0 in [i,j] ? max(i,j) : abs(i-j) + (max(i,j) % 3)

# Group Table
global const PauliTable = Dict(
    [(PauliPrimitive(i), PauliPrimitive(j)) => (ϵ(i,j), PauliPrimitive(k(i,j))) for (i,j) in product(0:3, 0:3)]
)

"""
N-Qubit Pauli Group
"""
struct Pauli{N}
    coeff::Number
    val::NTuple{N, PauliPrimitive}
end

function Pauli(pvec::PauliPrimitive...; coeff=1)
    #  Only allowing 4 possible coefficients
    @assert coeff in [1, -1, im, -im] "coeff can only be 1, -1, im, or -im"
    for p in pvec
        @assert p isa PauliPrimitive
    end
    Pauli(coeff, Tuple(pvec))
end

function Pauli(coeff::Number, p::PauliPrimitive)
    @assert coeff in [1, -1, im, -im]
    @assert p isa PauliPrimitive
    Pauli(coeff,(p,))
end

vec(s::Pauli) =  [s.val...]

##### Defining Pauli Vector Product  ###### 

# defining product of two Paulis using Table Lookup
*(p::PauliPrimitive, q::PauliPrimitive) = Pauli(PauliTable[p, q]...)

# Defining Pauli{1} times Pauli
function *(p::Pauli{1}, q::PauliPrimitive) 
    r = p.val[1] * q
    Pauli(r.val...; coeff= r.coeff * p.coeff,)
end

*(x::PauliPrimitive, y::Pauli{1}) = *(y, x)
*(x::Number, y::Pauli) = Pauli(y.coeff*x, y.val)
*(x::Pauli, y::Number) = *(y, x)

function *(p::Pauli{N}, q::Pauli{N}) where N
    r = p.val .* q.val
    Pauli(prod(getfield.(r, :coeff)), first.(getfield.(r, :val)))
end

#### Negative sign -Pauli(1, [...]) = Pauli(-1, [...])
-(x::Pauli) = Pauli(-x.coeff, vec(x))
abs(x::Pauli) = Pauli(abs(x.coeff), vec(x))
==(x::Pauli, y::Pauli) = isequal(x.coeff, y.coeff) & isequal(vec(x), vec(y))
    

# REPL representation for the Pauli Operator
function Base.show(io::IO, p::Pauli)
    if p.coeff in [-1, -im]
        prefix = "-"
    else
        prefix = ""
    end
    if p.coeff in [im, -im]
        prefix*="im"
    end
    print(io, prefix)
    [print(io, uppercase(last(string(i)))) for i in p.val]
end


macro p_str(p) 
    ops = eval.([Symbol("σ"*lowercase(v)) for v in p])
    return Pauli(ops...)
end

"""
Convenience macro for constructing a Pauli Operator.
@pauli XYZXII
"""
macro pauli(ex)
    s = string(ex)
    :(@p_str $s)
end