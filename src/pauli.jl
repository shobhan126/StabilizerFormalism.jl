# Pauli Primitive Enum
import Base: show, *, -, ==, abs, vec
import Base.Iterators.product

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
struct Pauli{N}
    coeff::Number
    val::NTuple{N, PauliPrimitive}
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

"""
    vec(s::Pauli)

Return the vector of PauliPrimitives that define the Pauli{N} type. 
"""
vec(p::Pauli) =  [p.val...]


##### Defining Pauli Vector Product  ###### 

# defining product of two Paulis using Table Lookup
*(p::PauliPrimitive, q::PauliPrimitive) = Pauli(PauliTable[p, q]...)

# Defining Pauli{1} times PauliPrimitive
function *(p::Pauli{1}, q::PauliPrimitive) 
    r = p.val[1] * q
    Pauli(r.val...; coeff= r.coeff * p.coeff,)
end

*(x::PauliPrimitive, y::Pauli{1}) = *(y, x)
*(x::Number, y::Pauli) = Pauli(y.coeff*x, y.val)
*(x::Pauli, y::Number) = *(y, x)

"""
Product between two Pauli Operators
"""
function *(p::Pauli{N}, q::Pauli{N}) where N
    r = p.val .* q.val
    Pauli(prod(getfield.(r, :coeff)), first.(getfield.(r, :val)))
end

#### Negative sign -Pauli(1, [...]) = Pauli(-1, [...])
-(x::Pauli) = Pauli(-x.coeff, vec(x))

"""
    abs(p::Pauli)
Return a new Pauli with the coefficient of `p` set to 1
"""
abs(p::Pauli) = Pauli(abs(p.coeff), vec(x))

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

"""
    @p_str("IXYZ")
    p"IXYZ..."
Macro to instantiate a Pauli{N} type from a string
"""
macro p_str(p) 
    ops = eval.([Symbol("σ"*lowercase(v)) for v in p])
    return Pauli(ops...)
end

"""
   @pauli XYZXII...

Convenience macro for constructing a Pauli{N} type.
"""
macro pauli(ex)
    s = string(ex)
    :(@p_str $s)
end


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

