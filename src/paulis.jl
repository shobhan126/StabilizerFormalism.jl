"""
Module defining Pauli Operators
    $(EXPORTS)
"""
module Paulis
    # Pauli Primitive Enum
    import Base: show, *, -, ==, abs, vec
    import Base.Iterators.product
    using SparseArrays: SparseMatrixCSC
    using DocStringExtensions

    include("paulitable.jl")

    """
    Create an N-Qubit Pauli Operator
        $(TYPEDEF)
    with fields:
        $(FIELDS)

    # Examples
    ```julia-repl
    julia> Pauli(:XIZI)
    XIZI
    ```
    """
    struct Pauli
        """
            true if coefficient has `-` sign
        """
        signbit::Bool
        """
            true if coefficient is imaginary
        """
        imagbit::Bool
        """
            Sparse Bit Representation of the Pauli Operator.
        """
        bits::SparseMatrixCSC
    end


    isneg(p::Pauli) = p.signbit


    """
    Create a Pauli Operator from its symplectic representation [xbits; zbits]
    
        $(FUNCTIONNAME)(xibts, zbits, neg=false, imag=false)
    
    """
    function Pauli(x::AbstractVector, z::AbstractVector, neg=false, imag=false)
        Pauli(neg, imag, SparseMatrixCSC(BitMatrix(hcat(x,z))))
    end


    """
    Create a Pauli Operator from its symplectic representation [xbits; zbits]
        
        $(FUNCTIONNAME)(xz, neg=false, imag=false)
    """
    Pauli(xz::AbstractVector, neg=false, imag=false) = Pauli(xz[1:Int(length(x)/2)], xz[Int(length(x)/2)+1:end], neg, imag)



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
    function Pauli(s::AbstractString)
        eval(:(@p_str $s))
    end


    """
    Bit Matrix Representation of Pauli Operator.

        $(FUNCTIONNAME)(p)
    """
    bits(p::Pauli) = p.bits
    xbits(p::Pauli) = bits(p)[:,1]
    zbits(p::Pauli) = bits(p)[:,2]

    """
    Take a adjoint of the pauli operator `iXYZ...`  => `-iXYZ...`.

        $(FUNCTIONNAME)(p::Pauli)
    """
    Base.adjoint(p::Pauli) = p.imagbit ? Pauli(p.signbit ⊻ true, p.imagbit, bits(p)) : p


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

    # Scalar Multiplication
    *(y::Pauli, x::Number) = Pauli((x < 0) ⊻ y.signbit , y.imagbit, y.bits)
    *(y::Pauli, x::Complex) = isequal(x,im) ? 
        Pauli(y.signbit ⊻ y.imagbit, ~y.imagbit, y.bits) : isequal(x,-im) ? 
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
    commuting(x::Pauli, y::Pauli) = (x * y).imagbit & (x.imagbit ⊻ y.imagbit)

    """
    Return true if two paulis anticommute
        $(FUNCTIONNAME)(x,y)
    """
    anticommuting(x::Pauli, y::Pauli) = ~((x * y).imagbit) & (x.imagbit ⊻ y.imagbit)

    # commuting(x::AbstractVector{<:Pauli}, y::AbstractVector{<:Pauli}) = [commuting.(i,y) for i in x]
    # anticommuting(x::AbstractVector{<:Pauli}, y::AbstractVector{<:Pauli}) = [anticommuting.(i,y) for i in x]

    # setting up iteration interface
    Base.length(x::Pauli) = 1
    Base.iterate(x::Pauli) = (x, nothing)
    Base.iterate(x::Pauli, n::Nothing) = nothing

    
    include("representations.jl")

    export Pauli, @p_str, @pauli
    export bits, xbits, zbits, symplectic, checkmatrix
    export ⊗, commuting, anticommuting
end