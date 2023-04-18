module StabilizerFormalism
    # should I load SparseArrays alongside this?
    # using SparseArrays 
    using DocStringExtensions
    include("paulis.jl")
    using .Paulis
    export AbstractPauli
    export Pauli, PauliTable, @p_str, @pauli, xbits, zbits
    export symplectic, checkmatrix
    include("stabilizer_group.jl")
    include("cliffords.jl")
    using .Cliffords: ⋊, Clifford
    export Cliffords
    
    # Maybe not pollute the namespace with clifford names
    # export  CNOT, H, S, X, Y, Z
end