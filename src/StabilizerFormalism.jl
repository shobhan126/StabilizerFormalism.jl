module StabilizerFormalism
    # should I load SparseArrays alongside this?
    # using SparseArrays 
    using DocStringExtensions
    include("paulis.jl")
    using .Paulis
    export Pauli, PauliTable, @p_str, @pauli
    export symplectic, checkmatrix
    include("cliffords.jl")
    export Cliffords
    using .Cliffords: ⋊, Clifford
    # Maybe not pollute the namespace with clifford names
    # export  CNOT, H, S, X, Y, Z
end