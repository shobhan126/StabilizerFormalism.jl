module StabilizerFormalism

    export Pauli, PauliTable, @p_str, @pauli
    export symplectic, checkmatrix

    include("pauli.jl")
    include("representations.jl")

end