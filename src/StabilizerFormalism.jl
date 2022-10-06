module StabilizerFormalism

    export Pauli, PauliTable, @p_str
    export symplectic, checkmatrix

    include("pauli.jl")
    include("representations.jl")

end