using Test
using StabilizerFormalism

@testset "StabilizerFormalism.jl" begin
    @testset "Paulis" begin include("test_paulis.jl") end
    @testset "Cliffords" begin include("test_cliffords.jl") end
end
