using Test
using StabilizerFormalism

@testset "StabilizerFormalism.jl" begin
    @testset "Paulis" begin include("test_paulis.jl") end
    @testset "Cliffords" begin include("test_cliffords.jl") end
    @testset "Basic Graph States" begin include("test_graph_state.jl") end
    @testset "Stabilizer Group" begin include("test_stabilizer_group.jl") end
end
