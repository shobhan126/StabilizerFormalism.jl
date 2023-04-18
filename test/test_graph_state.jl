using Test
using StabilizerFormalism: StabilizerGroup
using .Cliffords: CZ!
using StabilizerFormalism.Paulis
using StabilizerFormalism: GraphState

using Graphs
# Matrix from the example in Jabalizer
testGraphMatrix = Bool[ 
    1  1  1  1  1  1  0  0  0  0  0  0  0
    0  0  0  0  0  0  1  1  0  0  0  0  0
    0  0  0  0  0  0  1  0  1  0  0  0  0
    0  0  0  0  0  0  1  0  0  1  0  0  0
    0  0  0  0  0  0  1  0  0  0  1  0  0
    0  0  0  0  0  0  1  0  0  0  0  1  0
]
@testset "Graph to Stab" begin 
    gs = GraphState(6); foreach(i->add_edge!(gs, 1,i),2:6);
    zs = zbits.(StabilizerGroup(gs).stabs) |> x->hcat(x...) |> transpose
    # expcected output
    testzs = Bool[ 
    0  1  1  1  1  1
    0  1  0  0  0  0
    0  0  1  0  0  0
    0  0  0  1  0  0
    0  0  0  0  1  0
    0  0  0  0  0  1
    ];
    # rudimentary case; find more fancy examples eventually
    @test zs == testzs
end