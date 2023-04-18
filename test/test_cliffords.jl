using Test
using StabilizerFormalism
using .Cliffords
a = Pauli(:y) |> H(1)
x,y,z = Pauli.([:x, :y, :z])


@testset "Clifford Conjugation" begin
    x,y,z = Pauli.([:x, :y, :z]) 
    @test x |> H(1)  == z
    @test  x |> H(1) ⋊ H(1) == x # My Dumb Operator
    @test y |> H(1) == -y
    
end