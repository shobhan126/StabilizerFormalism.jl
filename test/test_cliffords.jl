using Test
using StabilizerFormalism
using .Cliffords
a = Pauli(:y) ⋊ H(1) ⋊ H(1) == Pauli(:y) 

x,y,z = Pauli.([:x, :y, :z])


@testset "Clifford Conjugation" begin
    x,y,z = Pauli.([:x, :y, :z]) 

    @test x ⋊ H(1)  == z
    @test x ⋊ H(1) ⋊ H(1) == x
    @test y ⋊ H(1) == -y
    
end