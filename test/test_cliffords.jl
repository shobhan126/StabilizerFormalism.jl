using Test
using StabilizerFormalism
using .Cliffords
a = Pauli(:y) |> H(1)
x,y,z = Pauli.([:x, :y, :z])


@testset "Clifford Conjugation" begin
    x,y,z = Pauli.([:x, :y, :z]) 

    # conjugation with Hadamard
    @test x |> H(1)  == z
    @test  x |> H(1) ⋊ H(1) == x # My Dumb Operator
    @test y |> H(1) == -y
    
    # conjugation with x
    @test x |> X(1) == x
    @test z |> X(1) == -z
    @test y |> X(1) == -y


    # conjugation with s
    @test x |> S(1) == y
    @test z |> S(1) == z
    @test y |> S(1) == x


    # conjugation with y
    @test x |> Y(1) == -x


    # cnot conjugation
    @test p"XI" |> CNOT(1,2) == p"XX"
    @test p"IX" |> CNOT(1,2) == p"IX"
    @test p"ZI" |> CNOT(1,2) == p"ZI"
    @test p"IZ" |> CNOT(1,2) == p"ZZ"
    @test p"IY" |> CNOT(1,2) == p"ZY"
    @test p"YI" |> CNOT(1,2) == p"YX"

end