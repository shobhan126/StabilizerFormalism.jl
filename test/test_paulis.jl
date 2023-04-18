using Test
using StabilizerFormalism
using StabilizerFormalism.Paulis


@testset "Pauli Operator Construction" begin
    
    x = @pauli XIYZXXXIXIIIZZZX

    @test x == Pauli(:XIYZXXXIXIIIZZZX)
    @test x == Pauli("XIYZXXXIXIIIZZZX")
    @test x == Pauli(lowercase("XIYZXXXIXIIIZZZX"))
end

x,y,z = Pauli.([:x,:y,:z])

@testset "Multiplicative Operations" begin
    # scalar multiplication
    imX = (im * x)
    im2X = im * imX
    @test imX.imagbit 
    @test ~im2X.imagbit & im2X.signbit
    @test -x == im2X 
    @test x * x == y * y == z * z
    @test x * y == im * Pauli(:z)
    @test y * z == im * Pauli(:x)
    @test z * x == im * Pauli(:y)

end

@testset "Tensor Operation" begin
    # test tensor product
    @test x ⊗ x == Pauli(:xx)
    @test im*Pauli(:x) ⊗ Pauli(:y) == im * Pauli(:xy)
end

@testset "Adjoint Operation" begin
    # check adjoint 
    @test (im * x)' == -im * x
end