using Test
using StabilizerFormalism

@testset "Stabilizer Construction" begin
    fiveQ = Pauli.([:XZZXI, :IXZZX, :XIXZZ, :ZXIXZ]) |> StabilizerGroup
    @test string.(fiveQ.stabs) == ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
end

# maybe define a generic array with Paulis.
# todo define adjoint operation and corresponding test.


