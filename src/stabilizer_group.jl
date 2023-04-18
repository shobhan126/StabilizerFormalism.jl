# Try storing as a the symplectic notation itself.
using StabilizerFormalism
import StabilizerFormalism: xbits, zbits, bits

mutable struct StabilizerGroup{N}
    stabs::Vector{Pauli{N}}
end

xbits(s::StabilizerGroup) = xbits.(s.stabs)
zbits(s::StabilizerGroup) = zbits.(s.stabs)
bits(s::StabilizerGroup) = hcat(xbits(s),zbits(s))

fiveQ = Pauli.([:XZZXI, :IXZZX, :XIXZZ, :ZXIXZ]) |> StabilizerGroup
xbits(fiveQ)
bits(fiveQ)