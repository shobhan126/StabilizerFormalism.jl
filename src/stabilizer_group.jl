# Try storing as a the symplectic notation itself.
import .Paulis: xbits, zbits, bits

mutable struct StabilizerGroup{N}
    stabs::Vector{Pauli{N}}
end

xbits(s::StabilizerGroup) = xbits.(s.stabs)
zbits(s::StabilizerGroup) = zbits.(s.stabs)
bits(s::StabilizerGroup) = hcat(xbits(s),zbits(s))

export StabilizerGroup
