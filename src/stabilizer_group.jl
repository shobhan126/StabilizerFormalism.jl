# Try storing as a the symplectic notation itself.
mutable struct StabilizerGroup{N}
    stabs::Vector{Pauli{N}}
end

xbits(s::StabilizerGroup) = xbits.(s.stabs)
zbits(s::StabilizerGroup) = zbits.(s.stabs)
bits(s::StabilizerGroup) = hcat(xbits(s),zbits(s))


"""
Convenience function for generating trivial state

"""
function zero_state(n) 
    p = StabilizerGroup([Pauli(zeros(Bool, n), zeros(Bool, n)) for i = 1:n])
    for i in 1:n
        p.stabs[i].zbits[i] = true
    end
    p
end


function zdiag(n) 
    p = StabilizerGroup([Pauli(zeros(Bool, n), zeros(Bool, n)) for i = 1:n])
    for i in 1:n
        p.stabs[i].zbits[i] = true
    end
    p
end


function xdiag(n) 
    p = StabilizerGroup([Pauli(zeros(Bool, n), zeros(Bool, n)) for i = 1:n])
    for i in 1:n
        p.stabs[i].xbits[i] = true
    end
    p
end

export StabilizerGroup, xdiag, zdiag
