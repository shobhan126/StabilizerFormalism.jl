
@enum PauliPrimitive σi=0 σx=1 σy=2 σz=3
# helpers to create group table
# anti-syemmtric tensor
ε = [
    0 1 -1;
    -1 0 1;
    1 -1 0
] * im

# ϵ(i::Int,j::Int) = (0 in [i, j]) | (i==j) ? 1 : ε[i,j]

global const paulibit = Dict(
    (false, false) => 0,
    (true, false) => 1,
    (true, true) => 2,
    (false, true) => 3,
)

ϵ(i::AbstractArray, j::AbstractArray) =  iszero(i) | iszero(j) | (i==j) ? 1 : ε[paulibit[i...], paulibits[j...]]

k(i,j) = (i == j) ? 0 : 0 in [i,j] ? max(i,j) : abs(i-j) + (max(i,j) % 3)

# Constructing Group Table
# global const PauliTable = Dict(
#     [(PauliPrimitive(i), PauliPrimitive(j)) => (ϵ(i,j), PauliPrimitive(k(i,j))) for (i,j) in product(0:3, 0:3)]
# )