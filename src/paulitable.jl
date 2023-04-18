const global paulibits = Dict(
    (false, false) => 1,
    (true, false) => 2,
    (true, true) => 3,
    (false, true) => 4,
)

const global psm = Bool[
    0 0 0 0
    0 0 0 1
    0 1 0 0
    0 0 1 0
]

const global pim = Bool[
    0 0 0 0
    0 0 1 1
    0 1 0 1
    0 1 1 0
]

signlogic(i, j, k, l) = @inbounds psm[paulibits[i, j], paulibits[k, l]]
imaglogic(i, j, k, l) = @inbounds pim[paulibits[i, j], paulibits[k, l]]