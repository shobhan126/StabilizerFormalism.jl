# TODO Autogenerate levi-civita dictionary indexed by tuples and keep it in the namespace.
const global paulibits = Dict(
    (false, false) => 0,
    (true, false) => 1,
    (true, true) => 2,
    (false, true) => 3,
)

const global paulibits2 = Dict(
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

signlogic(i, j, k, l) = @inbounds psm[paulibits2[i, j], paulibits2[k, l]]
imaglogic(i, j, k, l) = @inbounds pim[paulibits2[i, j], paulibits2[k, l]]

signlogic(i, j) = psm[paulibits2[i[1], i[2]], paulibits2[j[1], j[2]]]
imaglogic(i, j) = pim[paulibits2[i[1], i[2]], paulibits2[j[1], j[2]]]