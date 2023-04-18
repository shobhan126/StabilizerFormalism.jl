
using Graphs
using LoopVectorization: @turbo, @inbounds
export GraphState, add_edge!, add_vertex!, nv, ne, edges, src, dst
# Just use the simple graph for now; 
const GraphState = SimpleGraph
using .Cliffords

# add CZs for all the edges
# maybe collect the edges and do a parallel operation?
function StabilizerGroup(gs::SimpleGraph)
    stabs = xdiag(nv(gs))
    for e in edges(gs)
        # broadcast over all stabilizers
        CZ!(src(e), dst(e)).(stabs.stabs)
    end
    stabs
end

# adjacency matrix; just create a new function that initializes x!
# infact it is better to use an adjacency list! instead of a matrix for this
# no need for an adjacency matrix as the basic representation
# function adjm2state(A::AbstractMatrix{Bool})
#     n = size(A, 1)
#     s = xdiag(n) # start each qubit in X
#     # apply CZ with upper triangular adjacency matrix if there exists 
#     @inbounds for i in 1:n
#         for j in (i+1):n
#             A[i, j] ? CZ!(i, j).(s.stabs) : nothing
#         end
#     end
#     s
# end

# @time adjm2state(testGraphMatrix)

# checkmatrix(adjm2state(testGraphMatrix).stabs)