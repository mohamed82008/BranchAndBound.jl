__precompile__()

module BranchAndBound

using JuMP
using Nulls

include("types.jl")
include("dfs_algs.jl")
include("commit_uncommit.jl")
include("util.jl")

export BBProblem, BBAlgDFS, JuMPProblem

end
