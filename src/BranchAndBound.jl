__precompile__()

module BranchAndBound

using JuMP
using Nulls

include("types.jl")
include("dfs_simple.jl")
include("dfs_sb.jl")
include("commit_uncommit.jl")
include("util.jl")

export BBProblem, BB

end
