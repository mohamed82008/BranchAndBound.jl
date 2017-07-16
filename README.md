# BranchAndBound
## Overview

This package is under development. It currently only implements a simple depth first search (DFS) algorithm using a recursion-based implementation of the branch and bound algorithm.

For MIPs, the most fractional variable is used for branching by default and the larger bound is solved first, e.g. x >= 1 before x <= 0. The following work items are yet to be done:

1. Strong branching in the recursive implementation
2. Pseudo-cost branching in the recursive implementation
3. MIP standard cutting planes, e.g. Gromory's cut
4. Integration with jlSimplex.jl for MILPs
5. Specilization for network optimization
6. Randomness and restarts
7. A flattened implementation to allow free node selection and dynamic parallelism
8. Parallel implementation
9. Constraint propagation
10. Advanced node and variable scoring techniques

## Example

```
using BranchAndBound
using JuMP
using Clp
#using GLPKMathProgInterface

n = 10

volume = [rand()*20 for i in 1:n];
worth = [rand()*20 for i in 1:n];

m = Model(solver=ClpSolver(SolveType=0));
#m = Model(solver=GLPKSolverLP(method=:Simplex));
@variable(m, 0 <= x[1:n] <= 1, category=:Bin);
@constraint(m, sum(volume[i]*x[i] for i = 1:n) <= 20);
@objective(m, Max, sum(worth[i]*x[i] for i = 1:n));

bbproblem = BBProblem(m, mip=true);
@time my_solution, my_optimal_value = BBAlgDFS(bbproblem);
```
