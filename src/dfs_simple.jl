
function BB(root_problem::BBProblem{T,S,C,:Max}; alg=:dfs_sb) where {T,S,C<:AbstractChoice}
    if alg == :dfs_simple
        return dfs_simple(root_problem)
    elseif alg == :dfs_sb
        return dfs_sb(root_problem)
    end
end

"""
Applies a depth-first search branch and bound algorithm on `root_problem`.
Returns `solution::Nullable{S}`, `optimal_value::Float64`
If the problem is `:infeasible`, `solution` will be empty.
"""
function dfs_simple(root_problem::BBProblem{T,S,C,:Max}) where {T,S,C<:AbstractChoice}
    root_problem.solution, root_problem.optimal_value = solve(root_problem)

    if root_problem.solution.hasvalue
        root_problem.status = isfeasible(root_problem) ? (:feasible) : (:infeasible)
    end
    if root_problem.status == :infeasible && root_problem.optimal_value.value > root_problem.inner_bound.x
        dfs_simple_recur(root_problem)
    elseif root_problem.status == :feasible
        root_problem.inner_bound.x = root_problem.optimal_value.value
        root_problem.best_solution.x = root_problem.solution
        root_problem.outer_bound = root_problem.optimal_value.value
    end

    if root_problem.best_solution.x.hasvalue
        root_problem.status = :optimal
        root_problem.optimal_value = Nullable{Float64}(root_problem.inner_bound.x)
        root_problem.solution = root_problem.best_solution.x
    else
        root_problem.status = :infeasible
    end

    #Returns a (Nullable{S}, Float64)
    if root_problem.status == :infeasible
        println("Problem is infeasible.")
    end
    return root_problem.best_solution.x, root_problem.inner_bound.x
end    

function dfs_simple(root_problem::BBProblem{T,S,C,:Min}) where {T,S,C<:AbstractChoice}
    root_problem.solution, root_problem.optimal_value = solve(root_problem)

    if root_problem.solution.hasvalue
        root_problem.status = isfeasible(root_problem) ? (:feasible) : (:infeasible)
    end
    if root_problem.status == :infeasible && root_problem.optimal_value.value < root_problem.inner_bound.x
        dfs_simple_recur(root_problem)
    elseif root_problem.status == :feasible
        root_problem.inner_bound.x = root_problem.optimal_value.value
        root_problem.best_solution.x = root_problem.solution
        root_problem.outer_bound = root_problem.optimal_value.value
    end

    if root_problem.best_solution.x.hasvalue
        root_problem.status = :optimal
        root_problem.optimal_value = Nullable{Float64}(root_problem.inner_bound.x)
        root_problem.solution = root_problem.best_solution.x
    else
        root_problem.status = :infeasible
    end

    #Returns a (Nullable{S}, Float64)
    if root_problem.status == :infeasible
        println("Problem is infeasible.")
    end
    return root_problem.best_solution.x, root_problem.inner_bound.x
end    

function dfs_simple_recur(root_problem::BBProblem{T,S,C,:Max}) where {T,S,C<:AbstractChoice}
    choices = getchoices(root_problem)
    children_outer_bounds = Float64[]
    for choice in choices
        if root_problem.optimal_value.value > root_problem.inner_bound.x #Solve child only if it can improve the best solution found
            child_problem = make_child!(root_problem, choice)
            child_problem.solution, child_problem.optimal_value = solve(child_problem)
            if child_problem.solution.hasvalue
                child_problem.status = isfeasible(child_problem) ? (:feasible) : (:infeasible)
                child_problem.outer_bound = child_problem.optimal_value.value
                push!(children_outer_bounds, child_problem.outer_bound)
                if child_problem.status == :feasible && child_problem.optimal_value.value > root_problem.inner_bound.x
                    #Lower bound and best solution update
                    root_problem.inner_bound.x = child_problem.optimal_value.value
                    root_problem.best_solution.x = child_problem.solution
                elseif child_problem.status == :infeasible && child_problem.optimal_value.value > root_problem.inner_bound.x
                    dfs_simple_recur(child_problem)
                end
            end
            #Removing last choice added before moving on, because all problems share the same choices vector to save memory.
            pop_child!(root_problem)
        else
            root_problem.outer_bound = -Inf
            break
        end
    end

    if root_problem.outer_bound > -Inf
        if length(children_outer_bounds) > 0
            root_problem.outer_bound = max(children_outer_bounds...)
        else
            root_problem.outer_bound = -Inf
        end
    end
    return
end

function dfs_simple_recur(root_problem::BBProblem{T,S,C,:Min}) where {T,S,C<:AbstractChoice}
    choices = getchoices(root_problem)
    children_outer_bounds = Float64[]
    for choice in choices
        if root_problem.optimal_value.value < root_problem.inner_bound.x #Solve child only if it can improve the best solution found
            child_problem = make_child!(root_problem, choice)
            child_problem.solution, child_problem.optimal_value = solve(child_problem)
            if child_problem.solution.hasvalue
                child_problem.status = isfeasible(child_problem) ? (:feasible) : (:infeasible)
                child_problem.outer_bound = child_problem.optimal_value.value
                push!(children_outer_bounds, child_problem.outer_bound)
                if child_problem.status == :feasible && child_problem.optimal_value.value < root_problem.inner_bound.x
                    #Lower bound and best solution update
                    root_problem.inner_bound.x = child_problem.optimal_value.value
                    root_problem.best_solution.x = child_problem.solution
                elseif child_problem.status == :infeasible && child_problem.optimal_value.value < root_problem.inner_bound.x
                    dfs_simple_recur(child_problem)
                end
            end
            #Removing last choice added before moving on, because all problems share the same choices vector to save memory.
            pop_child!(root_problem)
        else
            root_problem.outer_bound = Inf
            break
        end
    end

    if root_problem.outer_bound < Inf
        if length(children_outer_bounds) > 0
            root_problem.outer_bound = min(children_outer_bounds...)
        else
            root_problem.outer_bound = Inf
        end
    end
    return
end

function make_child!(problem::BBProblem{T,S,C,Sense}, choice::C) where {T,S,C<:AbstractChoice,Sense}
    problem.relaxation.x.model.colVal = problem.solution.value
    child_problem = BBProblem{T,S,C,Sense}(problem.relaxation, 
        inner_bound = problem.inner_bound, 
        outer_bound = problem.optimal_value.value, 
        status = :unsolved, 
        choices = problem.choices, 
        solution = Nullable{S}(), 
        optimal_value = Nullable{Float64}(), 
        best_solution = problem.best_solution,
        isfeasible=problem.isfeasible,
        solve=problem.solve,
        getchoices=problem.getchoices)

    push!(child_problem.choices, choice)
    commit(problem, choice)
    return child_problem
end

function pop_child!(problem::BBProblem{T,S,C,Sense}) where {T,S,C<:AbstractChoice,Sense}
    choice = pop!(problem.choices)
    uncommit(problem, choice)
    return problem
end

function solvejump(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}) where {Sense}
    status = JuMP.solve(model(problem))
    if status == :Optimal
        return Nullable(model(problem).colVal), Nullable(model(problem).objVal)
    else
        return Nullable{Vector{Float64}}(), Nullable{Float64}()
    end
end

function getchoicesmip(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}) where {Sense}
    solution = problem.solution.value
    infeasvarinds = problem.infeasiblecols
    lbs = floor.(solution[infeasvarinds])
    resid = solution[infeasvarinds] .- lbs
    v,i = findmin(norm.(resid .- 0.5))
    ind = infeasvarinds[i]
    lb = lbs[i]
    ub = ceil(solution[ind])
    if ub < model(problem).colUpper[ind] - tol(problem)
        choice1 = JuMPVarChoice(problem.relaxation.x, ind, :(>=), ub)
    else
        choice1 = JuMPVarChoice(problem.relaxation.x, ind, :(=), ub)
    end
    if lb > model(problem).colLower[ind] + tol(problem)
        choice2 = JuMPVarChoice(problem.relaxation.x, ind, :(<=), lb)
    else
        choice2 = JuMPVarChoice(problem.relaxation.x, ind, :(=), lb)
    end
    return AbstractChoice{:JuMP}[choice1, choice2]
end

function isfeasiblemip(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}) where {Sense}
    solution = problem.solution.value
    varcat = problem.relaxation.x.varcat
    for i in 1:length(solution)
        if varcat[i] == :Bin
            if (solution[i] > tol(problem)) && (solution[i] < 1-tol(problem))
                push!(problem.infeasiblecols, i)
            end
        elseif varcat[i] == :Int
            resid = solution[i] - floor(solution[i])
            if (resid > tol(problem)) && (resid < 1-tol(problem))
                push!(problem.infeasiblecols, i)
            end
        end
    end
    if length(problem.infeasiblecols) == 0
        return true
    else
        return false
    end
end

function solve(problem::BBProblem{T,S,C,Sense}) where {T,S,C<:AbstractChoice,Sense}
    problem.solve.x(problem)
end
function isfeasible(problem::BBProblem{T,S,C,Sense}) where {T,S,C<:AbstractChoice,Sense}
    problem.isfeasible.x(problem)
end
function getchoices(problem::BBProblem{T,S,C,Sense}) where {T,S,C<:AbstractChoice,Sense}
    problem.getchoices.x(problem)
end
