"""
Implements a Strong Branching depth-first search branch and bound algorithm on `root_problem`.
Returns `solution::Nullable{S}`, `optimal_value::Float64`
If the problem is `:infeasible`, `solution` will be empty.
"""
function dfs_sb(root_problem::BBProblem{T,S,C,:Max}) where {T,S,C<:AbstractChoice}
    root_problem.solution, root_problem.optimal_value = solve(root_problem)
    if root_problem.solution.hasvalue
        root_problem.status = isfeasible(root_problem) ? (:feasible) : (:infeasible)
    end
    if root_problem.status == :infeasible && root_problem.optimal_value.value > root_problem.inner_bound.x
        bb_dfs_sb_recur(root_problem)
    elseif root_problem.status == :feasible
        root_problem.inner_bound.x = root_problem.optimal_value.value
        root_problem.best_solution.x = root_problem.solution
        root_problem.outer_bound = root_problem.optimal_value.value
    end

    if root_problem.best_solution.x.hasvalue
        root_problem.status = :optimal
        root_problem.optimal_value = Nullable{Float64}(root_problem.inner_bound.x)
        root_problem.outer_bound = root_problem.inner_bound.x
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

function bb_dfs_sb_recur(root_problem::BBProblem{T,S,C,:Max}) where {T,S,C<:AbstractChoice}
    selectedchildren, choices, children_outer_bounds = select_children_sb(root_problem)
    if length(selectedchildren) == 0
        root_problem.outer_bound = -Inf
        return
    end
    root_problem.outer_bound = max(children_outer_bounds...)
    for i in 1:length(selectedchildren)
        if children_outer_bounds[i] > root_problem.inner_bound.x
            child_problem = selectedchildren[i]
            choice = choices[i]
            revive_child!(child_problem, choice)
            bb_dfs_sb_recur(child_problem)
            pop_child!(root_problem)
        end
    end
    root_problem.outer_bound = -Inf
    return
end

function select_children_sb(root_problem::BBProblem{T,S,C,:Max}) where {T,S,C<:AbstractChoice}
    allpossiblechoices = getchoicesmip_sb(root_problem)
    ncs = length(allpossiblechoices[1])
    allchildrenproblems = BBProblem{T,S,C,:Max}[]
    allouterbounds = Float64[]
    for choiceset in allpossiblechoices
        for choice in choiceset
            child_problem = make_child!(root_problem, choice)
            child_problem.solution, child_problem.optimal_value = solve(child_problem)
            #Recording outer bounds of potential child node pairs
            if child_problem.solution.hasvalue
                child_problem.status = isfeasible(child_problem) ? (:feasible) : (:infeasible)
                #Lower bound and best solution update
                if child_problem.status == :feasible && child_problem.optimal_value.value > root_problem.inner_bound.x
                    root_problem.inner_bound.x = child_problem.optimal_value.value
                    root_problem.best_solution.x = child_problem.solution
                    child_problem.outer_bound = -Inf
                elseif child_problem.status == :infeasible && child_problem.optimal_value.value > root_problem.inner_bound.x
                    child_problem.outer_bound = child_problem.optimal_value.value
                else
                    child_problem.outer_bound = -Inf
                end
            else
                child_problem.outer_bound = -Inf
            end
            push!(allouterbounds, child_problem.outer_bound)
            push!(allchildrenproblems, pop_child!(child_problem))
        end
    end

    #Finding mean upper bound and choosing best variable to branch on
    meanubs = Float64[]
    for i in 1:ncs:length(allchildrenproblems)
        obs = allouterbounds[i:i+ncs-1]
        _count = sum(Int.(obs .> -Inf))
        if _count == 0
            push!(meanubs, NaN)
        else
            push!(meanubs, mean(obs .^(obs .> -Inf)))
        end
    end
    v,i = findmin(meanubs)
    if v == NaN
        return BBProblem{T,S,C,:Max}[], C[], Float64[]
    else
        selectedchildren, selectedchoices, selectedouterbounds = 
            BBProblem{T,S,C,:Max}[], C[], Float64[]
        for j in 1:ncs
            if allouterbounds[j+ncs*(i-1)] > root_problem.inner_bound.x
                push!(selectedchildren, allchildrenproblems[j+ncs*(i-1)])
                push!(selectedchoices, allpossiblechoices[i][j])
                push!(selectedouterbounds, allouterbounds[j+ncs*(i-1)])
            end
        end
        return selectedchildren, selectedchoices, selectedouterbounds
    end
end

function dfs_sb(root_problem::BBProblem{T,S,C,:Min}) where {T,S,C<:AbstractChoice}
    root_problem.solution, root_problem.optimal_value = solve(root_problem)
    if root_problem.solution.hasvalue
        root_problem.status = isfeasible(root_problem) ? (:feasible) : (:infeasible)
    end
    if root_problem.status == :infeasible && root_problem.optimal_value.value < root_problem.inner_bound.x
        bb_dfs_sb_recur(root_problem)
    elseif root_problem.status == :feasible
        root_problem.inner_bound.x = root_problem.optimal_value.value
        root_problem.best_solution.x = root_problem.solution
        root_problem.outer_bound = root_problem.optimal_value.value
    end

    if root_problem.best_solution.x.hasvalue
        root_problem.status = :optimal
        root_problem.optimal_value = Nullable{Float64}(root_problem.inner_bound.x)
        root_problem.outer_bound = root_problem.inner_bound.x
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

function bb_dfs_sb_recur(root_problem::BBProblem{T,S,C,:Min}) where {T,S,C<:AbstractChoice}
    selectedchildren, choices, children_outer_bounds = select_children_sb(root_problem)
    if length(selectedchildren) == 0
        root_problem.outer_bound = Inf
        return
    end
    root_problem.outer_bound = min(children_outer_bounds...)
    for i in 1:length(selectedchildren)
        if children_outer_bounds[i] < root_problem.inner_bound.x
            child_problem = selectedchildren[i]
            choice = choices[i]
            revive_child!(child_problem, choice)
            bb_dfs_sb_recur(child_problem)
            pop_child!(root_problem)
        end
    end
    root_problem.outer_bound = Inf
    return
end

function select_children_sb(root_problem::BBProblem{T,S,C,:Min}) where {T,S,C<:AbstractChoice}
    allpossiblechoices = getchoicesmip_sb(root_problem)
    ncs = length(allpossiblechoices[1])
    allchildrenproblems = BBProblem{T,S,C,:Min}[]
    allouterbounds = Float64[]
    for choiceset in allpossiblechoices
        for choice in choiceset
            child_problem = make_child!(root_problem, choice)
            child_problem.solution, child_problem.optimal_value = solve(child_problem)
            #Recording outer bounds of potential child node pairs
            if child_problem.solution.hasvalue
                child_problem.status = isfeasible(child_problem) ? (:feasible) : (:infeasible)
                #Lower bound and best solution update
                if child_problem.status == :feasible && child_problem.optimal_value.value < root_problem.inner_bound.x
                    root_problem.inner_bound.x = child_problem.optimal_value.value
                    root_problem.best_solution.x = child_problem.solution
                    child_problem.outer_bound = Inf
                elseif child_problem.status == :infeasible && child_problem.optimal_value.value < root_problem.inner_bound.x
                    child_problem.outer_bound = child_problem.optimal_value.value
                else
                    child_problem.outer_bound = Inf
                end
            else
                child_problem.outer_bound = Inf
            end
            push!(allouterbounds, child_problem.outer_bound)
            push!(allchildrenproblems, pop_child!(child_problem))
        end
    end

    #Finding mean upper bound and choosing best variable to branch on
    meanubs = Float64[]
    for i in 1:ncs:length(allchildrenproblems)
        obs = allouterbounds[i:i+ncs-1]
        _count = sum(Int.(obs .< Inf))
        if _count == 0
            push!(meanubs, NaN)
        else
            push!(meanubs, mean(obs .^(obs .< Inf)))
        end
    end
    v,i = findmin(meanubs)
    if v == NaN
        return BBProblem{T,S,C,:Min}[], C[], Float64[]
    else
        selectedchildren, selectedchoices, selectedouterbounds = 
            BBProblem{T,S,C,:Min}[], C[], Float64[]
        for j in 1:ncs
            if allouterbounds[j+ncs*(i-1)] < root_problem.inner_bound.x
                push!(selectedchildren, allchildrenproblems[j+ncs*(i-1)])
                push!(selectedchoices, allpossiblechoices[i][j])
                push!(selectedouterbounds, allouterbounds[j+ncs*(i-1)])
            end
        end
        return selectedchildren, selectedchoices, selectedouterbounds
    end
end

function revive_child!(problem::BBProblem{T,S,C,Sense}, choice::AbstractChoice) where {T,S,C<:AbstractChoice,Sense}
    push!(problem.choices, choice)
    commit(problem, choice)
    return
end

function getchoicesmip_sb(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}) where {Sense}
    solution = problem.solution.value
    infeasvarinds = problem.infeasiblecols
    lbs = floor.(solution[infeasvarinds])
    allpossiblechoices = Vector{AbstractChoice{:JuMP}}[]
    for i in 1:length(infeasvarinds)
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
        push!(allpossiblechoices, AbstractChoice{:JuMP}[choice1, choice2])
    end
    return allpossiblechoices
end
