############ Choices #############

abstract type AbstractChoice{T} end
abstract type AbstractJuMPVarChoice <: AbstractChoice{:JuMP} end
abstract type AbstractJuMPLinearConstrChoice <: AbstractChoice{:JuMP} end

immutable JuMPChoice <: AbstractChoice{:JuMP}
    subchoices::Vector{AbstractChoice{:JuMP}}
end
immutable JuMPVarChoiceLeq <: AbstractJuMPVarChoice
    varind::Int
    boundbefore::Float64
    boundafter::Float64
end
immutable JuMPVarChoiceGeq <: AbstractJuMPVarChoice
    varind::Int
    boundbefore::Float64
    boundafter::Float64
end
immutable JuMPVarChoiceEq <: AbstractJuMPVarChoice
    varind::Int
    varcat::Symbol
    boundsbefore::Tuple{Float64, Float64}
    boundafter::Float64
end
immutable JuMPNewLinearConstrChoice <: AbstractJuMPLinearConstrChoice
    ind::Int
    constr::JuMP.LinearConstraint
end
immutable JuMPStaticLinearConstrChoice <: AbstractJuMPLinearConstrChoice
    ind::Int
    coeffsbefore::Vector{Float64}
    constantbefore::Float64
    boundsbefore::Tuple{Float64,Float64}
    coeffsafter::Vector{Float64}
    constantafter::Float64
    boundsafter::Tuple{Float64,Float64}
end
immutable JuMPDynamicLinearConstrChoice <: AbstractJuMPLinearConstrChoice
    ind::Int
    varsbefore::Vector{Union{Symbol,Expr}}
    coeffsbefore::Vector{Float64}
    constantbefore::Float64
    boundsbefore::Tuple{Float64,Float64}
    varsafter::Vector{Union{Symbol,Expr}}
    coeffsafter::Vector{Float64}
    constantafter::Float64
    boundsafter::Tuple{Float64,Float64}
end

############ Branch and bound problem #############

"""
`BBProblem` is a type for branch and bound problems.
`T`: relaxation problem type, can be used to hold multiple representations of the problem where each is compatible with some algorithm/heuristic.
`S`: solution type.
`S`: choice type.
`Sense`: `:Max` or `:Min` for maximization or minimization.

`relaxation`: this is a reference to the problem after relaxing some of the constraints, e.g. integrality constraints.
`choices`: a vector of choices that have been committed to reach a current branch.
`status`: a Symbol that is `:unsolved` if the relaxation is not attempted yet or if is attempted and no solution was found.
    It is `:feasible` if the relaxation is solved and a solution feasible to the unrelaxed problem was found.
    It is `:infeasible` if the relaxation is solved and a solution that is feasible to relaxation but not to the unrelaxed problem was found, hence branching might be due.
`solution`: the solution found after solving the current branch's problem, is a Nullable.
`optimal_value`: the objective value associated with the solution found to the current branch.
`inner_bound`: lower bound in a maximization problem and upper bound in a minimization problem. This the objective value of the best solution found, shared between all the instances.
`outer_bound`: upper bound in a maximization problem and lower bound in a minimization problem.
`best_solution`: the best solution found during the search process, shared between all the instances.
"""
mutable struct BBProblem{T,S,C<:AbstractChoice,Sense}
    relaxation::Ref{T} #Shared among all problems, modified for each branch
    infeasiblecols::Vector{Int}
    inner_bound::Ref{Float64} #Shared among all problems
    outer_bound::Float64
    status::Symbol
    choices::Vector{C} #Shared in any one search path
    solution::Nullable{S}
    optimal_value::Nullable{Float64}
    best_solution::Ref{Nullable{S}} #Shared among all problems
    isfeasible::Ref{Function} #Shared among all problems
    solve::Ref{Function} #Shared among all problems
    getchoices::Ref{Function} #Shared among all problems

    function BBProblem{T,S,C,Sense}(relaxation::Ref{T}; 
        inner_bound::Ref{Float64} = (Sense == :Max ? Ref(-Inf) : Ref(Inf)), 
        outer_bound::Float64 = (Sense == :Max ? Inf : -Inf), 
        status::Symbol = :unsolved,
        choices::Vector{C}=Vector{C}(), 
        solution::Nullable{S}=Nullable{S}(), 
        optimal_value::Nullable{Float64}=Nullable{Float64}(),
        best_solution::Ref{Nullable{S}}=Ref(Nullable{S}()), 
        isfeasible::Ref{Function}=Ref{Function}((problem::BBProblem)->()), 
        solve::Ref{Function}=Ref{Function}((problem::BBProblem)->()),
        getchoices::Ref{Function}=Ref{Function}((problem::BBProblem)->())) where {T,S,C<:AbstractChoice,Sense}

        return new{T,S,C,Sense}(relaxation, Int[], inner_bound, outer_bound, status, choices, solution, 
            optimal_value, best_solution, isfeasible, solve, getchoices)
    end
end

############ JuMP problem #############

mutable struct JuMPProblem
    model::JuMP.Model
    varcat::Vector{Symbol}
    tol::Float64
end

########### Constructors ###########

function BBProblem(problem::T, SolutionType::DataType, ChoiceType::DataType, 
    Sense::Symbol, isfeasible, solve, getchoices) where {T}
    BBProblem{T,SolutionType,ChoiceType,Sense}(Ref{T}(problem), isfeasible=Ref{Function}(isfeasiblemip), 
        solve=Ref{Function}(solvejump), getchoices=Ref{Function}(getchoicesmip))
end

function JuMPProblem(model::JuMP.Model; tol::Float64=1E-6, mip=false)
    if mip
        varcat = copy(model.colCat)
        map((i)->((model.colCat[i] == :Int || model.colCat[i] == :Bin) && (model.colCat[i] = :Cont); return), 
            1:length(model.colCat))
    else
        varcat = model.colCat
    end        
    return JuMPProblem(model, varcat, tol)
end

BBProblem(problem::JuMPProblem) = BBProblem(problem, Vector{Float64}, AbstractChoice{:JuMP}, 
    sense(problem), isfeasible, solve, getchoices)
function BBProblem(model::JuMP.Model; mip=false)
    problem = JuMPProblem(model, mip=mip)
    return BBProblem(problem, Vector{Float64}, AbstractChoice{:JuMP}, sense(problem), 
        isfeasible, solve, getchoices)
end

JuMPChoice() = JuMPChoice(Vector{AbstractChoice{:JuMP}}())
JuMPChoice(problem::JuMPProblem, newlinconstrs::Vector{JuMP.LinearConstraint}) = 
    JuMPLinearConstrChoice(problem, newlinconstrs)
JuMPChoice(problem::JuMPProblem, vars::Vector{Union{Symbol,Expr}}, op::Symbol, bounds::Vector{Float64}) = 
    JuMPVarChoice(problem, vars, op, bounds)

function JuMPVarChoice(problem::JuMPProblem, varname::String, op::Symbol, bound::Float64)
    if op == :(<=)
        varind = findin(problem.model.colNames, varname)
        boundbefore = problem.model.colUpper[varind]
        return JuMPVarChoiceLeq(varind, boundbefore, bound)
    elseif op == :(>=)
        varind = findin(problem.model.colNames, varname)
        boundbefore = problem.model.colLower[varind]
        boundafter = bound
        return JuMPVarChoiceGeq(varind, boundbefore, bound)
    elseif op == :(=) || op == :(==)
        varind = findin(problem.model.colNames, varname)
        boundsbefore = (problem.model.colLower[varind], 
            problem.model.colUpper[varind])
        varcat = problem.model.colCat[varind]
        return JuMPVarChoiceEq(varind, varcat, boundsbefore, bound)
    end
end
function JuMPVarChoice(problem::JuMPProblem, varind::Int, op::Symbol, bound::Float64)
    if op == :(<=)
        boundbefore = problem.model.colUpper[varind]
        return JuMPVarChoiceLeq(varind, boundbefore, bound)
    elseif op == :(>=)
        boundbefore = problem.model.colLower[varind]
        return JuMPVarChoiceGeq(varind, boundbefore, bound)
    elseif op == :(=) || op == :(==)
        boundsbefore = (problem.model.colLower[varind], problem.model.colUpper[varind])
        varcat = problem.model.colCat[varind]
        return JuMPVarChoiceEq(varind, varcat, boundsbefore, bound)
    end
end
function JuMPVarChoice(problem::JuMPProblem, vars::Vector{Union{Symbol,Expr}}, op::Symbol, bounds::Vector{Float64})
    if op == :(<=)
        map((i) -> (JuMPVarChoiceLeq(vars[i], 
            problem.model.colUpper[problem.varindlookup[vars[i]]], bounds[i])), 
            1:length(vars))
    elseif op == :(>=)
        map((i) -> (JuMPVarChoiceGeq(vars[i], 
            problem.model.colLower[problem.varindlookup[vars[i]]], bounds[i])), 
            1:length(vars))
    elseif op == :(=) || op == :(==)
        subchoices = map(((i) -> (JuMPVarChoiceEq(vars[i], 
            problem.model.colCat[problem.varindlookup[vars[i]]],
            (problem.model.colLower[problem.varindlookup[vars[i]]], 
            problem.model.colUpper[problem.varindlookup[vars[i]]]), bounds[i]))), 
            1:length(vars))
    end
    return JuMPChoice(Vector{AbstractChoice{:JuMP}}(subchoices))
end
function JuMPLinearConstrChoice(problem::JuMPProblem, linconstr::JuMP.LinearConstraint)
    ind = length(problem.linconstr)+1
    return JuMPNewLinearConstrChoice(ind, linconstr)
end
function JuMPLinearConstrChoice(problem::JuMPProblem, newlinconstrs::Vector{JuMP.LinearConstraint})
    choice = JuMPChoice()
    offset = length(problem.linconstr)
    for i in 1:length(newlinconstrs)
        ind = offset+i
        append!(choice, JuMPNewLinearConstrChoice(ind, newlinconstrs[i]))
    end
    return choice
end
function JuMPLinearConstrChoice(problem::JuMPProblem, ind::Int, coeffsafter::Vector{Float64}, constantafter::Float64, boundsafter::Tuple{Float64,Float64})
    constr = problem.linconstr[ind]
    coeffsbefore = constr.terms.coeffs
    constantbefore = constr.terms.constant
    boundsbefore = (constr.lb, constr.ub)

    return JuMPLinearConstrChoice{:Old, :Static}(ind, coeffsbefore, constantbefore, 
        boundsbefore, coeffsafter, constantafter, boundsafter)
end
function JuMPLinearConstrChoice(problem::JuMPProblem, ind::Int, varsafter::Vector{Union{Symbol,Expr}}, coeffsafter::Vector{Float64}, constantafter::Float64, boundsafter::Tuple{Float64,Float64})
    constr = problem.linconstr[ind]
    varsbefore = map((i)->(parse(string(constr.terms.vars[i]))), 1:length(contr.terms.vars))
    coeffsbefore = constr.terms.coeffs
    constantbefore = constr.terms.constant
    boundsbefore = (constr.lb, constr.ub)

    return JuMPLinearConstrChoice{:Old, :Dynamic}(ind, varsbefore, coeffsbefore,
        constantbefore, boundsbefore, varsafter, coeffsafter, constantafter, boundsafter)
end
