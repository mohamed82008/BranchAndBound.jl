function commit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPChoice) where {Sense}
    for subchoice in choice.subchoices
        commit(problem, subchoice)
    end
end
function commit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPVarChoiceLeq) where {Sense}
    model(problem).colUpper[choice.varind] = choice.boundafter    
    return
end
function commit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPVarChoiceGeq) where {Sense}
    model(problem).colLower[choice.varind] = choice.boundafter
    return
end
function commit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPVarChoiceEq) where {Sense}
    _model = model(problem)
    _model.colVal[choice.varind] = choice.boundafter
    _model.colLower[choice.varind] = choice.boundafter
    _model.colUpper[choice.varind] = choice.boundafter
    _model.colCat[choice.varind] = :Fixed
    return
end
function commit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPNewLinearConstrChoice) where {Sense}
    push!(model(problem).linconstr, choice.constr)
    return
end
function commit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPStaticLinearConstrChoice) where {Sense}
    constr = model(problem).linconstr[choice.ind]
    constr.lb = choice.boundsafter[1]
    constr.ub = choice.boundsafter[2]
    constr.terms.coeffs = choice.coeffsafter
    constr.constant = choice.constantafter
    return
end
function commit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense},
    choice::JuMPDynamicLinearConstrChoice) where {Sense}
    constr = model(problem).linconstr[choice.ind]
    constr.lb = choice.boundsafter[1]
    constr.ub = choice.boundsafter[2]
    constr.terms.vars = eval.(choice.varsafter)
    constr.terms.coeffs = choice.coeffsafter
    constr.terms.constant = choice.constantafter
    return
end
function uncommit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPChoice) where {Sense}
    for subchoice in choice.subchoices
        uncommit(problem, subchoice)
    end
end
function uncommit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPVarChoiceLeq) where {Sense}
    model(problem).colUpper[choice.varind] = choice.boundbefore
    return
end
function uncommit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPVarChoiceGeq) where {Sense}
    model(problem).colLower[choice.varind] = choice.boundbefore
    return
end
function uncommit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPVarChoiceEq) where {Sense}
    _model = model(problem)
    _model.colLower[choice.varind] = choice.boundsbefore[1]
    _model.colUpper[choice.varind] = choice.boundsbefore[2]
    _model.colCat[choice.varind] = choice.varcat
    return
end
function uncommit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPNewLinearConstrChoice) where {Sense}
    constr = model(problem).linconstr[choice.ind]
    constr.lb = constr.ub = constr.terms.constant = 0.
    constr.terms.vars = JuMP.Variable[]
    constr.terms.coeffs = Float64[]
    return
end
function uncommit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPStaticLinearConstrChoice) where {Sense}
    constr = model(problem).linconstr[choice.ind]
    constr.lb = choice.boundsbefore[1]
    constr.ub = choice.boundsbefore[2]
    constr.terms.coeffs = choice.coeffsbefore
    constr.constant = choice.constantbefore
    return
end
function uncommit(problem::BBProblem{JuMPProblem,Vector{Float64},AbstractChoice{:JuMP},Sense}, 
    choice::JuMPDynamicLinearConstrChoice) where {Sense}
    _model = model(problem)
    constr = model(problem).linconstr[choice.ind]
    constr.lb = choice.boundsbefore[1]
    constr.ub = choice.boundsbefore[2]
    constr.terms.vars = JuMP.Variable.(_model, getindex.(problem.relaxation.x.varindlookup, choice.varsbefore))
    constr.terms.coeffs = choice.coeffsbefore
    constr.terms.constant = choice.constantbefore
    return
end

function commit(problem::BBProblem{T,S,C,Sense}, choice::C) where {T,S,C<:AbstractChoice,Sense}
    choice.commit(problem)
    return
end
function uncommit(problem::BBProblem{T,S,C,Sense}, choice::C) where {T,S,C<:AbstractChoice,Sense}
    choice.uncommit(problem)
    return
end
