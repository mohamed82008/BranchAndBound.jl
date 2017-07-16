import Base: push!, append!

sense(p::JuMPProblem) = p.model.objSense
function var(p::BBProblem{JuMPProblem, S, C, Sense}, var_sym, i::Int=0) where {S, C, Sense}
    if i == 0
        return p.relaxation.x.variables[var_sym]
    else
        return p.relaxation.x.variables[var_sym][i]
    end
end
function var(p::BBProblem{JuMPProblem, S, C, Sense}, var_sym, i::Tuple{Vararg{Int}}) where {S, C, Sense}
    if length(i) == 0
        return p.relaxation.x.variables[var_sym]
    else
        return p.relaxation.x.variables[var_sym][i...]
    end
end
var_sym(m::JuMPProblem, ind::Int) = Symbol(m.model.colNames[ind])
model(p::BBProblem{JuMPProblem, S, C, Sense}) where {S, C, Sense} = p.relaxation.x.model
tol(p::BBProblem{JuMPProblem, S, C, Sense}) where {S, C, Sense} = p.relaxation.x.tol
data(p::BBProblem) = p.relaxation.x.data

function push!(choice::JuMPChoice, subchoice::AbstractChoice{:JuMP})
    push!(choice.subchoices, subchoice)
    return
end
function append!(choice1::JuMPChoice, choice2::JuMPChoice)
    append!(choice1.subchoices, choice2.subchoices)
    return
end
