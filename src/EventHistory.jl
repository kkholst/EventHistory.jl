#=
 EventHistory.jl
 Event History Analysis in Julia

 Copyright (C) 2013-2018   Klaus K. Holst

 Authors: Klaus Kähler Holst <klaus@holst.it>
 Keywords: event history analysis, survival analysis,
           cox regression, competing risks models
=#

module EventHistory

using Distributions, DataFrames, StatsModels
import Base.show
import StatsBase: coef, coeftable, CoefTable, coefnames, confint, vcov, predict
import Calculus: deparse
import LinearAlgebra: pinv, diag
import StatsModels: Formula, ModelFrame, ModelMatrix
using StatsBase: StatisticalModel, RegressionModel
# using Compat
# import Compat.String

mutable struct EventHistoryModel <: RegressionModel
    model::AbstractString
    formula::Formula
    eventtype::DataType
    coef::Vector{Float64}
    coefmat::CoefTable
    IC::Matrix{Float64}
    invhess::Matrix{Float64}
    vcov::Matrix{Float64}
    opt::Vector
    grad::Vector{Float64}
    X::Matrix{Float64}
    eventtime::Matrix
    chaz::Matrix{Float64}
end

function show(io::IO, obj::EventHistoryModel)
    ## outcome = deparse(obj.call.args[2])
    ## covars = deparse(obj.call.args[3])
    ## print(io,"\nModel: ", obj.model,",", obj.eventtype, "; ",
    ##       string("$outcome ~ $covars"))
    print(io,"\nModel: ", obj.model, "; ", obj.formula,"\n")
    n = size(obj.eventtime,1)
    events::Int = sum(obj.eventtime[:,2])
    print(io,"\nn=",n,", events=",events,"\n\n")
    print(io,obj.coefmat)
end

include("event.jl")
include("phreg.jl") ## Cox Proportional Hazards Model
#include("lifetable.jl") ## Life-table construction
#include("aalen.jl") ## Aalens Additive Model
include("optim.jl") ## Optimization routines

export phreg, Event, predict, coef, vcov, coeftable,
    Time, Status, Cause, Entry
    ## lifetable, aalen

end # module
