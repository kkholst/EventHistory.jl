# EventHistory.jl
# Event History Analysis in Julia
#
# Copyright (C) 2013   Klaus K. Holst
#
# Authors: Klaus Kähler Holst <kkho@biostat.ku.dk>
# Keywords: event history analysis, survival analysis,
#           cox regression, competing risks models
##################################################

module EventHistory

using Distributions, DataFrames
import Base.show
import DataFrames: ModelFrame, ModelMatrix
import StatsBase: coef, coeftable, confint, vcov, predict
import Calculus: deparse
using StatsBase: StatisticalModel, RegressionModel

type EventHistoryModel <: RegressionModel
    model::String
    formula::Formula
    eventtype::DataType
    coef::Vector{Float64}
    coefmat::DataFrame
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

export phreg, Event, predict ## lifetable, aalen

end # module
