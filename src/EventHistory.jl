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
using StatsBase: StatisticalModel, RegressionModel

type EventHistoryModel <: RegressionModel
    model::String
    call::Expr
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
    print(io,"\nModel: ", obj.model,",", obj.eventtype, " ", obj.call)
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
