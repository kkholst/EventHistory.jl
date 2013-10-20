# EventHistory.jl
# Event History Analysis in Julia
#
# Copyright (C) 2013   Klaus K. Holst
# 
##################################################

module EventHistory

using Distributions, DataFrames, EventHistory
import Base.show

include("event.jl")
include("aalen.jl") ## Aalens Additive Model
include("phreg.jl") ## Cox Proportional Hazards Model
include("optim.jl") ## Optimization routines

export aalen, phreg, Event

end # module
