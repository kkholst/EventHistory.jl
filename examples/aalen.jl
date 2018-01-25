#=
Aalens additive model
=#


## https://github.com/JuliaData/DataFrames.jl
## http://juliastats.github.io/StatsBase.jl/stable/
 

using EventHistory
using DataFrames
using DataFramesMeta
using StatsModels
using CSV
using RCall

reval("data('sTRACE',package='timereg')");
reval("write.csv(sTRACE,file='/data/trace.csv',row.names=FALSE)");
dd = CSV.read("/data/trace.csv");
dd[:S] = Event([:time,:status], dd, x->(x.>1));

M = ModelFrame(@formula(S~1+age+sex),dd);
X = ModelMatrix(M).m
S = M.df[:,1] ## convert(Vector,M.df[:,1])
time = EventHistory.Time(S)
status = EventHistory.Status(S)

function aalen(X, t, status; entry=[], id=[], beta=[], opt...)
    
end

