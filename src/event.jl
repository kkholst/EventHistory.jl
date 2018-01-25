#=
Type 'EventClass' for Event History objects
Subtypes:
Surv:      (Right)-censored survival outcome
SurvInt:   Interval-censored outcome (Left,Right,Interval)
SurvTrunc: Right-censoring + left-truncation
CompRisk:  Competing risks model

constructor: Event
Access methods: Time, Status, Cause, ...
=#

const CensNot = 0
const CensLeft = 1
const CensRight = 2
const CensInt = 3

#=
Types/classes: Surv,CompRisk,...
=#

@compat abstract type EventClass end
#abstract EventInt <: EventClass

immutable Surv <: EventClass
    Time::Number  # Exit time
    Status::Bool  # Censoring status
end

immutable SurvTrunc <: EventClass
    Entry::Number # Entry time
    Time::Number  # Exit time
    Status::Bool  # Censoring status
end

immutable CompRisk <: EventClass
    Entry::Number # Entry time
    Time::Number  # Exit time
    Status::Bool  # Censoring status
    Cause::Number
    function CompRisk(Entry,Time,Cause) new(Entry,Time,Cause!=0,Cause) end
end

immutable SurvInt <: EventClass
    Time::Number     # Time 1
    Time2::Number    # Interval [Time;Time2]. Use Inf/-Inf for left/right censoring
    Status::Int      # Censoring (0:none,1:right,2:left,3:interval)
    function SurvInt(Time,Time2)
        if Time>Time2 error("time 1 larger than time 2") end
        if (Time==Time2)
            status = EventHistory.CensNot
        elseif Time2==Inf
            status = EventHistory.CensRight
        elseif Time==-Inf
            status = EventHistory.CensLeft
        else
            status = EventHistory.CensInt
        end
        new(Time,Time2,status)
    end
end

const Events = Vector{EventClass}
const Vec = Union{Vector,Matrix} 
const BoolVec = Union{Vector{Bool},Matrix{Bool},BitVector}


#=
show methods
=#

function show(io::IO, obj::EventClass)
    print(io, obj.Time, obj.Status>0 ? "":"+")
end

function show(io::IO, obj::SurvTrunc)
    print(io, "(", obj.Entry, ";", obj.Time, obj.Status>0 ? "":"+","]")
end

function show(io::IO, obj::CompRisk)
    print(io, "(", obj.Entry, ";", obj.Time, ":", obj.Status==0 ? "+" : obj.Cause,"]")
end

function show(io::IO, obj::SurvInt)
    if obj.Status==EventHistory.CensNot val=obj.Time
    elseif obj.Status==EventHistory.CensLeft val=string("(-Inf;",obj.Time2,"]")
    elseif obj.Status==EventHistory.CensRight val=string("[",Time,";Inf)")
    else val = string("[",obj.Time,";",obj.Time2,"]"); end
    print(io,val)
end


#=
Accessors
=#

# Meta-programming definitions of Time,Entry,Status,Cause access methods
for key = (:Time, :Entry, :Status, :Cause)
    for typ = (:Vector, :Matrix) #, :(DataFrames.DataVector))
        @eval function $(key){T<:EventClass}(e::$(typ){T})
            ##$(symbol(string("Event_",key)))(e)
            n = length(e)
            res = Array{typeof(e[1].$(key))}(n)
            for i=1:n
                res[i] = e[i].$(key)
            end
            res
        end
    end
end

export Time,Entry,Status,Cause


#=
Event constructor
=#

function Event(time::Vec,
               status::BoolVec)
    n = length(time)
    sz = size(time)
    if (length(sz)==1)
        E = Array{EventHistory.Surv}(n)
    else
        E = Array{EventHistory.Surv}(sz[1],sz[2])
    end
    for i=1:n
        E[i] = EventHistory.Surv(time[i],status[i])
        ## E[i] = EventHistory.Surv(time[i],convert(Bool,status[i]))
    end
    E
end

# function Event(time::Union{Vector,DataFrames.DataVector},
#                status::Union{Vector,DataFrames.DataVector},
#                method::Symbol=:comprisk)
#     Event(time,status,string(method))
# end

function Event(time::Vec,
               status::Vec,
               method::AbstractString)
    n = length(time)
    sz = size(time)
    ## interval
    if lowercase(method)=="interval"
        if (length(sz)==1)
            E = Array{EventHistory.SurvInt}(n)
        else
            E = Array{EventHistory.SurvInt}(sz[1],sz[2])
        end
        for i=1:n
            E[i] = EventHistory.SurvInt(time[i],status[i])
        end
        return E
    end
    ## comprisk:
    if (length(sz)==1)
        E = Array{EventHistory.CompRisk}(n)
    else
        E = Array{EventHistory.CompRisk}(sz[1],sz[2])
    end    
    for i=1:n
        E[i] = EventHistory.CompRisk(0,time[i],status[i])
    end
    E
end

function Event(entry::Vec,
               time::Vec,
               status::BoolVec)
    n = length(time)
    sz = size(time)
    if (length(sz)==1)
        E = Array{EventHistory.SurvTrunc}(n)
    else            
        E = Array{EventHistory.SurvTrunc}(sz[1],sz[2])
    end
    for i=1:n
        E[i] = EventHistory.SurvTrunc(entry[i],time[i],status[i])
    end
    E
end

function Event(entry::Vec,
               time::Vec,
               status::Vec)
    n = length(time)
    sz = size(time)
    # issurv = typeof(status[1])==Bool
    # if (issurv)
    #     E = Array{EventHistory.SurvTrunc}(sz[1],sz[2])
    #     for i=1:n
    #         E[i] = EventHistory.SurvTrunc(entry[i],time[i],status[i])
    #     end
    #     return E
    # end
    if (length(sz)==1)
        E = Array{EventHistory.CompRisk}(n)
    else            
        E = Array{EventHistory.CompRisk}(sz[1],sz[2])
    end
    for i=1:n
        E[i] = EventHistory.CompRisk(entry[i],time[i],status[i])
    end
    E
end


function Event(var::Vector{Symbol}, data::DataFrame, censdef::Function=x->x.>0)
    p = size(var,1)
    if p==2
        Status = data[var[2]]
    else
        Status = data[var[3]]
    end
    Status = censdef(Status)
    if p==2
        return Event(data[var[1]],Status)
    end
    return Event(data[var[1]],data[var[2]],Status)
end


#=
Misc
=#

# function transpose(x::EventHistory.EventClass)
#     x
# end
# export transpose

