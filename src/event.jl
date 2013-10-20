################################################################################
## Type 'EventClass' for Event History objects 
## Subtypes:
##   Surv:      (Right)-censored survival outcome
##   SurvInt:   Interval-censored outcome (Left,Right,Interval)
##   SurvTrunc: Right-censoring + left-truncation
##   CompRisk:  Competing risks model
##
## constructor: Event
## Access methods: Time, Status, Cause, ...
################################################################################

const CensNot = 0
const CensLeft = 1
const CensRight = 2
const CensInt = 3

###{{{ Types/classes: Surv,CompRisk,...

abstract EventClass
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
end
immutable SurvInt <: EventClass
    Time::Number     # Time 1
    Time2::Number    # Interval [Time;Time2]. Use Inf/-Inf for left/right censoring
    Status::Uint8 # Censoring (0:none,1:right,2:left,3:interval)
    function SurvInt(Time,Time2)
        if Time2 > Time error("Time 2 larger than time 1") end
        if (Time==Time2) 
            status = EventHistory.CensNot
        elseif Time2==Inf
            status = EventHistory.CensRight
        elseif Time==-Inf
            status = EventHistory.CensLeft
        else
            status = EventHistory.CensInt
            new(Time,Time2,status)
        end
    end
end
typealias Events Vector{EventClass}

###}}} Types

###{{{ show methods
function show(io::IO, obj::EventClass)
    print(io, obj.Time, obj.Status>0 ? "":"+")
end
function show(io::IO, obj::SurvTrunc)
    print(io, "(", obj.Entry, ";", obj.Time, obj.Status>0 ? "":"+","]")
end
function show(io::IO, obj::CompRisk)
    print(io, "(", obj.Entry, ";", obj.Time, ":", obj.Status==0 ? "+":obj.Status,"]")
end
function show(io::IO, obj::SurvInt)
    if obj.Status==EventHistory.CensNot val=obj.Time
    elseif obj.Status==EventHistory.CensLeft val=string("(-Inf;",obj.Time2,"]")
    elseif obj.Status==EventHistory.CensRight val=string("[",Time,";Inf)") 
        else val = string("[",Time,";",Time2,"]"); end
    print(io, "(", obj.Entry, ";", obj.Time, obj.Status>0 ? "":"+","]")
end
###}}} show methods

###{{{ Accessor
# Meta-programming definitions of Time,Entry,Status,Cause access methods
for key = (:Time, :Entry, :Status, :Cause)
    for typ = (:Vector, :DataVector) 
        @eval function $(key){T<:EventClass}(e::$(typ){T})
            ##$(symbol(string("Event_",key)))(e)
            n = size(e,1);
            res = Array(typeof(e[1].$(key)),n)
            for i=1:n
                res[i] = e[i].$(key)
            end
            res
        end
    end
end
      
export Time,Entry,Status,Cause
###}}} Accessor

###{{{ Event constructor
function Event(time::Union(Vector,DataVector),
               status::Union(Vector{Bool},DataVector{Bool}))
    n = size(time,1)    
    E = Array(EventHistory.Surv,n)
    for i=1:n
        E[i] = EventHistory.Surv(time[i],status[i])
        ## E[i] = EventHistory.Surv(time[i],convert(Bool,status[i]))
    end
    E
end

function Event(time::Union(Vector,DataVector),status::Union(Vector,DataVector))
    n = size(time,1)    
    E = Array(EventHistory.Surv,n)
    for i=1:n
        E[i] = EventHistory.CompRisk(0,time[i],status[i]==0,status[i])
    end
    E
end

function Event(entry::Union(Vector,DataVector), time::Union(Vector,DataVector), status::Union(Vector{Bool},DataVector{Bool}))
    n = size(time,1)
    E = Array(EventHistory.SurvTrunc,n)
    for i=1:n
        E[i] = EventHistory.SurvTrunc(entry[i],time[i],status[i])
    end
    E
end

function Event(entry::Union(Vector,DataVector), time::Union(Vector,DataVector), status::Union(Vector,DataVector))
    n = size(time,1)
    E = Array(EventHistory.CompRisk,n)
    for i=1:n
        E[i] = EventHistory.CompRisk(entry[i],time[i],status[i]==0,status[i])
    end
    E
end


function Event(var::Vector{Symbol}, data::DataFrame, censdef::Function=x->x.>0)
    p = size(var,1)
    if p==2
        Status = data[string(var[2])]
    else
        Status = data[string(var[3])]
    end
    if typeof(censdef)==Function
        Status = censdef(Status)
    end
    if p==2
        return Event(data[string(var[1])],Status)
    end
    return Event(data[string(var[1])],data[string(var[2])],Status)
end
###}}} Event constructor
