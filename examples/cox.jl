using EventHistory

##################################################
### EventClass
##################################################

## Right-censored event-time
stop   = [2,3,3];
status = [false,true,true];
e1 = Event(stop,status)

Time(e1)
Status(e1)

## Right-censored+left truncation
start  = [0,1,2];
e2 = Event(start,stop,status)

Entry(e2)

## Competing risks model
cause = [0,2,1];
e3 = Event(start,stop,cause)

Status(e3)
Cause(e3)

## Interval censoring
right  =  [2,3,Inf];
left =  [1,-Inf,1];
e4=Event(left,right,"interval")

## Formula syntax (see also below)
using DataFrames
d = DataFrame(start=start,stop=stop,status=status);
Event([:stop,:status],d)
Event([:start,:stop,:status],d)


##################################################
### Cox regression
##################################################

using RDatasets
ovarian = data("survival", "ovarian");
ovarian["group"] = ovarian["rx"]-1;
ovarian["S"] = Event([:futime,:fustat],ovarian);

e = phreg(:(S~age+group),ovarian)

##################################################
### Cox regression with left truncation
##################################################
d = DataFrame(start=[1,2,5,2,1,7,3,4,8,8],
              stop=[2,3,6,7,8,9,9,9,14,17],
              event=[1,1,1,1,1,1,1,0,0,0],
              x=[1,0,0,1,0,1,1,1,0,0]);
d["S"] = Event([:start,:stop,:event],d);
e = phreg(:(S~x),d)
