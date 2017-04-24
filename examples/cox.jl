import EventHistory
##reload("EventHistory")

##################################################
### EventClass
##################################################

## Right-censored event-time
stop   = [2,3,3];
status = [false,true,true];
e1     = EventHistory.Event(stop,status)

EventHistory.Time(e1)
EventHistory.Status(e1)

## Right-censored+left truncation
start = [0,1,2];
e2    = EventHistory.Event(start,stop,status)

EventHistory.Entry(e2)

## Competing risks model
cause = [0,2,1];
e3    = EventHistory.Event(start,stop,cause)

EventHistory.Status(e3)
EventHistory.Cause(e3)

## Interval censoring
right = [2,3,Inf];
left  = [1,-Inf,1];
e4    = EventHistory.Event(left,right,"interval")

## Formula syntax (see also below)
d = DataFrames.DataFrame(start=start,stop=stop,status=status);
EventHistory.Event([:stop,:status],d)
EventHistory.Event([:start,:stop,:status],d)


##################################################
### Cox regression
##################################################
using RDatasets
ovarian = dataset("survival", "ovarian");
ovarian[:Group] = ovarian[:Rx].-1;
S = EventHistory.Event([:FUTime,:FUStat], ovarian);
ovarian[:S] = EventHistory.Event([:FUTime,:FUStat], ovarian);

mm = EventHistory.phreg(@formula(S~Age+Group), ovarian)

## Prediction
EventHistory.predict(mm, surv=false, X=[0 0]); ## Baseline
s56 = EventHistory.predict(mm, X=[56 1], order=true); ## Survival probabilities age 40, group 1
EventHistory.predict(mm, X=[56 0], time=[100,400,600]); ## ... at time 100,400,600
EventHistory.predict(mm, X=[56 1; 56 0], time=[600,100,400]) ## ... both groups

using Winston
plot(s56[:,1], s56[:,2])

s = EventHistory.predict(mm, X=[56 1; 56 0], order=true)
pr = DataFrame(Time=[s[:,1];s[:,1]], S=[s[:,2];s[:,3]], Group=repeat(["Group1","Group2"], inner=size(s,1)))

using Gadfly
p = Gadfly.plot(pr, x="Time", y="S", color="Group",
                Geom.step, Geom.point,
                Guide.ylabel("Survival probability"), Guide.title("Age 56"))

draw(PNG("surv.png",7inch,7inch),p)

##################################################
### Cox regression with left truncation
##################################################
d = DataFrame(start=[1,2,5,2,1,7,3,4,8,8],
              stop=[2,3,6,7,8,9,9.5,10,14,17],
              event=[1,1,1,1,1,1,1,0,0,0],
              x=[1,0,0,1,0,1,1,1,0,0],
              z=[1.0,0,2.0,0,3.0,0,4.0,0,5.0,0]);
d[:S] = EventHistory.Event([:start,:stop,:event], d);
mm = EventHistory.phreg(@formula(S~x), d)

mm = EventHistory.phreg(@formula(S~x+z), d)




