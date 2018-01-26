import EventHistory
using DataFrames
using StatsModels
##using EventHistory
##reload("EventHistory")

d = DataFrame(start=[1,2,5,2,1,7,3,4,8,8],
              stop=[2,3,6,7,8,9,9.5,10,14,17],
              event=[1,1,1,1,1,1,1,0,0,0],
              x=[1,0,0,1,0,1,1,1,0,0],
              z=[1.0,0,2.0,0,3.0,0,4.0,0,5.0,0]);
d[:S] = EventHistory.Event([:start,:stop,:event], d);
mm = EventHistory.phreg(@formula(S~x), d)
