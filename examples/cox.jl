##################################################
### Cox regression
##################################################

using RDatasets
using EventHistory
ovarian = data("survival", "ovarian");
ovarian["group"] = ovarian["rx"]-1;
ovarian["S"] = Event([:futime,:fustat],ovarian);

e = EventHistory.phreg(:(S~age+group),ovarian)

##################################################
### Cox regression with left truncation
##################################################
d = DataFrame(start=[1,2,5,2,1,7,3,4,8,8],
              stop=[2,3,6,7,8,9,9,9,14,17],
              event=[1,1,1,1,1,1,1,0,0,0],
              x=[1,0,0,1,0,1,1,1,0,0])
d["S"] = EventHistory.Event([:start,:stop,:event],d);
e = EventHistory.phreg(:(S~x),d)
