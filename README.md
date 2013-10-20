# EventHistory

Event History Analysis for the Julia Language
  
## Event class

Right-censored event-time
```julia
> stop   = [2,3,3];
> status = [false,true,true];
> e1 = Event(stop,status)
3-element Array{Surv,1}:
 2+
 3 
 3 
> Time(e1)
3-element Array{Int64,1}:
 2
 3
 3
> Status(e1)
3-element Array{Bool,1}:
 false
  true
  true
```

Right-censored+left truncation
```julia
> start  = [0,1,2];
> e2 = Event(start,stop,status)
3-element Array{SurvTrunc,1}:
 (0;2+]
 (1;3] 
 (2;3] 
> Entry(e2)
3-element Array{Int64,1}:
 0
 1
 2
```

Competing risks model
```julia
> cause = [0,2,1];
> e3 = Event(start,stop,cause)
3-element Array{CompRisk,1}:
 (0;2:0]
 (1;3:+]
 (2;3:+]
> Status(e3)
3-element Array{Bool,1}:
  true
 false
 false
> Cause(e3)
3-element Array{Int64,1}:
 0
 2
 1
```

Interval censoring
```julia
> right  =  [2,3,Inf];
> left =  [1,-Inf,1];
> e4=Event(left,right,"interval")
3-element Array{SurvInt,1}:
 [1.0;2.0] 
 (-Inf;3.0]
 [Time;Inf)
```

Formula syntax (see also below)
```julia
> using DataFrames
> d = DataFrame(start=start,stop=stop,status=status);
> Event([:stop,:status],d)
3-element Array{Surv,1}:
 2+
 3 
 3 
> Event([:start,:stop,:status],d)
3-element Array{SurvTrunc,1}:
 (0;2+]
 (1;3] 
 (2;3] 
```

## Cox regression

## Examples

### Cox regression

Ovarian cancer example (randomized trial)
```julia
> using RDatasets
> using EventHistory
> ovarian = data("survival", "ovarian");
> ovarian["group"] = ovarian["rx"]-1;
> ovarian["S"] = Event([:futime,:fustat],ovarian);

> e = EventHistory.phreg(:(S~age+group),ovarian)

Model: Cox,Surv :(~(1))

n=26, events=12

2x4 DataFrame:
         Estimate      S.E.  dU^-1/2    P-value
[1,]     0.147327 0.0488846 0.046147 0.00258032
[2,]    -0.803973  0.633937 0.632049   0.204718


```

### Cox regression, Left truncation+right censoring
Simple example from the `survival` R-package
```julia
> d = DataFrame(start=[1,2,5,2,1,7,3,4,8,8],
                stop=[2,3,6,7,8,9,9,9,14,17],
                event=[1,1,1,1,1,1,1,0,0,0],
                x=[1,0,0,1,0,1,1,1,0,0])
> d["S"] = EventHistory.Event([:start,:stop,:event],d);
> e = EventHistory.phreg(:(S~x),d)

Model: Cox,Surv :(~(1))

n=10, events=7

1x4 DataFrame:
          Estimate     S.E.  dU^-1/2  P-value
[1,]    -0.0211052 0.838301 0.795177 0.979914

```

## Installation 

Get it from https://github.com/kkholst/EventHistory.jl
```julia
Pkg.add("EventHistory")
```

Note dependencies on: `DataFrames`,`Distributions`

## TODO

- Baseline estimates
- Stratified analysis
- Additive models
- Frailty models
- Residuals
- Ties
- ...
