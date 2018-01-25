using Base.Test, DataFrames, EventHistory, StatsModels

function test_show(x)
    io = IOBuffer()
    show(io, x)
end

d = DataFrame(start = [1,2,5,2,1,7,3,4,8,8],
              stop  = [2,3,6,7,8,9,9,9,14,17],
              event = [1,1,1,1,1,1,1,0,0,0],
              x     = [1,0,0,1,0,1,1,1,0,0]);
d[:S] = EventHistory.Event([:start,:stop,:event], d);
mm = EventHistory.phreg(@formula(S~x), d)


@test EventHistory.coef(mm) ≈ [-0.02110521] atol=1e-8
@test mm.invhess ≈ [0.6323063] atol=1e-6 

@test mm.grad[1]<1e-6

p = EventHistory.predict(mm, X=[0; 1], time=[6]; order=true);
@test isapprox(p[2:3], [0.3522254,0.3599854], atol=1e-6)

d = DataFrame(start = [1,2,5,2,1,7,3,4,8,8],
              stop  = [2,3,6,7,8,9,9.5,10,14,17],
              event = [1,1,1,1,1,1,1,0,0,0],
              x     = [1,0,0,1,0,1,1,1,0,0],
              z     = [1.0,0,2.0,0,3.0,0,4.0,0,5.0,0]);
d[:S] = EventHistory.Event([:start,:stop,:event], d);
mm = EventHistory.phreg(@formula(S~x+z), d)
@test length(coef(mm)) == 2
