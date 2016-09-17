using Base.Test, DataFrames, EventHistory

function test_show(x)
    io = IOBuffer()
    show(io, x)
end


d = DataFrame(start=[1,2,5,2,1,7,3,4,8,8],
              stop=[2,3,6,7,8,9,9,9,14,17],
              event=[1,1,1,1,1,1,1,0,0,0],
              x=[1,0,0,1,0,1,1,1,0,0]);
d[:S] = EventHistory.Event([:start,:stop,:event],d);
mm = EventHistory.phreg(S~x,d)


@test_approx_eq_eps EventHistory.coef(mm) [-0.02110521] 1e-8
@test_approx_eq_eps mm.invhess 0.6323063 1e-6 

@test mm.grad[1]<1e-6

p = EventHistory.predict(mm,X=[0; 1], time=[6]; order=true);
@test_approx_eq_eps p[2:3] [0.3522254,0.3599854] 1e-6
