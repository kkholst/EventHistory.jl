using DataFrames
using StatsModels
fm = @formula(Z ~ X + Y)
df = DataFrame(X = randn(10), Y = randn(10), Z = randn(10));
mf = ModelFrame(@formula(Z ~ X + Y), df);
