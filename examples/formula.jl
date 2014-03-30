
fm = Z ~ X + Y
df = DataFrame(X = randn(10), Y = randn(10), Z = randn(10));
mf = ModelFrame(Z ~ X + Y, df);
