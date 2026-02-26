using ExaModels
N=5
M=10
c = ExaCore(backend = nothing)
x = variable(c, N, M)
objective(c, x[i] * x[i] + x[i] * x[i] for i = 2:N, j = 1:M)
dump(ExaModel(c; prod = true))
