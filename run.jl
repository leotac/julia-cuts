include("ls.jl")
using Gurobi, CPLEX
ls.solveULS("test.dat",solver=GurobiSolver(Cuts=0),valid=true)
ls.solveULS("test.dat",solver=CplexSolver(CPX_PARAM_CUTSFACTOR=1),valid=true)

