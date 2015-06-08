module ls
using JuMP, CPLEX, Gurobi, MathProgBase
include("utils.jl")

const TOL = 1e-6

function separate(T::Int64, sumd::Array{Float64, 2}, z_val, q_val, z, q)
   S = zeros(Bool, T)
   for l in 1:T
      fill!(S, false)
      lhsvalue = 0.  #q(L\S) + sum{d[j:l]*z[j] for j in S}
      empty = true
      for j in 1:l
         if q_val[j] > sumd[j,l]*z_val[j] + TOL
            S[j] = true
            empty = false
            lhsvalue += sumd[j,l]*z_val[j]
         else
            lhsvalue += q_val[j]
         end
      end
      if empty
         continue
      end
      if lhsvalue < sumd[1,l] - TOL
         lhs = sum(q[1:l])
         for j = (1:T)[S]
            lhs += sumd[j,l]*z[j] - q[j]
         end
         return lhs - sumd[1,l]
      end
   end
   return
end

function solveULS(path::String; relax=false, solver=CplexSolver(), valid::Bool = true) 
   T, c, h, K, d = readULS(path)
   m = Model(solver = solver)
   @defVar(m, z[1:T], Bin)
   @defVar(m, q[i = 1:T] >= 0)
   @defVar(m, s[1:T] >= 0)

   @setObjective(m, Min, sum{c[t]*q[t] + K[t]*z[t] + h[t]*s[t], t=1:T})
   @addConstraint(m, activation[t = 1:T], q[t] <= sum(d[t:T])*z[t])
   @addConstraint(m, balance[t = 1:T], (t>1?s[t-1]:0) + q[t] == s[t] + d[t])

   #precompute sum(d[j:l])
   sumd = zeros(Float64, T, T)
   for l = 1:T, j = 1:l
      sumd[j,l] = sum(d[j:l])
   end

   separationtime = 0.
   separations = 0
   called = 0
   function lSgenerator(cb)
      called += 1
      tt = time()
      z_val = getValue(z)
      q_val = getValue(q)
      expr = separate(T, sumd, z_val, q_val, z, q)
      if expr != nothing
         @addUserCut(cb, expr >= 0)
      end
      separationtime += time()-tt
      separations += 1
   end 

   if valid
      addCutCallback(m, lSgenerator)
   end

   status = solve(m)
   @printf("Objective value: %.2f\n", getObjectiveValue(m))
   @printf("Separation time: %.2f ms\n", separationtime*1000)
   println("Separated: $separations")
   try
      println("Nodes: ", CPLEX.getnodecnt(m.internalModel))
   catch y
      println("Nodes: ", Gurobi.get_node_count(MathProgBase.getrawsolver(m.internalModel)))
   end
   status
end
end
