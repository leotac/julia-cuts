#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <ilcplex/ilocplex.h>
#define TOL 1e-6

using namespace std;

struct Inst
{
   int T;
   vector<double> c;
   vector<double> h;
   vector<double> K;
   vector<double> d;
};

Inst readULS(string path)
{
   ifstream file(path);
   Inst inst;
   file >> inst.T;
   inst.c = vector<double>(inst.T);
   for(int t=0; t<inst.T; ++t)
      file >> inst.c[t];
   inst.h = vector<double>(inst.T);
   for(int t=0; t<inst.T; ++t)
      file >> inst.h[t];
   inst.K = vector<double>(inst.T);
   for(int t=0; t<inst.T; ++t)
      file >> inst.K[t];
   inst.d = vector<double>(inst.T);
   for(int t=0; t<inst.T; ++t)
      file >> inst.d[t];
   file.close();
   return inst;
};


int calls = 0;
int separated = 0;
double separationtime = 0;

// User cut callback 
ILOUSERCUTCALLBACK3(lsCutCallback, vector<vector<double>>&, sumd, IloNumVarArray, z, IloNumVarArray, q)
{
   // Skip the separation if not at the end of the cut loop
   if ( !isAfterCutLoop() )
      return;

   auto t1 = std::chrono::high_resolution_clock::now();
   ++calls;
   IloEnv masterEnv = getEnv();

   int T = sumd.size();
   IloNumArray z_val = IloNumArray(masterEnv, T);
   IloNumArray q_val = IloNumArray(masterEnv, T);
   getValues(z_val, z);
   getValues(q_val, q);

   vector<int> S = vector<int>(T, 0);
   for(int l = 0; l < T; ++l)
   {
      fill(S.begin(), S.begin() + l, 0);

      double lhsvalue = 0;  //q(L\S) + sum{d[j:l]*z[j] for j in S}
      bool empty = true;
      for(int j = 0; j <= l; ++j)
      {
         if(q_val[j] > sumd[l][j]*z_val[j] + TOL)
         {
            S[j] = 1;
            empty = false;
            lhsvalue += sumd[l][j]*z_val[j];
         }
         else
            lhsvalue += q_val[j];
      }

      if(empty)
         continue;
      else
      {
         if(lhsvalue < sumd[l][0] - TOL)
         {
            IloExpr lhs = IloExpr(masterEnv);
            for(int j = 0; j <= l; ++j)
            {
               if(S[j] == 1)
                  lhs += sumd[l][j]*z[j];
               else
                  lhs += q[j];
            }
            add(lhs >= sumd[l][0], IloCplex::UseCutPurge); //UseCutFilter, UseCutForce
            //cout << "Cut " << separated << ": " << (lhs >= sumd[l][0]) << endl;
            ++separated;
            lhs.end();      
            abortCutLoop();
            break;
         }
      }
   }   
   z_val.end();
   q_val.end();
   auto t2 = chrono::high_resolution_clock::now();
   separationtime += chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
   return;
};


int main(int argc, char** argv)
{
   if(argc != 2)
   {  cout << "Pass .dat file as argument." << endl;
      return 1;
   }
   string path = argv[1];
   Inst inst = readULS(path);
   int T = inst.T;
   vector<double>& c = inst.c;
   vector<double>& h = inst.h;
   vector<double>& K = inst.K;
   vector<double>& d = inst.d;

   IloEnv env;
   IloModel model = IloModel(env);
   IloObjective objective = IloMinimize(env);
 
   IloNumVarArray z = IloNumVarArray(env);
   IloNumVarArray q = IloNumVarArray(env);
   IloNumVarArray s = IloNumVarArray(env);
  
   /* Add variables */
   char var_name[255];
   for (int t = 0; t < T; ++t)
   {
      IloNumVar var;
      snprintf(var_name, 255, "z_%d", t+1);
      if (false)
         var = IloNumVar(env, 0, 1, IloNumVar::Float, var_name);
      else
         var = IloNumVar(env, 0, 1, IloNumVar::Bool, var_name);
      objective.setLinearCoef(var, K[t]);
      z.add(var);
   
      snprintf(var_name, 255, "q_%d", t+1);
      var = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, var_name);
      objective.setLinearCoef(var, c[t]);
      q.add(var);
  
      snprintf(var_name, 255, "s_%d", t+1);
      var = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, var_name);
      objective.setLinearCoef(var, h[t]);
      s.add(var);
   
      // Activation 
      IloRange con = IloRange(env, -IloInfinity, 0);
      con.setLinearCoef(q[t], 1.0);
      double subd = 0;
      for(int i = t; i < T; ++i)
         subd += d[i];
      con.setLinearCoef(z[t], -subd);
      model.add(con);
     
      // Balance equations
      con = IloRange(env, d[t], d[t]);
      con.setLinearCoef(q[t], 1.0);
      if(t>0)
         con.setLinearCoef(s[t-1], 1.0);
      con.setLinearCoef(s[t], -1.0);
      model.add(con);

   }

   model.add(objective);

   //precompute sumd[l][j] := sum(d[j:l])
   vector<vector<double>> sumd = vector<vector<double>>(T);   
   for(int l = 0; l < T; ++l)
   {
      sumd[l] = vector<double>(l + 1);
      sumd[l][l] = d[l];
      for(int j = l-1; j >= 0; --j)
         sumd[l][j] = sumd[l][j+1] + d[j];
   }


   IloCplex cplex = IloCplex(model);
   cplex.use(lsCutCallback(env, sumd, z, q));
   cplex.setParam(IloCplex::CutsFactor, 1.0);
   cplex.setParam(IloCplex::MIPDisplay, 2);
   cplex.solve();

   cout << "Separation time: "<< setprecision(2) << separationtime << " milliseconds" << endl;
   cout << "Separated: " << separated << endl;
   cout << "Called: " << calls << endl;
   cout << "Nodes: " << cplex.getNnodes() << endl;
   env.end();
   return 0;
}

