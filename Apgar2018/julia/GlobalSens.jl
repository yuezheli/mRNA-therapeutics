# This script create global sensitivity analysis for the LNP-mRNA model

include("Apgar2018.jl")

using GlobalSensitivity, QuasiMonteCarlo
using DifferentialEquations.EnsembleAnalysis


# solve the system
sol = solve(prob, Tsit5(), reltol = 1e-2, progress=true, progress_steps=60, saveat = 300);
t = collect(sol.t);
default_mRNA_AUC = sum([sol.u[i][5] for i in 1:length(sol.t)]); 

# set up global sensitivity analysis
## focus global sensitivity analysis on ka, kl, dmRNA
lb = [1.17e-7, 1.93e-7, 1.07e-7];
ub = [1.17e-3, 1.93e-3, 1.07e-3];
N = 500 ; # sampling points
A,B = QuasiMonteCarlo.generate_design_matrices(N,lb,ub,SobolSample());

# define function for global sensitivity analysis
## use mRNA AUC for readout
f1 = function(p)
    # pass parameters into the ODE problem
    function prob_func(prob, i, repeat)
        # copy the default pars
        ApgarP = deepcopy(prob.p);
        ApgarP[4] = p[1, i]; # change ka value
        ApgarP[7] = p[2, i]; # change kl value
        ApgarP[8] = p[3, i]; # change dmRNA value
        prob2 = remake(prob,p=ApgarP);
        return prob2;
    end
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func);
    sol = solve(ensemble_prob, EnsembleThreads(); saveat=t, trajectories=size(p,2), reltol = 1e-2);
    
    out = zeros(1,size(p,2));
    for i in 1:size(p,2)
        # extract cytosolic mRNA mass
        cytomRNA = [sol[i].u[j][5] for j in 1:length(sol[i].t)];
        AUC_mRNA = sum(cytomRNA);
        out[1,i] = AUC_mRNA;
    end
    return out; 
end

@time sobol_result_new = gsa(f1,Sobol(),A,B, batch=true);

# plot result
using Plots

p1 = bar(["ka","kl", "dmRNA"],sobol_result_new.ST[1,:],title="Total Order AUC(mRNA)",legend=false);
p2 = bar(["ka","kl", "dmRNA"],sobol_result_new.S1[1,:],title="First Order AUC(mRNA)",legend=false);
plot(p1,p2)