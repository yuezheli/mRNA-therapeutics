using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra, Catalyst;

@parameters t, kw, k12, k21, ka, ke, de, kl, dmRNA, kt, dUGTc, ktbg; 
@variables LNP(t), LNPp(t), LNPa(t), LNPe(t), mRNAc(t), protein(t); 

rx = []; 

# LNP clearance
push!(rx, Reaction(kw, [LNP], nothing)); # central compartment
push!(rx, Reaction(kw, [LNPp], nothing)); # peripheral compartment

# LNP distribution
push!(rx, Reaction(k12, [LNP], [LNPp])); # central to peripheral
push!(rx, Reaction(k21, [LNPp], [LNP])); # peripheral to central 

# LNP attachment to hepatocytes
push!(rx, Reaction(ka, [LNP], [LNPa]));

# cytosol, mRNA
push!(rx, Reaction(ke, [LNPa], [LNPe])); # mRNA-LNP internalization
push!(rx, Reaction(de, [LNPe], nothing)); # endosomal degredation
push!(rx, Reaction(kl, [LNPe], [mRNAc])); # endosomal escape
push!(rx, Reaction(dmRNA, [mRNAc], nothing)); # mRNA degredation

# cytosol, protein
push!(rx, Reaction(kt*mRNAc, nothing, [protein])); # trans-gene protein synthesis
push!(rx, Reaction(ktbg, nothing, [protein])); # endogenous protein synthesis
push!(rx, Reaction(dUGTc, [protein], nothing)); # protein degredation

@named rs = ReactionSystem(rx, t, [LNP, LNPp, LNPa, LNPe, mRNAc, protein], [kw, k12, k21, ka, ke, de, kl, dmRNA, kt, dUGTc, ktbg]); 

# define parameters 
pars = vcat(Pair(kw, 2.41e-5), Pair(k12, 4.79e-5), Pair(k21, 2.65e-7), Pair(ka, 1.17e-5), Pair(ke, 7.7e-5), Pair(de, 9.32e-5), 
            Pair(kl, 1.93e-5), Pair(dmRNA, 1.07e-5), Pair(kt, 17.73), Pair(dUGTc, 6.76e-6), Pair(ktbg, 1e-5)); 


# initial condition 
dosing = 0.3;  # LNP dose; 0.3mg.kg-1
animal_weight = 0.4; # assume Gunn rats weight 400g
moleweight_LNP = 1; # estimation based on https://pubs.acs.org/doi/10.1021/mp500367k; mg.nmol-1
dose = dosing * animal_weight / moleweight_LNP;

initLNP = [LNP => dose, LNPp => 0, LNPa => 0, LNPe => 0, mRNAc => 0, protein => 0]; 

# define time span
tspan = (0., 60*60*24*3); 

# define ODE
ode = convert(ODESystem, rs);

# check the system is fine
using Latexify; 
latexify(ode)

# solve the equation
prob = ODEProblem(ode, initLNP, tspan, pars);
sol = solve(prob, Tsit5(), reltol = 1e-2, progress=true, progress_steps=60, saveat = 300);

# read in simulation result from mrgsolve
using CSV, DataFrames;
mrgsolveresult = CSV.read("../data/ForJuliaValidation.csv", DataFrame);

# plot results
plotLNP = plot(sol.t/3600, [sol.u[i][1] for i in 1:length(sol.t)], label="LNP, central compartment",  linewidth = 2);
plot!(mrgsolveresult.hour, mrgsolveresult.LNP, label = "LNP, mrgsolve", linestyle = :dot, linewidth = 4);
xlims!(0.0, 72);
xticks!([0.0, 24, 48, 72], ["0", "24", "48", "72"]);
xlabel!("Time (hr)");
ylabel!("LNP in plasma (nmol)"); 


plotmRNA = plot(sol.t/3600, [sol.u[i][5] for i in 1:length(sol.t)], label="cytosolic mRNA",  linewidth = 2);
plot!(mrgsolveresult.hour, mrgsolveresult.mRNAc, label = "mRNAc, mrgsolve", linestyle = :dot, linewidth = 4);
xlims!(0.0, 72);
xticks!([0.0, 24, 48, 72], ["0", "24", "48", "72"]);
xlabel!("Time (hr)");
ylabel!("cytosolic mRNA (nmol)"); 


plotProtein = plot(sol.t/3600, [sol.u[i][6] for i in 1:length(sol.t)], label="protein",  linewidth = 2);
plot!(mrgsolveresult.hour, mrgsolveresult.protein, label = "protein, mrgsolve", linestyle = :dot, linewidth = 4, legend=:bottomright);
xlims!(0.0, 72);
xticks!([0.0, 24, 48, 72], ["0", "24", "48", "72"]);
xlabel!("Time (hr)");
ylabel!("cytosolic protein (nmol)"); 

plotd = plot(plotLNP, plotmRNA, plotProtein, layout = @layout [a b c]);
plot!(size = (1200, 400)); 

savefig(plotd,"../img/julia_Apgar2018.png");