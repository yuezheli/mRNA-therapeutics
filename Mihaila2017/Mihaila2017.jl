using DifferentialEquations, Plots, CSV, DataFrames, GlobalSensitivity, ComponentArrays;
using Parameters: @unpack; 

# define the model
function mihaila!(du, u, p, t)
    @unpack k1, k2, k3, k4, k5, k6, k7, K8, k9, E, R = p;
    @unpack N, S, SR, SRM, M = u; 

    # define some parameters
    Avogadro = 6.02e23 # Avogadro constant
    Vintra = 1.4e-12 # intracellular compartment volume; unit L
    k8 = (K8/ Avogadro) * 1e9/ Vintra; # unit converted from copies.h-1 -> nM.h-1

    # ODE on endosomal LNP (N)
    du.N = k1*E - k2*N - k3*N;

    # ODE on free siRNA (S)
    du.S = k2*N - k5*S - k4*S*R; 

    # ODE on Ago2-bound siRNA (SR)
    du.SR = k4*S*R - k6*M*SR;

    # ODE on active RISC mRNA (SRM)
    du.SRM = k6*M*SR - k7*SRM;

    # ODE on mRNA (M)
    du.M = k8 - k9*M - k7*SRM; 

    return nothing;
end

# convert molecule numbers to nM
function InitConv(MolCount = 10, vinterest = 1.4e-12, constantinterest = 6.02e23)
    # unit: molecule count; output: nM
    # input volume: unit L
    return(MolCount/ constantinterest/ vinterest * 1e9);
end


constP = ComponentArray(
    k1 = 0.005/3.6e-7, 
    k2 = 5e-4, k3 = 3, k4 = 0.001, k5 = 0.03, k6 = 0.1, k7 = 7.2, K8 = 100, k9 = 1, 
    E = 10, R = InitConv(1e4)); 


# define time frame for simulation
tspan = (0.0, 21); 

# define initial condition
u0 = ComponentArray(N = 0, S = 0, SR = 0, SRM = 0, M = InitConv(100));

# solve ODE
mihailaODE = ODEProblem(mihaila!, u0, tspan, constP);
mihailaResult = solve(mihailaODE, Tsit5(), saveat = 0.1);

# read data from mrgsolve
mrgsolveresult = CSV.read("data/JuliaValidation.csv", DataFrame);

# plot simulation result
plotmRNA = plot(mihailaResult.t, mihailaResult[:M], label = "mRNA"); 
plot!(mrgsolveresult[!, :time], mrgsolveresult[!, :M], label = "mRNA, mrgsolve", linestyle = :dot, linewidth = 4); 
xlabel!("Time (h)"); 
ylabel!("mRNA mass (nmol)"); 

plotLNP = plot(mihailaResult.t, mihailaResult[:N], label = "endosomal LNP"); 
plot!(mrgsolveresult[!, :time], mrgsolveresult[!, :N], label = "LNP, mrgsolve", linestyle = :dot, linewidth = 4, legend=:bottomright); 
xlabel!("Time (h)"); 
ylabel!("endosomal LNP mass (nmol)"); 

plotsiRNA = plot(mihailaResult.t, mihailaResult[:S], label = "free siRNA"); 
plot!(mrgsolveresult[!, :time], mrgsolveresult[!, :S], label = "LNP, mrgsolve", linestyle = :dot, linewidth = 4, legend=:bottomright); 
xlabel!("Time (h)"); 
ylabel!("free siRNA mass (nmol)"); 

plotd = plot(plotLNP, plotmRNA, plotsiRNA, layout = @layout [a b c]);

savefig(plotd,"img/julia_Mihaila2017.png");