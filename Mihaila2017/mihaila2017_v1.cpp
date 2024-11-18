[PROB]

Model from Mihaila et al., 2017
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415968/


[SET]

[CMT]

E // extracellular LNP
N // endosomal LNP
S // free siRNA
R // RISC complexes
SR // Ago2-bound siRNA
SRM // active RISC mRNA
M // mRNA

[PARAM]

k1 = 0.005 // LNP crossing the plasma membrane; L.h-1
k2 = 5e-4 // endosomal escape/ unpacking; L.h-1
k3 = 3 // lysosomal degredation; L.h-1
k4 = 0.001 // siRNA loading to RISC; L.nM-1.h-1
k5 = 0.03 // degredation of siRNA in the cytoplasma; L.h-1
k6 = 0.1 // formation of active RISC with target mRNA; L.nM-1.h-1
k7 = 7.2 // cleavage of target mRNA by RISC; L.h-1
K8 = 100 // transcription of mRNA; copies.h-1
k9 = 1 // degredation of mRNA; L.h-1

Vextra = 3e-4 // extracellular compartment volume; unit L
Vintra = 1.4e-12 // intracellular compartment volume; unit L

// constant conversion
Avogadro = 6.02e23


[ODE]

dxdt_E = 0;
dxdt_N = k1*E - k2*N - k3*N;
dxdt_S = k2*N - k5*S - k4*S*R; 
dxdt_SR = k4*S*R - k6*M*SR;
dxdt_SRM = k6*M*SR - k7*SRM;
dxdt_M = K8 - k9*M - k7*SRM; 
