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

K2 = 5e-4 // endosomal escape/ unpacking; L.h-1
K3 = 3 // lysosomal degredation; L.h-1
K4 = 0.001 // siRNA loading to RISC; L.nM-1.h-1
K5 = 0.03 // degredation of siRNA in the cytoplasma; L.h-1
K6 = 0.1 // formation of active RISC with target mRNA; L.nM-1.h-1
K7 = 7.2 // cleavage of target mRNA by RISC; L.h-1
K8 = 100 // transcription of mRNA; copies.h-1
K9 = 1 // degredation of mRNA; L.h-1

Vextra = 3e-4 // extracellular compartment volume; unit L
Vintra = 1.4e-12 // intracellular compartment volume; unit L

// constant conversion
Avogadro = 6.02e23

[MAIN]
// convert all unit related to copies to L.h-1
double k8 = (K8/ Avogadro) * 1e9/ Vintra; // unit converted from copies.h-1 -> nM.h-1

// convert all unit to h-1
double k2 = K2/ Vintra; 
double k3 = K3/ Vintra; 
double k4 = K4/ Vintra; 
double k5 = K5/ Vintra; 
double k6 = K6/ Vintra; 
double k7 = K7/ Vintra; 
double k9 = K9/ Vintra; 

[ODE]

dxdt_E = 0;
dxdt_N = k1*E  - k2*N - k3*N;
dxdt_S = k2*N - k5*S - k4*S*R; 
dxdt_SR = k4*S*R - k6*M*SR;
dxdt_SRM = k6*M*SR - k7*SRM;
dxdt_M = k8 - k9*M - k7*SRM; 

[TABLE]

capture RNAcount = M * 1e-9 * Vintra * Avogadro; 