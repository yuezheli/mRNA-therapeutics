[PROB]

Model from Banks et al., 2003
https://www.nature.com/articles/3302076

Parameter used in this file is based on HeLa cells

[SET]

delta = 0.1, end = 10

[CMT]

M // plasmid in medium
C // plasmid in cytosol
N // plasmid in nucleus

[PARAM]
// all parameters in this section has unit h-1
k1 = 7.45e-4
k2 = 1.24
k3 = 6.68e-1
k4 = 4.5e-1

// volumes, all parameters in this section has unit L-1 (converted from uL)
Vc = 3.98E-7 // cell cytoplasm volume
Vm = 1e-3 // volume of medium
Vn = 7.22e-9 // cell nucleus medium

[MAIN]

// compute both diffusive and active rate in the model
double k_cn_diff = k4 * Vn; 
double k_cn_act = k3 * Vc/Vn - k4; 

[ODE]

dxdt_M = -k1*M + k2*C; 
dxdt_C = k1*M - k2*C - k3*C + k4*N;
dxdt_N = k3*C - k4*N;

[CAPTURE]
k_cn_diff, k_cn_act