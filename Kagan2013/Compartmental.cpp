[PROB]

Compartmental from Kagan et al., 2013
https://pubmed.ncbi.nlm.nih.gov/23793994/

[SET]

delta = 0.1, end = 24

[CMT]

// liposomal form
A_LIP_c
A_LIP_p

// nonliposomal form
A_c
A_p

// dummy variable
A_clearance

[PARAM]

FR = 1.83

krelease = 3.5e-3 // unit: h-1

wt = 0.25 // rate weight, unit: kg

dose = 0.8 // dose = 0.8 mg.kg-1

V_c= 8.52e-1 // unit: L.kg-1
k12 = 4.92e-1 // unit: h-1
k21 = 1.53e-1 // unit: h-1
kel = 1.6e-1  // unit: h-1

L_V_c = 7.07e-2 // unit: L.kg-1
L_k12 = 3.53e-1 // unit: h-1
L_k21 = 2.75e-1 // unit: h-1
L_kel = 8.65e-2 // unit: h-1

[MAIN]

double Vc = V_c * wt; // central volume, nonliposomal, unit: L-1
double L_Vc = L_V_c * wt; // central volume, liposomal, unit: L-1

// initial value of the drug
A_LIP_c_0 = dose * wt * (1-FR/100) ; 
A_c_0 = dose * wt * FR/100;

[ODE]

dxdt_A_LIP_c = - L_k12*A_LIP_c + L_k21*A_LIP_p - L_kel*A_LIP_c - krelease*A_LIP_c; 
dxdt_A_LIP_p = L_k12*A_LIP_c - L_k21*A_LIP_p - krelease*A_LIP_p; 

dxdt_A_c = krelease*A_LIP_c - kel*A_c - k12*A_c + k21*A_p;
dxdt_A_p = k12*A_c - k21*A_p + krelease*A_LIP_p;

dxdt_A_clearance = L_kel*A_LIP_c + kel*A_c;

[TABLE]

capture totaldrug = A_LIP_c + A_LIP_p + A_c + A_p + A_clearance;

capture C_LIP_AmB = A_LIP_c/ L_Vc;
capture C_nonlip_AmB = A_c/ Vc;
capture C_AmB = C_LIP_AmB + C_nonlip_AmB; 
