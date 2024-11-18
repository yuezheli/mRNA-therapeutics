[PROB]

Simplified from Kagan et al., 2013
https://pubmed.ncbi.nlm.nih.gov/23793994/

Only contain the PBPK model for AmB in the nonliposomal form 

[SET]

delta = 1, end = 24

[CMT]
//---- Nonliposomal compartments ----//
// plasma
A_pl

// GI
A_gi

// heart
A_ht

//spleen
A_sp_vas
A_sp_exv
A_sp_deep

// liver
A_li

// kidney
A_kd_vas
A_kd_exv
A_kd_deep

// lung
A_lu

// remainder
A_rm_vas
A_rm_exv

// drug clearance; dummy variable
A_clear

[PARAM]

//---- cardiac output ----//
// Qco = 53/1000 * 24 // unit: L.h-1; value from De Buck et al., 2007; https://dmd.aspetjournals.org/content/35/10/1766.short
Qco = (100-43.7)/100 * 14.1 * 0.25^0.75; // unit: L.h-1, rat HCT value from https://www.mdpi.com/1422-0067/19/9/2824/htm TABLE S1

f_u_pl = 0.11 // fraction of unbound in plasma, rats

// fraction of flows compared to cardiac output, rat values
q_frac_li = 18.3/100 // the volume that leaves liver
q_frac_kd = 14.1/100
q_frac_sp = 1/100
q_frac_gi = 14.3/100
q_frac_ht = 4.9/100
// q_frac_lu = 100/100

//---- Tissue volume ----//

wt = 0.25 // assuming rat weight = 0.25kg 

// note tissue volume is a percentage of body weight
// assume the volume unit should be L
// the calculation matches rat tissue volume from https://dmd.aspetjournals.org/content/35/10/1766
v_frac_li = 3.66/100
v_frac_kd = 0.73/100
v_frac_sp = 0.2/100
v_frac_gi = 2.7/100
v_frac_ht = 0.33/100
v_frac_lu = 0.5/100
v_frac_blood = 7.4/100 // value from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243854/

//---- fraction of vascular space ----//
vs_frac_li = 21/100
vs_frac_kd = 16/100
vs_frac_sp = 22/100
vs_frac_gi = 19/100
vs_frac_ht = 26/100
vs_frac_lu = 36/100
vs_frac_mu = 4/100
vs_frac_rm = 4/100


//---- nonliposomal AmB ----//
Kpgi  = 10.7
Kplu = 34.3
Kpht = 2
Kpli = 33

CL_li = 6e-2 // unit: L.h-1
f_u_kd = 7.26e-3
PSkd = 7.5e-2 // unit: L.h-1
Kakd = 2.58e-1 //unit: h-1
Kdkd = 1.29e-3 // unit: h-1
CL_kd = 1e-1
f_u_sp = 1.71e-3
PSsp = 5.98e-1 // unit: L.h-1
Kasp = 9.05e-1 // unit: h-1
Kdsp = 4.56e-3 // unit: h-1
f_u_rm = 1.6e-2
PSrm = 5.47e-1 // unit: L.h-1
CL_rm = 1.58e-1 // unit: L.h-1


// dosing amount
dose = 20 // unit: mg.kg-1; rat data

[MAIN]

// initial value of the drug
A_pl_0 = dose * wt ; 


// flow through organs
// double Qha = Qco * q_frac_ha; 
double Qli = Qco * q_frac_li;
double Qkd = Qco * q_frac_kd;
double Qsp = Qco * q_frac_sp;
double Qgi = Qco * q_frac_gi;
double Qht = Qco * q_frac_ht;

// double Qli = Qha + Qsp + Qgi; 
double Qha = Qli - Qsp - Qgi;

// double Qrm = Qco * (1 - q_frac_ha - q_frac_kd - q_frac_sp - q_frac_gi - q_frac_ht); 
double Qrm = Qco * (1 - q_frac_li - q_frac_kd - q_frac_ht);

//---- tissue volume ----/

double v_frac_rm = 1 - v_frac_blood - v_frac_li - v_frac_kd - v_frac_sp - v_frac_gi - v_frac_ht - v_frac_lu; 

double V_pl = wt * v_frac_blood; 
double V_li = wt * v_frac_li;
double V_kd = wt * v_frac_kd;
double V_sp = wt * v_frac_sp;
double V_gi = wt * v_frac_gi;
double V_ht = wt * v_frac_ht;
double V_lu = wt * v_frac_lu;
double V_rm = wt * v_frac_rm; 


//liver
double V_li_vas = wt * v_frac_li * vs_frac_li;
double V_li_exv = wt * v_frac_li * (1 - vs_frac_li);

// kidney
double V_kd_vas = wt * v_frac_kd * vs_frac_kd;
double V_kd_exv = wt * v_frac_kd * (1 - vs_frac_kd);

// spleen
double V_sp_vas = wt * v_frac_sp * vs_frac_sp;
double V_sp_exv = wt * v_frac_sp * (1 - vs_frac_sp);

// GI
double V_gi_vas = wt * v_frac_gi * vs_frac_gi;
double V_gi_exv = wt * v_frac_gi * (1 - vs_frac_gi);

// heart
double V_ht_vas = wt * v_frac_ht * vs_frac_ht;
double V_ht_exv = wt * v_frac_ht * (1 - vs_frac_ht);

// lung
double V_lu_vas = wt * v_frac_lu * vs_frac_lu;
double V_lu_exv = wt * v_frac_lu * (1 - vs_frac_lu);

// rest of the body
double V_rm_vas = wt * v_frac_rm * vs_frac_rm;
double V_rm_exv = wt * v_frac_rm * (1 - vs_frac_rm);

[ODE]
// list of nonliposomal concentrations
double C_pl = A_pl/ V_pl;
double C_gi = A_gi/ V_gi;
double C_ht = A_ht/ V_ht; 
double C_sp_vas = A_sp_vas/ V_sp_vas;
double C_sp_exv = A_sp_exv/ V_sp_exv;
double C_li = A_li/ V_li;
double C_kd_vas = A_kd_vas/ V_kd_vas;
double C_kd_exv = A_kd_exv/ V_kd_exv;
double C_lu = A_lu/ V_lu;
double C_rm_vas = A_rm_vas/ V_rm_vas; 
double C_rm_exv = A_rm_exv/ V_rm_exv;


//---- nonliposomal compartment ----//

// plasma
dxdt_A_pl = Qco*C_lu/Kplu - Qco*C_pl ; 

// GI
dxdt_A_gi = Qgi*C_pl - Qgi*C_gi/Kpgi  ; 

// heart
dxdt_A_ht = Qht*C_pl - Qht*C_ht/Kpht ;

// spleen
dxdt_A_sp_vas = Qsp*C_pl  - Qsp*C_sp_vas - PSsp*(f_u_pl * C_sp_vas - f_u_sp * C_sp_exv); 
dxdt_A_sp_exv = PSsp*(f_u_pl * C_sp_vas - f_u_sp * C_sp_exv) - Kasp*f_u_sp*C_sp_exv*V_sp_exv  + Kdsp*A_sp_deep ;
dxdt_A_sp_deep = Kasp*f_u_sp*C_sp_exv*V_sp_exv - Kdsp*A_sp_deep;

// kidney
dxdt_A_kd_vas = Qkd*C_pl - Qkd*C_kd_vas - PSkd*(f_u_pl * C_kd_vas - f_u_kd * C_kd_exv) - CL_kd*f_u_pl*C_kd_vas ;
dxdt_A_kd_exv = PSkd*(f_u_pl * C_kd_vas - f_u_kd * C_kd_exv) - Kakd*f_u_kd*C_kd_exv*V_kd_exv + Kdkd*A_kd_deep ;
dxdt_A_kd_deep = Kakd*f_u_kd*C_kd_exv*V_kd_exv - Kdkd*A_kd_deep;


// liver
dxdt_A_li = Qha*C_pl + Qsp*C_sp_vas + Qgi*C_gi/Kpgi - Qli*C_li/Kpli - CL_li*f_u_pl*C_li/Kpli ; 

// lung
dxdt_A_lu = Qli*C_li/Kpli + Qht*C_ht/Kpht + Qkd*C_kd_vas + Qrm*C_rm_vas - Qco*C_lu/Kplu ;

// remainder
dxdt_A_rm_vas = Qrm*C_pl - Qrm*C_rm_vas - PSrm*(f_u_pl*C_rm_vas - f_u_rm*C_rm_exv) ;
dxdt_A_rm_exv = PSrm*(f_u_pl*C_rm_vas - f_u_rm*C_rm_exv) - CL_rm*f_u_rm*C_rm_exv ;

// clearance; dummy compartment; added to track total drug mass
dxdt_A_clear = CL_li*f_u_pl*C_li/Kpli + CL_kd*f_u_pl*C_kd_vas + CL_rm*f_u_rm*C_rm_exv;

[TABLE]

// check for mass balance
capture KDtot = A_kd_vas + A_kd_exv + A_kd_deep;
capture SPtot = A_sp_vas + A_sp_exv + A_sp_deep;
capture RMtot =  A_rm_vas + A_rm_exv; 

capture totaldrug = A_pl + A_gi + A_ht + KDtot + A_li + SPtot + A_lu + RMtot + A_clear;


[capture]
C_pl,Qco
