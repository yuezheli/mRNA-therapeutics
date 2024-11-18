[PROB]

Model taken from Kagan et al., 2013
https://pubmed.ncbi.nlm.nih.gov/23793994/

This file only contains the model for liposomal form of drug; all terms related to the small molecule drug is being removed
  
Remove special dynamics on liver and spleen.

[SET]

delta = 0.1, end = 24

[CMT]

//---- Liposomal compartments ----//
// plasma
A_pl_LIP 
  
// GI
A_gi_vas_LIP
A_gi_exv_LIP
  
// heart
A_ht_vas_LIP
A_ht_exv_LIP
  
// spleen
A_sp_vas_LIP
A_sp_exv_LIP
  
// liver
A_li_vas_LIP
A_li_exv_LIP
  
// kidney
A_kd_vas_LIP
A_kd_exv_LIP
  
// lung
A_lu_vas_LIP
A_lu_exv_LIP
  
// remainder
A_rm_vas_LIP
A_rm_exv_LIP
  
//---- Dummy compartments ----//
// drug clearance; 
A_clear
  
  
  [PARAM]
//---- cardiac output (mus) ----//
Qco = (100-45)/100 * 14.1 * (24/1000)^0.75; // unit: L.h-1, rat HCT value from https://www.mdpi.com/1422-0067/19/9/2824/htm TABLE S1

// fraction of flows compared to cardiac output
q_frac_li = 16.1/100 // the volume that leaves liver
q_frac_kd = 9.1/100
q_frac_sp = 1.125/100
q_frac_gi = 12.87/100
q_frac_ht = 6.6/100

//---- Tissue volume (mus) ----//

wt = 24/1000 // assuming mouse weight = 24g 

// note tissue volume is a percentage of body weight
// assume the volume unit should be L
// the calculation matches rat tissue volume from https://dmd.aspetjournals.org/content/35/10/1766
v_frac_li = 5.49/100
v_frac_kd = 1.67/100
v_frac_sp = 0.35/100
v_frac_gi = 4.22/100
v_frac_ht = 0.5/100
v_frac_lu = 0.73/100
v_frac_blood = 4.9/100 // value from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243854/

//---- fraction of vascular space (mus) ----//
vs_frac_li = 21/100
vs_frac_kd = 24/100
vs_frac_sp = 17/100
vs_frac_gi = 19/100
vs_frac_ht = 26/100
vs_frac_lu = 50/100
vs_frac_mu = 4/100
vs_frac_rm = 4/100

//---- liposomal uptake ----//
// unit: L.h-1; MOUSE PARAMETER
UPgi = 2.04e-4
UPsp = 5.95e-5
UPli = 4.62e-4
UPrm = 1.97e-5
UPlu = 0
UPkd = 0
UPht = 0 

// siRNA-LNP clearance rate
CLplasma = 0.0035
CLtissue = 0.0035 // unit: h-1; temporarily used as the same value

[MAIN]

// flow through organs
// double Qha = Qco * q_frac_ha; 
double Qli = Qco * q_frac_li;
double Qkd = Qco * q_frac_kd;
double Qsp = Qco * q_frac_sp;
double Qgi = Qco * q_frac_gi;
double Qht = Qco * q_frac_ht;

double Qha = Qli - Qsp - Qgi;
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
// list of liposomal concentrations
double C_pl_LIP = A_pl_LIP/ V_pl;

// GI
double C_gi_vas_LIP = A_gi_vas_LIP/ V_gi_vas;
double C_gi_exv_LIP = A_gi_exv_LIP/ V_gi_exv; 

// heart
double C_ht_vas_LIP = A_ht_vas_LIP/ V_ht_vas; 
double C_ht_exv_LIP = A_ht_exv_LIP/ V_ht_exv;

// spleen
double C_sp_vas_LIP = A_sp_vas_LIP/ V_sp_vas;
double C_sp_exv_LIP = A_sp_exv_LIP/ V_sp_exv;

// liver
double C_li_vas_LIP = A_li_vas_LIP/ V_li_vas;
double C_li_exv_LIP = A_li_exv_LIP/ V_li_exv;

// kidney
double C_kd_vas_LIP = A_kd_vas_LIP/ V_kd_vas;
double C_kd_exv_LIP = A_kd_exv_LIP/ V_kd_exv;

// lung
double C_lu_vas_LIP = A_lu_vas_LIP/ V_lu_vas; 
double C_lu_exv_LIP = A_lu_exv_LIP/ V_lu_exv;

// remainder
double C_rm_vas_LIP = A_rm_vas_LIP/ V_rm_vas;
double C_rm_exv_LIP = A_rm_exv_LIP/ V_rm_exv;


//---- Liposomal compartments ----//

// plasma
dxdt_A_pl_LIP = Qco * C_lu_vas_LIP - Qco*C_pl_LIP - CLplasma*C_pl_LIP*V_pl;

// GI
dxdt_A_gi_vas_LIP = Qgi*C_pl_LIP - Qgi * C_gi_vas_LIP - UPgi*C_gi_vas_LIP - CLplasma*C_gi_vas_LIP*V_gi_vas;
dxdt_A_gi_exv_LIP = UPgi*C_gi_vas_LIP - CLtissue*C_gi_exv_LIP*V_gi_exv; 

// heart
dxdt_A_ht_vas_LIP = Qht*C_pl_LIP - Qht * C_ht_vas_LIP - UPht*C_ht_vas_LIP - CLplasma*C_ht_vas_LIP*V_ht_vas;
dxdt_A_ht_exv_LIP = UPht*C_ht_vas_LIP - CLtissue*C_ht_exv_LIP*V_ht_exv; 

// spleen
dxdt_A_sp_vas_LIP = Qsp*C_pl_LIP - Qsp * C_sp_vas_LIP - UPsp*C_sp_vas_LIP - CLplasma*C_sp_vas_LIP*V_sp_vas;
dxdt_A_sp_exv_LIP = UPsp*C_sp_vas_LIP - CLtissue*C_sp_exv_LIP*V_sp_exv; 

// liver
dxdt_A_li_vas_LIP = Qha*C_pl_LIP + Qsp * C_sp_vas_LIP + Qgi * C_gi_vas_LIP - Qli * C_li_vas_LIP  - UPli*C_li_vas_LIP - CLplasma*C_li_vas_LIP*V_li_vas;
dxdt_A_li_exv_LIP = UPli*C_li_vas_LIP - CLtissue*C_li_exv_LIP*V_li_exv;

// kidney
dxdt_A_kd_vas_LIP = Qkd*C_pl_LIP - Qkd * C_kd_vas_LIP - UPkd*C_kd_vas_LIP - CLplasma*C_kd_vas_LIP*V_kd_vas;
dxdt_A_kd_exv_LIP = UPkd*C_kd_vas_LIP - CLtissue*C_kd_exv_LIP*V_kd_exv; 

// lung
dxdt_A_lu_vas_LIP = Qli * C_li_vas_LIP + Qht * C_ht_vas_LIP + Qkd * C_kd_vas_LIP + Qrm * C_rm_vas_LIP - Qco * C_lu_vas_LIP - UPlu*C_lu_vas_LIP - CLplasma*C_lu_vas_LIP*V_lu_vas;
dxdt_A_lu_exv_LIP = UPlu*C_lu_vas_LIP - CLtissue*C_lu_exv_LIP*V_lu_exv; 

// remainder
dxdt_A_rm_vas_LIP = Qrm*C_pl_LIP - Qrm * C_rm_vas_LIP - UPrm*C_rm_vas_LIP - CLplasma*C_rm_vas_LIP*V_rm_vas;
dxdt_A_rm_exv_LIP = UPrm*C_rm_vas_LIP - CLtissue*C_rm_exv_LIP*V_rm_exv; 



// dummy variable 
dxdt_A_clear = CLtissue*C_gi_exv_LIP*V_gi_exv + CLtissue*C_ht_exv_LIP*V_ht_exv + CLtissue*C_sp_exv_LIP*V_sp_exv + 
               CLtissue*C_li_exv_LIP*V_li_exv + CLtissue*C_kd_exv_LIP*V_kd_exv + CLtissue*C_lu_exv_LIP*V_lu_exv + CLtissue*C_rm_exv_LIP*V_rm_exv + 
               CLplasma*C_pl_LIP*V_pl + CLplasma*C_gi_vas_LIP*V_gi_vas + CLplasma*C_ht_vas_LIP*V_ht_vas + CLplasma*C_sp_vas_LIP*V_sp_vas + 
               CLplasma*C_li_vas_LIP*V_li_vas + CLplasma*C_kd_vas_LIP*V_kd_vas + CLplasma*C_lu_vas_LIP*V_lu_vas + CLplasma*C_rm_vas_LIP*V_rm_vas; 

[CAPTURE]
C_pl_LIP, 
C_gi_vas_LIP, C_gi_exv_LIP, 
C_ht_vas_LIP, C_ht_exv_LIP, 
C_sp_vas_LIP, C_sp_exv_LIP, 
C_li_vas_LIP, C_li_exv_LIP, 
C_kd_vas_LIP, C_kd_exv_LIP, 
C_lu_vas_LIP, C_lu_exv_LIP, 
C_rm_vas_LIP, C_rm_exv_LIP
  
[TABLE]

capture Atot = A_pl_LIP + A_gi_vas_LIP + A_gi_exv_LIP + A_ht_vas_LIP + A_ht_exv_LIP + 
  A_sp_vas_LIP + A_sp_exv_LIP + A_li_vas_LIP + A_li_exv_LIP + A_kd_vas_LIP + A_kd_exv_LIP + 
  A_lu_vas_LIP + A_lu_exv_LIP + A_rm_vas_LIP + A_rm_exv_LIP + A_clear; 