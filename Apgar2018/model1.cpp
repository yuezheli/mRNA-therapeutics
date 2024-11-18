[PROB]

Implementation of the Apgar et al., 2018. 
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6391595/

[SET]

delta = 10

[CMT]
// mass of species; unit in nmol
LNP       // LNP in plasma
LNPp      // LNP in peripheral tissues
LNPa      // LNP attached to hepatocyte
LNPe      // endocytosed LNP
mRNAc     // cytoplasmic mRNA
UGTc      // cytoplasmic UGT
Bil       // bilirubin
Bil_UGTc  // Bilirubin-UGTc complex
MGT       // monglucouronide bilirubin
DGT       // diglucouronide bilirubin

// dummy variables; A.U.
sBil      // the high production of bilirubin
running   // the input for junior rats

[PARAM]

//----- parameters that are the same between rats and human -----//

// LNP dynamics
kw = 2.41E-5    // first order elimination of LNP; s-1
k12 = 4.79E-5   // Vc -> Vp distribution rate; s-1
k21 = 2.65E-7   // Vp -> Vc distribution rate; s-1
ka = 1.17E-5    // LNP attachment to hepatocyte; s-1
ke = 7.7E-5     // LNP endocytosis; s-1
de = 9.32E-5    // LNP degradation in endosome; s-1

// mRNA related parameters
kl = 1.93E-5      // endosomal escape rate, s-1
dmRNA = 1.07E-5   // mRNA degradation; s-1
kt = 17.73        // translation rate; s-1

// protein
dUGTc = 6.76E-6  // cytoplasmic protein degradation rate; s-1
kon = 1E-5       // protein binding on rate; nmol-1
koff = 0.2589    // protein unbinding off rate; s-1
kcat = 0.0011    // enzyme kcat rate for glucuronidation; s-1


//----- parameters that are different between rats and human -----//
kclearBil = 3.5E-6 // bilirubin clearance rate
kclearMGT = 3.5E-5 // elimination of monoglucuronide bilirubin; s-1
kclearDGT = 3.5E-5 // elimination of diglucuronide bilirubin; s-1

Vc = 0.0078        // plasma volume; L-1

moleweight_bili = 5.85E-4 // bilirubin molecular weight; mg.nmol-1
moleweight_UGT3 = 0.175 // UGT molecular weight, 175kDa; mg.nmol-1

//----- parameters that is not clear -----//
ksyn = 1E-5        // basal synthesis of bilirubin; this parameter needs to be adjusted; s-1
ksynhigh = 1E-5    // high synthesis of bilirubin
ktbg = 1E-5        // endogenous translation of UGT
kelSbil = 1E-5     // production of bilirubin

Vhep = 5.85E-3     // calculation based on google doc description; L-1

// value for sBil to control whether the rats is an adult
init_sBil = 1;
init_running = 0;

//----- parameters related to dosing -----//
animal_weight = 0.4 // assumes Gunn rat weight ~ 400g; https://pubmed.ncbi.nlm.nih.gov/16487915/
dosing = 0.3        // dosing amount; mg/kg

moleweight_LNP = 1 // estimation based on https://pubs.acs.org/doi/10.1021/mp500367k; mg.nmol-1

[MAIN]

// set initial value
sBil_0 = init_sBil; 
running_0 = init_running;

// dosing regiment
double dose = dosing * animal_weight / moleweight_LNP; // convert dose to nmol; mass

LNP_0 = dose; 

[ODE]

// mass balance of LNP
dxdt_LNP = -kw * LNP - k12 * LNP + k21 * LNPp - ka * LNP; 
dxdt_LNPp = k12 * LNP - k21 * LNPp - kw * LNPp; 
dxdt_LNPa = ka * LNP - ke * LNPa; 
dxdt_LNPe = ke * LNPa - de * LNPe - kl * LNPe; 

// mass balance of mRNA (stoichiometry may be problematic)
dxdt_mRNAc = kl * LNPe - dmRNA * mRNAc - kt * mRNAc + kt * mRNAc; 
  
// mass balance of protein
// dxdt_UGTc =  ktbg + kt * mRNAc - dUGTc * UGTc - kon * Bil * UGTc + koff * Bil_UGTc + kcat * Bil_UGTc;
// dxdt_Bil_UGTc =  kon * Bil * UGTc - dUGTc * Bil_UGTc - kcat * Bil_UGTc; 
dxdt_UGTc =  ktbg + kt * mRNAc - dUGTc * UGTc - kon/Vhep * Bil * UGTc + koff * Bil_UGTc + kcat * Bil_UGTc;
dxdt_Bil_UGTc =  kon/Vhep * Bil * UGTc - dUGTc * Bil_UGTc - kcat * Bil_UGTc;
dxdt_MGT = kcat * Bil_UGTc - kclearMGT * MGT - kcat * UGTc * MGT; 
dxdt_DGT = kcat * UGTc * MGT - kclearDGT * DGT; 
dxdt_Bil= ksyn + ksynhigh * sBil - kclearBil * Bil - kon/Vhep * Bil * UGTc + dUGTc * Bil_UGTc;
dxdt_sBil = - ksynhigh * sBil  + ksynhigh * sBil - kelSbil * sBil * running  ; // high synthesis of Bilirubin
dxdt_running = - kelSbil * sBil * running;
// issues: check kon, koff, kcat related terms and their units

[TABLE]

capture TotalBilirubin = (Bil + MGT + DGT) ; // unit: nmol/L
capture BilirubinBlood = Bil * moleweight_bili / (10 * Vc); // unit: mg/dL
capture BiliProd = ksyn + ksynhigh * sBil;  // unit: nmol
capture PlasmaDrug = LNP;
capture mRNA = mRNAc + LNPa + LNPe;
capture Enzyme = UGTc + Bil_UGTc; 

