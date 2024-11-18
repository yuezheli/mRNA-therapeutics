[PROB]

Simplied model from Apgar et al., 2018. 
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6391595/

what is included: 
1. LNP dynamics
2. mRNA dynamics
3. protein expression

[SET]

delta = 10, end = 60*60*24*3

[CMT]
// mass of species; unit in nmol
LNP       // LNP in plasma
LNPp      // LNP in peripheral tissues
LNPa      // LNP attached to hepatocyte
LNPe      // endocytosed LNP
mRNAc     // cytoplasmic mRNA
protein   // cytoplasmic protein

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


Vc = 0.0078        // plasma volume; L-1

//----- parameters that is not clear -----//
ktbg = 1E-5        // endogenous translation of UGT

Vhep = 5.85E-3     // calculation based on google doc description; L-1

//----- parameters related to dosing -----//
animal_weight = 0.4 // assumes Gunn rat weight ~ 400g; https://pubmed.ncbi.nlm.nih.gov/16487915/
dosing = 0.3        // dosing amount; mg/kg

moleweight_LNP = 1 // estimation based on https://pubs.acs.org/doi/10.1021/mp500367k; mg.nmol-1

[MAIN]
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
dxdt_protein =  ktbg + kt * mRNAc - dUGTc * protein;


[TABLE]

capture PlasmaDrug = LNP;
capture mRNA = mRNAc + LNPa + LNPe;

