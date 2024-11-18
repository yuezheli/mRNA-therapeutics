[PROB]

This model file creates the model published in Varga et al., 2001
https://pubmed.ncbi.nlm.nih.gov/11708880/

[GLOBAL]


[SET]

delta = 10, end = 60*24*3 // time in minutes

[CMT]
// complex stands for rAAV (or lipid packages) + plasmid
Complex_internal
Complex_cytoplasmic
ComplexBound_cytoplasmic
ComplexBound_NPC
ComplexBound_nuclear
Complex_nuclear

// plasmid is the trans gene
Plasmid_nuclear
Plasmid_cytoplasmic
PlasmidBound_cytoplasmic
PlasmidBound_NPC
PlasmidBound_nuclear

// vector is the rAAN protein capsid
Vector_nuclear
Vector_cytoplasmic
VectorBound_cytoplasmic
VectorBound_NPC
VectorBound_nuclear

// other
X_Plasmid_cytoplasmic 
Protein

[PARAM]

// k_internalization = 1 // placeholder for the debugging

k_escape = 1e-2 // endosomal escape; unit min-1
k_unpack = 1e9 // vector unpacking; unit min-1
k_bind = 2e-3 // formation of nuclear import protein bound vector and/ or plasmid; unit min-1
k_NPC = 1e3 // nuclear pore association; unit min-1
k_in = 3e-3 // nuclear pore import; unit min-1
k_dissociation = 1e-3 // import protein dissociation within the nucleus; unit min-1
k_degredation = 5e-3 // plasmid degredation; unit min-1
k_protein = 1e-2 // protein production; unit min-1

ComplexTotal = 9e14 // total number of plasmid

[ODE]

// cytoplasmic

double k_internalization = ComplexTotal * exp(-SOLVERTIME);

dxdt_Complex_internal = k_internalization - k_escape*Complex_internal; 

dxdt_Complex_cytoplasmic = k_escape*Complex_internal - k_bind*Complex_cytoplasmic - k_unpack*Complex_cytoplasmic;

dxdt_Vector_cytoplasmic = k_unpack*Complex_cytoplasmic - k_bind*Vector_cytoplasmic; // note no VECTOR DEGREDATION

dxdt_VectorBound_cytoplasmic = k_bind*Vector_cytoplasmic - k_NPC*VectorBound_cytoplasmic; // added equation

dxdt_Plasmid_cytoplasmic = k_unpack*Complex_cytoplasmic - k_bind*Plasmid_cytoplasmic - k_degredation*Plasmid_cytoplasmic;

dxdt_PlasmidBound_cytoplasmic = k_bind*Plasmid_cytoplasmic - k_NPC*PlasmidBound_cytoplasmic;

dxdt_ComplexBound_cytoplasmic = k_bind*Complex_cytoplasmic - k_NPC*ComplexBound_cytoplasmic; 

dxdt_X_Plasmid_cytoplasmic = k_degredation*Plasmid_cytoplasmic - k_bind*X_Plasmid_cytoplasmic; // note: the second term is really weird

// NPC-related status

dxdt_ComplexBound_NPC = k_NPC*ComplexBound_cytoplasmic - k_in*ComplexBound_NPC; 

dxdt_PlasmidBound_NPC = k_NPC*PlasmidBound_cytoplasmic - k_in*PlasmidBound_NPC;

dxdt_VectorBound_NPC = k_NPC*VectorBound_cytoplasmic - k_in*VectorBound_NPC; // added equation

// nuclear

dxdt_ComplexBound_nuclear = k_in*ComplexBound_NPC - k_dissociation*ComplexBound_nuclear;

dxdt_PlasmidBound_nuclear = k_in*PlasmidBound_NPC - k_dissociation*PlasmidBound_nuclear;

dxdt_Complex_nuclear = k_dissociation*ComplexBound_nuclear - k_unpack*Complex_nuclear;

dxdt_VectorBound_nuclear = k_in*VectorBound_NPC - k_dissociation*VectorBound_nuclear; 

dxdt_Vector_nuclear = k_unpack*Complex_nuclear + k_dissociation*VectorBound_nuclear; 

dxdt_Plasmid_nuclear = k_unpack*Complex_nuclear + k_dissociation*PlasmidBound_nuclear;

// protein synthesis

dxdt_Protein = k_protein*Plasmid_nuclear; 

[CAPTURE]

k_internalization

[TABLE]

capture total_plasmid_nuclear = ComplexBound_NPC + ComplexBound_nuclear + Complex_nuclear + Plasmid_nuclear + PlasmidBound_NPC + PlasmidBound_nuclear; 

capture total_plasmid_cytoplasmic = Complex_cytoplasmic + ComplexBound_cytoplasmic + Plasmid_cytoplasmic + PlasmidBound_cytoplasmic;

capture total_plasmid = total_plasmid_nuclear + total_plasmid_cytoplasmic;