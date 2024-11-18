[PROB]

This model file creates the model published in Varga et al., 2005
https://www.nature.com/articles/3302495

Note this implementation does not track vector dynamics. In addition, the tracking on degredaded plasmid is adjusted to reflect the biology. 

Furthermore, extracellular complex (i.e. vector + plasmid) is modeled as a dummy variable to track complex internalization.

[SET]

delta = 10, end = 60*24*3 // time in minutes

[CMT]
Complex_extracellular
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

Protein
X_plasmid // plasmid that has been degraded

[PARAM]
//--- parameters that are vector specific ---//
// (default, Ad5)
k_bind_uptake = 6e-3
k_deg_vesicle = 2e-2 // lysosomal degredaion
k_escape = 1.6e-2 // endosomal escape; unit min-1
k_bind_vector = 1e-1
k_unpack = 1e-2 // vector unpacking; unit min-1

//--- parameters that holds throughout ---//
// inherited parameters
k_bind_plasmid = 2e-3 // formation of nuclear import protein bound vector; unit min-1
k_degredation = 5e-3 // plasmid degredation; unit min-1
k_NPC = 1e3 // nuclear pore association; unit min-1
k_in = 3e-3 // nuclear pore import; unit min-1
k_dissociation = 1e-3 // import protein dissociation within the nucleus; unit min-1

//--- other parameters ---//
k_protein = 1e-2 // protein production; unit min-1

[ODE]

dxdt_Complex_extracellular = -k_bind_uptake*Complex_extracellular; 

// cytoplasmic
dxdt_Complex_internal = k_bind_uptake*Complex_extracellular - k_escape*Complex_internal - k_deg_vesicle*Complex_internal; 

dxdt_Complex_cytoplasmic = k_escape*Complex_internal - k_bind_vector*Complex_cytoplasmic - k_unpack*Complex_cytoplasmic;

dxdt_Plasmid_cytoplasmic = k_unpack*Complex_cytoplasmic - k_bind_plasmid*Plasmid_cytoplasmic - k_degredation*Plasmid_cytoplasmic;

dxdt_PlasmidBound_cytoplasmic = k_bind_plasmid*Plasmid_cytoplasmic - k_NPC*PlasmidBound_cytoplasmic;

dxdt_ComplexBound_cytoplasmic = k_bind_vector*Complex_cytoplasmic - k_NPC*ComplexBound_cytoplasmic; 

// NPC-related status

dxdt_ComplexBound_NPC = k_NPC*ComplexBound_cytoplasmic - k_in*ComplexBound_NPC; 

dxdt_PlasmidBound_NPC = k_NPC*PlasmidBound_cytoplasmic - k_in*PlasmidBound_NPC;

// nuclear

dxdt_ComplexBound_nuclear = k_in*ComplexBound_NPC - k_dissociation*ComplexBound_nuclear;

dxdt_PlasmidBound_nuclear = k_in*PlasmidBound_NPC - k_dissociation*PlasmidBound_nuclear;

dxdt_Complex_nuclear = k_dissociation*ComplexBound_nuclear - k_unpack*Complex_nuclear;

dxdt_Plasmid_nuclear = k_unpack*Complex_nuclear + k_dissociation*PlasmidBound_nuclear;

// protein synthesis

dxdt_Protein = k_protein*Plasmid_nuclear; 

// degredaded plasmid

dxdt_X_plasmid = k_deg_vesicle*Complex_internal + k_degredation*Plasmid_cytoplasmic; 

[TABLE]

capture total_plasmid_nuclear =  ComplexBound_NPC + ComplexBound_nuclear + Complex_nuclear + Plasmid_nuclear + PlasmidBound_NPC + PlasmidBound_nuclear; 

capture total_plasmid_cytoplasmic = Complex_internal + Complex_cytoplasmic + ComplexBound_cytoplasmic + Plasmid_cytoplasmic + PlasmidBound_cytoplasmic;

capture total_plasmid = total_plasmid_nuclear + total_plasmid_cytoplasmic;

capture MassBalancePlasmid = Complex_extracellular + total_plasmid + X_plasmid; 