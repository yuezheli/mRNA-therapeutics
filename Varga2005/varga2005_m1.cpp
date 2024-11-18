[PROB]

This model file is adapted from Varga et al., 2005
https://www.nature.com/articles/3302495

All plasmid unpacking in cytosol is dropped. Also nuclear import only happens to intact plasmid + vector.

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

Protein

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
k_NPC = 1e3 // nuclear pore association; unit min-1
k_in = 3e-3 // nuclear pore import; unit min-1
k_dissociation = 1e-3 // import protein dissociation within the nucleus; unit min-1

//--- other parameters ---//
k_protein = 1e-2 // protein production; unit min-1

[ODE]

dxdt_Complex_extracellular = -k_bind_uptake*Complex_extracellular; 

// cytoplasmic
dxdt_Complex_internal = k_bind_uptake*Complex_extracellular - k_escape*Complex_internal - k_deg_vesicle*Complex_internal; 

dxdt_Complex_cytoplasmic = k_escape*Complex_internal - k_bind_vector*Complex_cytoplasmic;

dxdt_ComplexBound_cytoplasmic = k_bind_vector*Complex_cytoplasmic - k_NPC*ComplexBound_cytoplasmic; 

// NPC-related status

dxdt_ComplexBound_NPC = k_NPC*ComplexBound_cytoplasmic - k_in*ComplexBound_NPC; 

// nuclear

dxdt_ComplexBound_nuclear = k_in*ComplexBound_NPC - k_dissociation*ComplexBound_nuclear;

dxdt_Complex_nuclear = k_dissociation*ComplexBound_nuclear - k_unpack*Complex_nuclear;

dxdt_Plasmid_nuclear = k_unpack*Complex_nuclear;

// protein synthesis

dxdt_Protein = k_protein*Plasmid_nuclear; 