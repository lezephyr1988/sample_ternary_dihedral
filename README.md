# sample_ternary_dihedral
Sampling dihedral for non-crash E2/POI-contacting conformation in E3 ternary complex

    DESCRIPTION:
        randomly sampling flexible region for E3/PROTAC/POI complex for identify finds non-crach E3/POI-contacting conformation.
    
    PARAMETERS:
        motion_part, the object containing RBX1, E2 and UB
        fixed_part, the object containing PROTAC, VHL, POI, Cullin, NEDD8
        flexible_residue_string, list of flexible residues in motion part (use ; as seperator, for example "resi 2;resi 3;resi 5")
        E2_residue, the catalytic Cysteine of E2 (motion_part)
        POI_residue, the ubiquination site (lysine) of POI (fixed_part)
        crash_mask, residues in motion_part that are ignored in calculating crash
        crash_radius, the minimal distance between motion_part and fixed_part
        max_crash_number, the maximal acceptable number of crash. It could be zero for restrictive calculation, and it also could be relax for further manually refinement
        E2_POI_distance, the maximal distance between side chain of E2_residue and POI_residue
        number, number for randomly sampling (about 10000 to 100000 for real project)
        initialize_degree, 0-360. 0 for starting search from the starting structure, 360 for fully disturb the starting structure and starting from random conformation
        file_prefix,    string, prefix for writing sampled conformations meet requirment.
        
   
    RETURNS:
        list of distance for conformations meet requirment
        
    Example:
        After opening the sample_ternary_dihedral.pse, load the sample_ternary_dihedral.pml file in pymol, and then run:
        -------------
        sample_ternary_dihedral("motion_part","fixed_part","motion_part & chain E & resi 36;motion_part & chain E & resi 37;motion_part & chain E & resi 38;motion_part & chain E & resi 39;motion_part & chain E & resi 40","motion_part & chain F & resi 85", "fixed_part & chain B & resi 364","motion & chain E & resi 23-40", 3.0,200, 15.0, 100, 0, "sample")   
        -------------
        If any conformation meet the requirment is captured, it will be output to sample_N.pdb
        
        Note1: For quick demo, the E2_POI_distance in the example is 15.0, the running number is only 100, the initialize_degree is set to 0. 
        Note2: All proteins in the motion_part should be connect by a dummy bond (this could make all protein in the set move together). It the example, it could be set by "bond last motion_part& chain E, first motion_part& chain F", and the bond could be removed by "unbond last motion_part& chain E, first motion_part& chain F"
          
    AUTHOR:
        Yaxia Yuan, 2021.  
