# -*- coding: utf-8 -*-
#
# sample_ternary_dihedral.py -- sample dihedral for non-crash E2/POI-contacting conformation
#
import random
def sample_ternary_dihedral( motion_part, fixed_part, flexible_residue_string, E2_residue, POI_residue,crash_mask, crash_radius, max_crash_number, E2_POI_distance, number, initialize_degree, file_prefix):
    """
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
    
    """
    distance_list=[]
    angle_list=[]
    current_angle_list=[]
    meet_number=0
    flexible_residue=flexible_residue_string.split(";")
    angle_scale=0.05 #initally: 18 degree rotation
    #add backbone hydrogen/set name
    for residue in flexible_residue:
        cmd.h_add(residue+" & name N")
        cmd.alter("hydrogen within 1.5 of "+residue+" &name N",'name="H"')

        angle=cmd.get_dihedral(residue+" & name H",residue+" & name N", residue+" & name CA",residue+" & name C")
        angle_list.append(angle)
        current_angle_list.append((1-random.random())*initialize_degree+angle)
        
        angle=cmd.get_dihedral(residue+" & name N",residue+" & name CA", residue+" & name C",residue+" & name O")
        angle_list.append(angle)
        current_angle_list.append((1-random.random())*initialize_degree+angle)
        
    for i in range(number):
        #randomly change phi/psi for all residues in the flexible residue list
        current_num=0
        for residue in flexible_residue:
            current_angle_list[current_num]+=(1.0-random.random())*360.0*angle_scale
            current_angle_list[current_num+1]+=(1.0-random.random())*360.0*angle_scale
            cmd.set_dihedral(residue+" & name H",residue+" & name N", residue+" & name CA",residue+" & name C",current_angle_list[current_num])
            cmd.set_dihedral(residue+" & name N",residue+" & name CA", residue+" & name C",residue+" & name O",current_angle_list[current_num+1])
            current_num+=2
        # check whether it meet the E2/POI distance
        distance=cmd.get_distance(E2_residue+" & name SG",POI_residue+" & name NZ")
        if distance <= E2_POI_distance:
            #check whether it has crash
            crash_select="("+motion_part+") and not ("+crash_mask+")"
            fixed_select="("+fixed_part+") and not ("+POI_residue+")"
            crash_number=cmd.select("("+crash_select+" within "+str(crash_radius)+" of "+fixed_part+")")
            if crash_number<=max_crash_number: #ok
                print("Conformation",i+1,"Distance",distance,"Success")
                distance_list.append(distance)
                meet_number+=1
                #save result
                cmd.save(file_prefix+"_"+str(meet_number)+".pdb",motion_part+" or "+fixed_part)
                angle_scale=0.1 #36degree for starting from new conformation
            else:
                print("Conformation",i+1,"Distance",distance,"Crashed",crash_number)
                angle_scale=0.005 #1.8degree for slow change
               
        else:
            print("Conformation",i+1,"Distance",distance)
            angle_scale=0.05*(1-E2_POI_distance/distance)+0.001       #change the rotation level, if it near to the desired distance, the rotation will become slow
            
    #reset dihedral
    reset_number=0
    for residue in flexible_residue:
        cmd.set_dihedral(residue+" & name H",residue+" & name N", residue+" & name CA",residue+" & name C",angle_list[reset_number])
        cmd.set_dihedral(residue+" & name N",residue+" & name CA", residue+" & name C",residue+" & name O",angle_list[reset_number+1])
        reset_number+=2
    print("Total",meet_number,"conformations outputed")
    for i in range(len(distance_list)):
        print(file_prefix+"_"+str(i+1)+".pdb","Distance",distance_list[i])

#example           
#sample_ternary_dihedral("motion_part","fixed_part","motion_part & chain E & resi 36;motion_part & chain E & resi 37;motion_part & chain E & resi 38;motion_part & chain E & resi 39;motion_part & chain E & resi 40","motion_part & chain F & resi 85", "fixed_part & chain B & resi 364","motion & chain E & resi 23-40", 3.0,200, 12.0, 100,0, "sample")
cmd.extend("sample_ternary_dihedral", sample_ternary_dihedral)