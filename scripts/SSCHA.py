#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys,os
from numpy import *
import numpy as np
import time
import os
import contextlib

import ase

from ase.calculators.espresso import Espresso
from ase.visualize import view


import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.Cluster

from utils_sscha_mlip.MTP_environment import MTP_Environment
from utils_sscha_mlip.Cluster_Management import Cluster_Management
import utils


#IMPORTANT PARAMETERS FOR THE WORKFLOW


GENERATE_ENS = True                                #if True, it will generate new displaced configurations. if False it will expect to find configurations already generated
T            = 0                                   #TEMPERATURE OF THE SIMULATION
N_RANDOM     = 4                                   #Number of population for each configuration
POPULATION   = 1                                   #Initial population of the simulation
POP_MAX      = 10                                  #Final population of the simulation
PREFIX       = 'PdCuH_'                            #PREFIX OF THE DYNAMICAL MATRICES
MATDYN       = PREFIX + str(POPULATION-1) + '.dyn' #NAME OF THE MATRICES
SUPERCELL    = (2,2,2)                             #size of the supercell 
NQIRR        = 6                                   # Number of irreducible matrices
ENS_FOLDER   = 'TMP_gen'                           #FOLDER FOR THE FORCES,PRESSURES,AND CONFIGURATIONS
KL_VAR       = 0.4                                 #Kong-Liu parameter (see SSCHA)
MAX_KA       = 500                                 #Maximum number of steps for the SSCHA
MIMIM_STRUCTURE = False                            #Minimizes the atomic position if True

RELAX_FLAG          = False                        #Relaxes the cell if True
static_bulk_modulus = 10                           #Guess for the bulk modulus in GPa
target_press        = 40                           #Target pressure in GPa
fix_cell_shape      = False                        #perform an isotropic variation of the cell if True



folder     = os.getcwd() + "/SSCHA"                #Folder where to execute the dft calculations


header = utils.SCF_FILE                           #Header of the Quantum Espresso scf inputs

CALCULATOR = "ESPRESSO"                           # either "ESPRESSO" or "VASP" This flag chooses the program for the self consistent calculations
ADDRESS    = "Your cluster address "              #Your cluster address if you want to run your DFT calculations on a cluster from a local machine
marker     = "PdCuH"    #str(os.urandom(2).hex()) #Header for the scf jobs   

Submitter  = Cluster_Management(PREFIX,POPULATION,ENS_FOLDER,ADDRESS,folder, marker,N_RANDOM)   #Starting an environment for the Cluster management. Check the library for further details. By default this routine will run a local calculation. 
Submitter.base_script_VASP = utils.base_script_VASP           #Assings some slurm scripts for the scf calculations. Check utils.py file
Submitter.base_script_QE = utils.base_script_QE               #Assings some slurm scripts for the scf calculations. Check utils.py file
utils.create_scripts(folder,CALCULATOR)                       #generates folders 

TRAINED_FLAG = False                                                                   #Place true only if the SSCHA run interrupted after the MTP potential has already been trained
GAMMA        = 200                                                                     #Value of the Gamma selet for the active learning protocol
MLIP_PATH    = '/projects/academic/ezurek/francesco/mlip-2-master/bin/mlp'             #path to the MLIP2 execitable
MTP_handler  = MTP_Environment(POPULATION,ENS_FOLDER,GAMMA,ADDRESS,folder,MLIP_PATH)   #Prepares an handler to manage the MLIP interface
MTP_handler.TRAIN_FLAG = "test_tail"
MTP_handler.nprocs = 1                                                                 #Number of processors for the parallel division of the MLIP calculation 







#IMPORTANT VARIABLES THAT NEED TO BE STORED THROUGHOUT THE EXECUTION OF THE SCRIPT
# energies (array of energy) contains all the dft/Mlp energies that need to be combined


dyn = CC.Phonons.Phonons(MATDYN, NQIRR) #LOAD THE DYN MAT IN MEMORY
if POPULATION == 1:
    dyn.Symmetrize()                       #IN FIRST STEP ONLY: APPLIES SUM RULE
    dyn.ForcePositiveDefinite()            #IN FIRST STEP ONLY: FLIPS NEGATIVE FREQUENCIES
    dyn.save_qe( PREFIX + str(POPULATION-1) + ".dyn")                 #iN FIRST STEP ONLY: SAVES NEW DYN 
ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL) #LOADS THE DYN IN THE SSCHA PROGRAM (class ens)
minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
BFGS = sscha.Optimizer.UC_OPTIMIZER(minimizer.dyn.structure.unit_cell)
static_bulk_modulus /= sscha.SchaMinimizer.__evA3_to_GPa__
BFGS.alpha = 1 / (3 * static_bulk_modulus * minimizer.dyn.structure.get_volume())



while POPULATION < POP_MAX:

    MATDYN = PREFIX + str(POPULATION-1) + '.dyn'
    Submitter.POPULATION = POPULATION
    MTP_handler.POPULATION = POPULATION
    
    #################################################################################
    dyn = CC.Phonons.Phonons(MATDYN, NQIRR) #LOAD THE DYN MAT IN MEMORY
    print(dyn.GetSupercell())               #CHECK THE SIZE OF SUPERCELL
    if POPULATION == 1:
        dyn.Symmetrize()                       #IN FIRST STEP ONLY: APPLIES SUM RULE
        dyn.ForcePositiveDefinite()            #IN FIRST STEP ONLY: FLIPS NEGATIVE FREQUENCIES
        dyn.save_qe( PREFIX + str(POPULATION-1) + ".dyn")                 #iN FIRST STEP ONLY: SAVES NEW DYN 
    
    w_s, pols = dyn.DiagonalizeSupercell()  #GETS EIGEN- VALIES AND VECTORS
    print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in  w_s]))
    ##################################################################################
    
    
    ##################################################################################
    ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL) #LOADS THE DYN IN THE SSCHA PROGRAM (class ens)
    if GENERATE_ENS == True:
        ens.generate(N_RANDOM)       #GENERATES THE CONFIGURATIONS BASED ON DYN
        #view(ens.structures[0].get_ase_atoms()) #VIEWS FIRST CONFIGURATION
        # Save the ensemble
        ens.save( ENS_FOLDER + str(POPULATION), POPULATION)  #SAVES THE CONFIGURATIONS ON FILE
    else:
        ens.load(data_dir=ENS_FOLDER + str(POPULATION), population=POPULATION,N=N_RANDOM,load_noncomputed_ensemble=True)
    ##################################################################################



    ##################################################################################
    #PREPARES CFG EMPTY FILES
    if MTP_handler.GAMMA > 0: MTP_handler.Generate_CFG()
    ##########################################################################################################################
    
    ##################################################################################
    #LOADS AND SUBMIT JOBS ON THE CLUSTER
    
    #CHECKS MINDIST JUST TO BE SURE NOTHING IS GOING WRONG
    if MTP_handler.GAMMA > 0: MTP_handler.Calc_GRADE()
    print("DONE")
   #####################################################################################
  
   
    
    print("creating the gamma table") 
    ###################################################################################
    if MTP_handler.GAMMA > 0: MTP_handler.Fill_GAMMA_Table()
    else: MTP_handler.Fill_GAMMA_Table_0(N_RANDOM)
    #########################################################################################
   
    
    
    ######################################################################################
    if not os.path.exists("scfin"):
        os.mkdir("scfin")
    #########################################################################################
    
    
    #########################################################################################
    #THIS CELL GENERATES THE INPUTS FOR THE QUANTUM ESPRESSO CALCULATIONS
    
    typical_espresso_header = header.format(ens.structures[0].N_atoms) 
    
    # Now we need to read the scf files
    if CALCULATOR == "ESPRESSO": Submitter.Generate_SCF(typical_espresso_header)
    elif CALCULATOR == "VASP": saved_ordering = Submitter.Generate_POSCAR()
    ##############################################################################################

    
    ##############################################################################################
    #LOADS AND SUBMIT JOBS ON THE CLUSTER
    if GENERATE_ENS == True:
        if CALCULATOR == "ESPRESSO": Submitter.Send_to_Folders(MTP_handler.conf_table_gamma)
        if CALCULATOR == "VASP": Submitter.Send_to_Folders_VASP(MTP_handler.conf_table_gamma)
        print("DONE")
    #################################################################################################
    
    print("going into queue and resubmit")
    #################################################################################################
    if CALCULATOR == "ESPRESSO": Submitter.Queue_and_resubmit(MTP_handler.conf_table_gamma)
    if CALCULATOR == "VASP": Submitter.Queue_and_resubmit_VASP(MTP_handler.conf_table_gamma)
    #######################################################################################################################################################
    
    
    #########################################################################################################################################################
    #LOADS THE forces, pressures AND energies IN THE FILES AND MEMORIES FROM THE QE INPUTS
    
    #As written, we must convert the total energy of the supercell in Ry, the forces in Ry/Bohr, and the stress in Ry/Bohr^3. Luckily quantum espresso already gives these quantities in the correct units, but be careful when using different calculators. This problem does not arise when using automatic calculators, as the SSCHA and ASE will cooperate to convert the units to the correct one. Now we will parse the Quantum ESPRESSO output looking for the energy, the forces, and the stress tensor.
    if CALCULATOR == "ESPRESSO": energies = Submitter.read_SCF()
    if CALCULATOR == "VASP": energies = Submitter.read_OUTCAR(saved_ordering)  
    print("DONE")
    ########################################################################################################################################################
    
    
    
    #########################################################################################################################################################
    #ROUTINE TO GO FROM SSCHA TMP FILE TO MLIP CFG FILE
    
    if TRAINED_FLAG == False and MTP_handler.GAMMA > 0:   
        MTP_handler.Compile_Trainingset()
        MTP_handler.Train_submit()

        print("DONE")
    
    if MTP_handler.GAMMA > 0: TOTAL_NUMBER_OF_MLIP_CONF = MTP_handler.read_CFG(energies)
    else: TOTAL_NUMBER_OF_MLIP_CONF = 0
    ######################################################################################################################################################################
    
        
   












   #######################################################################################################################################################################
    #HERE WE RUN THE MINIMIZATION
    
    
    #ens.update_weights(dyn, T=0) # CHECK THIS COMMAND ON DOCUMENTATION
    ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
    ens.load(data_dir=ENS_FOLDER + str(POPULATION), population=POPULATION,N=N_RANDOM)   #THIS IS TO USE IF YOU ALREADY HAVE THE CONFIGURATIONS AND DON'T WANT TO GENERATE NEW ONES
    
    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ens) #LOADS THE ROUTINE FOR THE MINIMIZATION
    
    minimizer.minim_struct = MIMIM_STRUCTURE
    
    #minimizer.neglect_symmetries = False
    
    minimizer.min_step_dyn = 0.01 #Values around 1 are good
    minimizer.min_step_struc = 0.1
    #minimizer.precond_dyn = False
    #minimizer.root_representation = "root2"
    
    minimizer.kong_liu_ratio = KL_VAR# Usually 0.5 is a good value
    
    minimizer.meaningful_factor = 0.0000001
    minimizer.max_ka = MAX_KA
    #minimizer.minimization_algorithm = "cgrf"
    
    # FUNCTION TO SAVE DATA###############################
    all_frequencies = []
    def add_current_frequencies(minimizer):
        w, p = minimizer.dyn.DiagonalizeSupercell()
        all_frequencies.append(w)
        # In this way the file will be updated at each step
        np.savetxt("all_frequencies.dat", all_frequencies)
    #######################################################
    
    
    IOdata = sscha.Utilities.IOInfo()
    IOdata.SetupSaving("minim")
    
    #original_stdout = sys.stdout #CHANGE STANDARD OUTPUT
    
    #with open('minimi.out', 'w') as f:
    #sys.stdout = f  #BECAUSE I LIKE TO PRINT ON FILE
    minimizer.init() #INITIALIZES MINIMIZATION
    minimizer.run(custom_function_post = IOdata.CFP_SaveAll) #RUN MINIMIZATION
    minimizer.finalize() #PRINT FINAL RESULT
    #sys.stdout = original_stdout # RESETS TO THE STANDARD OUTPUT
    
    #minimizer.plot_results(save_filename="dynamic.dat")  #MAKES SOME PLOTS AND SAVES DATA
    
    final_results_file = open("Final_Results_File", "w")
    with contextlib.redirect_stdout(final_results_file):
        minimizer.finalize()
    final_results_file.close()
 
    
    #HIGHLIGHTS IMPORTANT RESULTS
    print("The total free energy per unit cell is:", minimizer.get_free_energy(), " Ry")
    print("The total stress tensor is [Ry/bohr^3]:")
    print(minimizer.get_stress_tensor()[0])
    print("And the stochastic error on the stress tensor is:")
    print(minimizer.get_stress_tensor()[1])
    print("The stocastic error of the free energy instead, was:", minimizer.get_free_energy(return_error = True)[1], " Ry")
    
    #MAKES A FILE FOR THE FREQUENCIES READABLE BY XMGRACE
    res = os.system('cat "minim.freqs" | nl > "freq.dat" ')
    
    ########################################################################################################################################################################
    
    #view(minimizer.dyn.structure.get_ase_atoms()) #SHOWS FINAL STRUCTURE
    
    minimizer.dyn.save_qe(PREFIX + str(POPULATION) + ".dyn") #SAVES NEW MATRICES ON FILE
    
    #PRINTS FREQUENCY VARIATION
    #w_old, p_old = ens.dyn_0.DiagonalizeSupercell() 
    #w_new, p_new = minimizer.dyn.DiagonalizeSupercell()
    #print(" Old frequencies |  New frequencies")
    #print("\n".join(["{:16.4f} | {:16.4f}  cm-1".format(w_old[i] * CC.Units.RY_TO_CM, w_new[i] * CC.Units.RY_TO_CM) for i in range(len(w_old))]))
    
    ##########################################################################################################################################################################
    
    ##########################################################################################################################################################################
#    if HESSIAN == True:
#        dyn = CC.Phonons.Phonons(PREFIX + str(POPULATION) + ".dyn", NQIRR)
#        ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
#        ens.load(data_dir=ENS_FOLDER + str(POPULATION), population=POPULATION,N=N_RANDOM)   #THIS IS TO USE IF YOU ALREADY HAVE THE CONFIGURATIONS AND DON'T WANT TO GENERATE NEW ONES
#        minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
#        relax = sscha.Relax.SSCHA(minimizer,N_configs = N_RANDOM)
#    
#        print("I am doing the hessian")
#        free_energy_hessian = relax.minim.ensemble.get_free_energy_hessian(include_v4 = False)
#        # Now we can save the free energy hessian like it was a dynamical matrix
#        free_energy_hessian.save_qe("hessian_3d_" + str(N_RANDOM) + ".dyn")
#    
#    # We can print its eigenvalues to see if there are imaginary phonons
#        w, p = free_energy_hessian.DiagonalizeSupercell()
#        print(w * CC.Units.RY_TO_CM)
    ############################################################################################################################################################################
    
    
    minim_file = open("minim.dat","r")
    lines = minim_file.readlines()
    grad_for = float(lines[-1].split()[2])
    grad_pos = float(lines[-1].split()[4])
    minim_file.close()
    
    ###############################################################################################################################################################################
    #DOES CLEAN UP AND READIES EVERYTHING FOR THE NEXT ITERATION
    if not os.path.exists("LOG_FOLDER"):
        os.mkdir("LOG_FOLDER")
        
    if not os.path.exists("LOG_FOLDER/STEP"+str(POPULATION)):
        os.mkdir("LOG_FOLDER/STEP"+str(POPULATION))
        
    
    res = os.system('mv scfin minim* MLIP* dynamic.dat freq.dat TMP_gen'+ str(POPULATION) + ' dyn_start* dyn_end* ' + str(PREFIX) + str(POPULATION-1) + ".dyn* LOG_FOLDER/STEP"+str(POPULATION))
    
    #################################################################################################################################################################################
    
    
    
    ##################################################################################################################################################################################################
    if RELAX_FLAG == True: 
        fix_volume = False
    
        
        stress_tensor, stress_err = minimizer.get_stress_tensor()
        
        final_pressure_file = open("Pressure.dat", "a")
        final_pressure_file.write(str(POPULATION) + " ")
        for i in range(0,3):
            for j in range(0,3):
                final_pressure_file.write(str(stress_tensor[i][j]*sscha.SchaMinimizer.__RyBohr3_to_evA3__ * sscha.SchaMinimizer.__evA3_to_GPa__) + " ")
        final_pressure_file.write("\n")
        final_pressure_file.close()
                                                                                 
        
        stress_tensor *= sscha.SchaMinimizer.__RyBohr3_to_evA3__
        stress_err *=  sscha.SchaMinimizer.__RyBohr3_to_evA3__
        
        # Get the pressure
        Press = np.trace(stress_tensor) / 3
        
        
        
        
        # Get the volume
        Vol = minimizer.dyn.structure.get_volume()
        
        # Get the Helmoltz-Gibbs free energy
        target_press_evA3 = target_press / sscha.SchaMinimizer.__evA3_to_GPa__
        helmoltz = minimizer.get_free_energy() * sscha.SchaMinimizer.__RyToev__
        gibbs = helmoltz + target_press_evA3 * Vol - minimizer.eq_energy
     
        final_DATA_file = open("Final_data.dat", "a")
        final_DATA_file.write(str(POPULATION) + " " + str(N_RANDOM) + " " + str(TOTAL_NUMBER_OF_MLIP_CONF) + " " + str(T) + " " + str(target_press) + " " + str(helmoltz) + " " + str(gibbs) + " " + str(minimizer.eq_energy) + "\n")
        final_DATA_file.close()
    
    
    
        # Prepare a mark to underline which quantity is actually minimized by the
        # Variable relaxation algorithm if the helmoltz free energy (in case of fixed volume)
        # Or the Gibbs free energy (in case of fixed pressure)
        mark_helmoltz = ""
        mark_gibbs = ""
        if fix_volume:
            mark_helmoltz = "<--"
        else:
            mark_gibbs = "<--"
        
        # Extract the bulk modulus from the cell minimization
        new_bulkmodulus = 1 / (3 * BFGS.alpha * minimizer.dyn.structure.get_volume())
        new_bulkmodulus *= sscha.SchaMinimizer.__evA3_to_GPa__
        
        # Print the enthalpic contribution
        message = """
         ======================
         ENTHALPIC CONTRIBUTION
         ======================
        
         P = {:.4f} GPa   V = {:.4f} A^3
        
         P V = {:.8e} eV
        
         Helmoltz Free energy = {:.10e} eV {}
         Gibbs Free energy = {:.10e} eV {}
         Zero energy = {:.10e} eV
        
         """.format(target_press , Vol,target_press_evA3 * Vol, helmoltz, mark_helmoltz, gibbs, mark_gibbs, minimizer.eq_energy)
        print(message)
        # print " ====================== "
        # print " ENTHALPIC CONTRIBUTION "
        # print " ====================== "
        # print ""
        # print "  P = %.4f GPa    V = %.4f A^3" % (target_press , Vol)
        # print ""
        # print "  P V = %.8e eV " % (target_press_evA3 * Vol)
        # print ""
        # print " Helmoltz Free energy = %.8e eV " % helmoltz,
        # if fix_volume:
        #     print "  <-- "
        # else:
        #     print ""
        # print " Gibbs Free energy = %.8e eV " % gibbs,
        # if fix_volume:
        #     print ""
        # else:
        #     print "  <-- "
        # print " (Zero energy = %.8e eV) " % self.minim.eq_energy
        # print ""
        
        # Perform the cell step
        I = np.array([
         [ 1.000000, 0.0000000,  0.00000000 ],
         [ 0.000000, 1.0000000,  0.00000000 ],
         [ 0.000000, 0.0000000,  1.00000000 ],
        ])
        if fix_cell_shape == True:
            # Use a isotropic stress tensor to keep the same cell shape
            cell_gradient = I * (Press - target_press_evA3)
            print("pressure check:",cell_gradient,Press,target_press_evA3)
        else:
            cell_gradient = (stress_tensor - I *target_press_evA3)
            print("pressure check:",cell_gradient,Press,target_press_evA3)
        
        new_uc = minimizer.dyn.structure.unit_cell.copy()
        BFGS.UpdateCell(new_uc,  cell_gradient, fix_volume)
        
        # Strain the structure and the q points preserving the symmetries
        minimizer.dyn.AdjustToNewCell(new_uc)
        #self.minim.dyn.structure.change_unit_cell(new_uc)
        
        message = """
         Currently estimated bulk modulus = {:8.3f} GPa
         (Note: this is just indicative, do not use it for computing bulk modulus)
        
         """.format(new_bulkmodulus)
        
        
        
        print(message)
        
        
        print (" New unit cell:")
        print (" v1 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[0,0], new_uc[0,1], new_uc[0,2]))
        print (" v2 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[1,0], new_uc[1,1], new_uc[1,2]))
        print (" v3 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[2,0], new_uc[2,1], new_uc[2,2]))
        
        print ()
        print ("Check the symmetries in the new cell:")
        sys.stdout.flush()
        qe_sym = CC.symmetries.QE_Symmetry(minimizer.dyn.structure)
        qe_sym.SetupQPoint(verbose = True)
        
        print ("Forcing the symmetries in the dynamical matrix.")
        fcq = np.array(minimizer.dyn.dynmats, dtype = np.complex128)
        qe_sym.SymmetrizeFCQ(fcq, minimizer.dyn.q_stars, asr = "custom")
        for iq,q in enumerate(minimizer.dyn.q_tot):
            minimizer.dyn.dynmats[iq] = fcq[iq, :, :]
        
        
            
            
        minimizer.dyn.save_qe(PREFIX + str(POPULATION) + ".dyn")
    ##########################################################################################
    
    
    ############################
    ############################
    ############################
    POPULATION = POPULATION +1
    GENERATE_ENS = True  #RESET THE VARIABLE TO DEFAULT FOR NEXT CYCLE
    TRAINED_FLAG = False
    ############################
    ############################
    ############################




