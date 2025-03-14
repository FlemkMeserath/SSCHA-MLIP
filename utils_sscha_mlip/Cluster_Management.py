import numpy as np
import sys,os
import time
import math

class Cluster_Management:

    def __init__(self,prefix,POP=1,ENS_FOL="TMP_gen1",address="./",fol="./SSCHA", marker="SSCHA",n_rand=1):
        self.__Ry_to_eV__    = 13.6057039763
        self.__Bohr_to_A__   = 0.529177210903
        self.__e_charge__    = 1.60217663
        self.POPULATION      = POP
        self.ENS_FOLDER      = ENS_FOL
        self.ADDRESS         = address
        self.folder          = fol
        self.LOCAL           = True
        self.PREFIX          = prefix
        self.PATH            = fol
        self.MARKER          = marker
        self.TIMER           = 300
        self.N_RANDOM        = n_rand
        self.MAX             = 500 
        self.base_script_QE = ""
        self.base_script_VASP = ""

 
        
        

    def Generate_SCF(self,typical_espresso_header):
        all_scf_files = [os.path.join(self.ENS_FOLDER + str(self.POPULATION), f) for f in os.listdir(self.ENS_FOLDER + str(self.POPULATION)) if f.startswith("scf_")]
        
        
        # We will generate the input file in a new directory
        
        for file in all_scf_files:
            # Now we are cycling on the scf_ files we found.
            # We must extract the number of the file
            # The file is the string "data_ensemble_manual/scf_population1_X.dat"
            # Therefore the X number is after the last "_" and before the "." character
            # We can split before the string file at each "_", isolate the last part "X.dat"
            # and then split it again on "." (obtaining ["X", "dat"]) and select the first element
            # then we convert the "X" string into an integer
            number = int(file.split("_")[-1].split(".")[0])
            
            # We decide the filename for the espresso input
            # We will call it run_calculation/espresso_run_X.pwi
            filename = os.path.join("scfin", self.PREFIX + "{}.scf.in".format(number))
            
            # We start writing the file
            with open(filename, "w") as f:
                # We write the header
                f.write(typical_espresso_header)
                
                # Load the scf_population_X.dat file
                ff = open(file, "r")
                structure_lines = ff.readlines()
                ff.close()
                
                # Write the content on the espresso_run_X.pwi file
                # Note in the files we specify the units for both the cell and the structure [Angstrom]
                f.writelines(structure_lines) 
            f.close() 
            






    def Generate_POSCAR(self):
        all_scf_files = [os.path.join(self.ENS_FOLDER + str(self.POPULATION), f) for f in os.listdir(self.ENS_FOLDER + str(self.POPULATION)) if f.startswith("scf_")]
    
    
        # We will generate the input file in a new directory
    
        for file in all_scf_files:
            # Now we are cycling on the scf_ files we found.
            # We must extract the number of the file
            # The file is the string "data_ensemble_manual/scf_population1_X.dat"
            # Therefore the X number is after the last "_" and before the "." character
            # We can split before the string file at each "_", isolate the last part "X.dat"
            # and then split it again on "." (obtaining ["X", "dat"]) and select the first element
            # then we convert the "X" string into an integer
            number = int(file.split("_")[-1].split(".")[0])
    
            # We decide the filename for the espresso input
            # We will call it run_calculation/espresso_run_X.pwi
            filename = os.path.join("scfin", self.PREFIX + "{}.POSCAR".format(number))
    
            # We start writing the file
            with open(filename, "w") as f:
                # We write the header
    
                # Load the scf_population_X.dat file
                ff = open(file, "r")
                structure_lines = ff.readlines()
                ff.close()
    
                file_lenght = len(structure_lines)
                f.write("\n")
                f.write("1.0000\n")
                for i in range(0,3): f.write(structure_lines[i+1])
    
                atm_t = []
                atm_n = []
                atom_ordering = [None]*(file_lenght-6)
                atm_t.append(structure_lines[6].split()[0])
                atm_n.append(1)
    
                for i in range(7,file_lenght):
                    flag = True
                    for j in range(0,len(atm_t)):
                        if atm_t[j] == structure_lines[i].split()[0]:
                            atm_n[j] = atm_n[j] + 1
                            flag = False
                    if flag == True:
                        atm_t.append(structure_lines[i].split()[0])
                        atm_n.append(1)
                for i in atm_t:
                    f.write(i + str("  "))
                f.write("\n")
                for i in atm_n:
                    f.write(str(i) + str("  "))
                f.write("\nCartesian\n")
    
    
    
                count_n = 0
                for i in atm_t:
                    for j in range(6,file_lenght):
                        if i == structure_lines[j].split()[0]:
                            f.write(structure_lines[j].split()[1] + str("  ") + structure_lines[j].split()[2] + str("  ") + structure_lines[j].split()[3] + str("\n"))
                            atom_ordering[j-6] = count_n
                            count_n += 1
    
                # Write the content on the espresso_run_X.pwi file
                # Note in the files we specify the units for both the cell and the structure [Angstrom]
            f.close()
        return atom_ordering
        
        
        
        
        
        
        
        
    def Send_to_Folders(self,conf_table_gamma):
        current_dir = os.getcwd()
        res = os.system("> scfin/list.dat")
        for i in range(0,len(conf_table_gamma)):
            tmp_line = "echo '#" +str(i+1) + "# " + self.PREFIX + str(conf_table_gamma[i][0]) + "' >> scfin/list.dat"
            res = os.system(tmp_line)   
    
        Submitter_file = open("Submitter.sh", "w")
        if int(self.MAX) > int(len(conf_table_gamma)):
            Submitter_file.write('#!/bin/bash -l\n\nmax=' + str(len(conf_table_gamma)) + '\ntop=' + str(len(conf_table_gamma)) +  '\nPREFNAME="' + self.MARKER + '" \n')
        else:
            Submitter_file.write('#!/bin/bash -l\n\nmax=' + str(self.MAX) + '\ntop=' + str(len(conf_table_gamma)) +  '\nPREFNAME="' + self.MARKER + '" \n')    
        Submitter_file.write(self.base_script_QE)
        Submitter_file.close()
        res = os.system('chmod +x Submitter.sh')
            
        if self.LOCAL == False:    
            res = os.system('scp -q scfin/list.dat Submitter.sh ' + self.ADDRESS + ':"' + self.PATH + '/"')
            res = os.system('scp -q scfin/' + self.PREFIX +'*.in ' + self.ADDRESS + ':"' + self.PATH + '/Ins"')
            res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; ./Submitter.sh"')
        else:
            res = os.system('cp  scfin/list.dat Submitter.sh ' + self.PATH + '/')
            res = os.system('cp  scfin/' + self.PREFIX + '*.in ' + self.PATH + '/Ins')
            res = os.system('cd ' + current_dir)
            print("submitting into send to folder")
            res = os.system('cd ' + self.PATH + '/; ./Submitter.sh')        
            res = os.system('cd ' + current_dir)
        res = os.system('rm Submitter.sh')
    





    
    def Send_to_Folders_VASP(self,conf_table_gamma):
        current_dir = os.getcwd()
        res = os.system("> scfin/list.dat")
        for i in range(0,len(conf_table_gamma)):
            tmp_line = "echo '#" +str(i+1) + "# " + self.PREFIX + str(conf_table_gamma[i][0]) + "' >> scfin/list.dat"
            res = os.system(tmp_line)
    
        Submitter_file = open("Submitter.sh", "w")
        if int(self.MAX) > int(len(conf_table_gamma)):
            Submitter_file.write('#!/bin/bash -l\n\nmax=' + str(len(conf_table_gamma)) + '\ntop=' + str(len(conf_table_gamma)) +  '\nPREFNAME="' + self.MARKER + '" \n')
        else:
            Submitter_file.write('#!/bin/bash -l\n\nmax=' + str(self.MAX) + '\ntop=' + str(len(conf_table_gamma)) +  '\nPREFNAME="' + self.MARKER + '" \n')
        Submitter_file.write(self.base_script_VASP)
        Submitter_file.close()
        res = os.system('chmod +x Submitter.sh')
    
        if self.LOCAL == False:
            res = os.system('scp -q scfin/list.dat Submitter.sh ' + self.ADDRESS + ':"' + self.PATH + '/"')
            res = os.system('scp -q scfin/' + self.PREFIX +'*.POSCAR ' + self.ADDRESS + ':"' + self.PATH + '/Ins"')
            res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; ./Submitter.sh"')
        else:
            res = os.system('cp  scfin/list.dat Submitter.sh ' + self.PATH + '/')
            res = os.system('cp  scfin/' + self.PREFIX + '*.POSCAR ' + self.PATH + '/Ins')
            res = os.system('cd ' + current_dir)
            print("submitting into send to folder")
            res = os.system('cd ' + self.PATH + '/; ./Submitter.sh')
            res = os.system('cd ' + current_dir)
        res = os.system('rm Submitter.sh')     
        



        
        
        
    def Queue_and_resubmit(self,conf_table_gamma):
        current_dir = os.getcwd()
        p=True
       
        while(p==True):
            
            print("Starting cycle.")
            print("Downloading queue info.")
            if self.LOCAL == False:
                res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; ./queue_write.sh"')
            else:
                res = os.system('cd ' + current_dir) 
                res = os.system('cd ' + self.PATH + '/; ./queue_write.sh')
                res = os.system('cd ' + current_dir) 
            print("Downloading output files.")
            if self.LOCAL == False:
    
                res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/Outs; tar -czf archive.tar.gz *.out"')
                res = os.system('scp -q ' + self.ADDRESS + ':"' + self.PATH + '/Outs/archive.tar.gz" ./scfin/' )
                res = os.system('scp -q ' + self.ADDRESS + ':"' + self.PATH + '/queue.tmp" ./')
    
            else:
                res = os.system('cd ' + current_dir) 
                res = os.system('cd ' + self.PATH + '/Outs; tar -czf archive.tar.gz *.out')
                res = os.system('cd ' + current_dir) 
                res = os.system('cp ' + self.PATH + '/Outs/archive.tar.gz ./scfin/' )
                res = os.system('cp ' + self.PATH + '/queue.tmp ./')
                        
            res = os.system('tar -xf ./scfin/archive.tar.gz -C ./scfin/')       
    
            queue_file = open("queue.tmp","r")
            queue = queue_file.read()
            print(queue)
            still_in_queue = queue.find(self.MARKER)
            
            if still_in_queue != -1:
                print("Still in queue.")
                res = os.system('rm queue.tmp')
                res = os.system('rm ./scfin/archive.tar.gz')
                print("Cleaning on cluster.")
    
                if self.LOCAL == False:
                    res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; rm queue.tmp; rm Outs/archive.tar.gz "')
                else:
                    res = os.system('cd '+ self.PATH + '/; rm queue.tmp; rm Outs/archive.tar.gz')           
    
                ttt = time.localtime()
                print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(self.TIMER/60) + " minutes (" + str(self.TIMER/3600) + " hours)")
                time.sleep(self.TIMER)
                
            else:   
                res = os.system('rm queue.tmp')
                res = os.system('rm ./scfin/archive.tar.gz')
                print("Cleaning up on cluster.")
                if self.LOCAL == False:
                    res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; rm queue.tmp; rm Outs/archive.tar.gz "')
                else:
                    res = os.system('cd ' + current_dir)
                    res = os.system('cd '+ self.PATH + '/; rm queue.tmp; rm Outs/archive.tar.gz ')
                    res = os.system('cd ' + current_dir)
                
                j=0
                list_file = open("./scfin/list.dat","w")
                for i in range(1,len(conf_table_gamma)+1):
                    filename = "./scfin/" + self.PREFIX + str(conf_table_gamma[i-1][0]) + ".scf.out"
                    try:
                        open_output = open(filename, 'r')
                        output_file_stored = open_output.read()
                        check_var = output_file_stored.find("JOB DONE")
                        open_output.close()
                    except:
                        check_var = -1
                        
                    if check_var == -1:
                        j = j+1
                        list_file.write("#"+str(j)+"# "+ self.PREFIX + str(conf_table_gamma[i-1][0]) + "\n")
                list_file.close()
                
                if j>0:
                    print("Resubmitting jobs.")
                    Submitter_in = open("Submitter.sh","w")
                    if j < self.MAX:
                        Submitter_in.write('#!/bin/bash -l\n\nmax=' + str(j) + '\ntop=' + str(j) +  '\nPREFNAME="' + self.MARKER + '"\n')
                    else:
                        Submitter_in.write('#!/bin/bash -l\n\nmax=' + str(self.MAX) + '\ntop=' + str(j) +  '\nPREFNAME="' + self.MARKER + '"\n')
                    Submitter_in.write(self.base_script_QE)
                    Submitter_in.close()
                    res = os.system('chmod +x Submitter.sh')

                    if self.LOCAL == False:    
                        print('scp -q scfin/' + self.PREFIX + '*.in ' + self.ADDRESS + ':"' + self.PATH + '/Ins"')
                        res = os.system('scp -q scfin/list.dat Submitter.sh ' + self.ADDRESS + ':"' + self.PATH + '/"')
                        res = os.system('scp -q scfin/' + self.PREFIX + '*.in ' + self.ADDRESS + ':"' + self.PATH + '/Ins"')
                        res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; ./Submitter.sh"')
                    else:
                        res = os.system('cp  scfin/list.dat Submitter.sh ' + self.PATH + '/')
                        res = os.system('cp  scfin/' + self.PREFIX +'*.in ' + self.PATH + '/Ins')
                        res = os.system('cd ' + current_dir) 
                        res = os.system('cd ' + self.PATH + '/; ./Submitter.sh') 
                        res = os.system('cd ' + current_dir)       
                    res = os.system('rm Submitter.sh')
                else:
                    p=False
                if j > 0:
                    ttt = time.localtime()
                    print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(self.TIMER/60) + " minutes (" + str(self.TIMER/3600) + " hours)")
                    time.sleep(self.TIMER)
        if self.LOCAL == False:
            res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; rm Job* *.out; rm Ins/*; rm Outs/*"') 
        else:
            res = os.system('rm '+ self.PATH + '/Job*; rm '+ self.PATH + '/*.out; rm '+ self.PATH + '/Ins/*; rm '+ self.PATH + '/Outs/*')
            
    
    

    
    def Queue_and_resubmit_VASP(self,conf_table_gamma):
        current_dir = os.getcwd()
        p=True
    
        while(p==True):
    
            print("Starting cycle.")
            print("Downloading queue info.")
            if self.LOCAL == False:
                res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; ./queue_write.sh"')
            else:
                res = os.system('cd ' + current_dir)
                res = os.system('cd ' + self.PATH + '/; ./queue_write.sh')
                res = os.system('cd ' + current_dir)
            print("Downloading output files.")
            if self.LOCAL == False:
    
                res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/Outs; tar -czf archive.tar.gz *.OUTCAR"')
                res = os.system('scp -q ' + self.ADDRESS + ':"' + self.PATH + '/Outs/archive.tar.gz" ./scfin/' )
                res = os.system('scp -q ' + self.ADDRESS + ':"' + self.PATH + '/queue.tmp" ./')
    
            else:
                res = os.system('cd ' + current_dir)
                res = os.system('cd ' + self.PATH + '/Outs; tar -czf archive.tar.gz *.OUTCAR')
                res = os.system('cd ' + current_dir)
                res = os.system('cp ' + self.PATH + '/Outs/archive.tar.gz ./scfin/' )
                res = os.system('cp ' + self.PATH + '/queue.tmp ./')
    
            res = os.system('tar -xf ./scfin/archive.tar.gz -C ./scfin/')
    
            queue_file = open("queue.tmp","r")
            queue = queue_file.read()
            print(queue)
            still_in_queue = queue.find(self.MARKER)
    
            if still_in_queue != -1:
                print("Still in queue.")
                res = os.system('rm queue.tmp')
                res = os.system('rm ./scfin/archive.tar.gz')
                print("Cleaning on cluster.")
    
                if self.LOCAL == False:
                    res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; rm queue.tmp; rm Outs/archive.tar.gz "')
                else:
                    res = os.system('cd '+ self.PATH + '/; rm queue.tmp; rm Outs/archive.tar.gz')
    
                ttt = time.localtime()
                print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(self.TIMER/60) + " minutes (" + str(self.TIMER/3600) + " hours)")
                time.sleep(self.TIMER)
    
            else:
                res = os.system('rm queue.tmp')
                res = os.system('rm ./scfin/archive.tar.gz')
                print("Cleaning up on cluster.")
                if self.LOCAL == False:
                    res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; rm queue.tmp; rm Outs/archive.tar.gz "')
                else:
                    res = os.system('cd ' + current_dir)
                    res = os.system('cd '+ self.PATH + '/; rm queue.tmp; rm Outs/archive.tar.gz ')
                    res = os.system('cd ' + current_dir)
    
                j=0
                list_file = open("./scfin/list.dat","w")
                for i in range(1,len(conf_table_gamma)+1):
                    filename = "./scfin/" + self.PREFIX + str(conf_table_gamma[i-1][0]) + ".OUTCAR"
                    try:
                        open_output = open(filename, 'r')
                        output_file_stored = open_output.read()
                        check_var = output_file_stored.find("General timing and accounting informations for this job:")
                        open_output.close()
                    except:
                        check_var = -1
    
                    if check_var == -1:
                        j = j+1
                        list_file.write("#"+str(j)+"# "+ self.PREFIX + str(conf_table_gamma[i-1][0]) + "\n")
                list_file.close()
    
                if j>0:
                    print("Resubmitting jobs.")
                    Submitter_in = open("Submitter.sh","w")
                    if j < self.MAX:
                        Submitter_in.write('#!/bin/bash -l\n\nmax=' + str(j) + '\ntop=' + str(j) +  '\nPREFNAME="' + self.MARKER + '"\n')
                    else:
                        Submitter_in.write('#!/bin/bash -l\n\nmax=' + str(self.MAX) + '\ntop=' + str(j) +  '\nPREFNAME="' + self.MARKER + '"\n')
                    Submitter_in.write(self.base_script_VASP)
                    Submitter_in.close()
                    res = os.system('chmod +x Submitter.sh')

                    if self.LOCAL == False:
                        print('scp -q scfin/' + self.PREFIX + '*.in ' + self.ADDRESS + ':"' + self.PATH + '/Ins"')
                        res = os.system('scp -q scfin/list.dat Submitter.sh ' + self.ADDRESS + ':"' + self.PATH + '/"')
                        res = os.system('scp -q scfin/' + self.PREFIX + '*.POSCAR ' + self.ADDRESS + ':"' + self.PATH + '/Ins"')
                        res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; ./Submitter.sh"')
                    else:
                        res = os.system('cp  scfin/list.dat Submitter.sh ' + self.PATH + '/')
                        res = os.system('cp  scfin/' + self.PREFIX +'*.POSCAR ' + self.PATH + '/Ins')
                        res = os.system('cd ' + current_dir)
                        res = os.system('cd ' + self.PATH + '/; ./Submitter.sh')
                        res = os.system('cd ' + current_dir)
                    res = os.system('rm Submitter.sh')
                else:
                    p=False
                if j > 0:
                    ttt = time.localtime()
                    print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(self.TIMER/60) + " minutes (" + str(self.TIMER/3600) + " hours)")
                    time.sleep(self.TIMER)
        if self.LOCAL == False:
            res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+ self.PATH + '/; rm Job* *.out; rm Ins/*; rm Outs/*"')
        else:
            res = os.system('rm '+ self.PATH + '/Job*; rm '+ self.PATH + '/*.out; rm '+ self.PATH + '/Ins/*; rm '+ self.PATH + '/Outs/*')
            





    def read_SCF(self):
        directory = "scfin"
        output_filenames = [f for f in os.listdir(directory) if f.endswith(".out")] # We select only the output files
        output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly
        energies = np.zeros(self.N_RANDOM)
        
        for file in output_files:
            # Get the number of the configuration.
            id_number = int(file.split(".")[0].split(self.PREFIX)[-1]) # The same as before, we need the to extract the configuration number from the filename
            # Load the file
            ff = open(file, "r")
            lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
            ff.close()
            
            # Lets look for the energy (in espresso the first line that starts with !)
            # next is used to find only the first occurrence
            energy_line = next(l for l in lines if len(l) > 0 if l.split()[0] == "!")
            
            # Lets collect the energy (the actual number is the 5th item on the line, but python indexes start from 0)
            # note, also the id_number are saved starting from 1
            energies[id_number - 1] = float(energy_line.split()[4])
            
            # Now we can collect the force
            # We need the number of atoms
            nat_line = next( l for l in lines if len(l) > 0 if l.split()[0] == "number" and l.split()[2] == "atoms/cell" )
            nat = int(nat_line.split()[4])
            
            # Now allocate the forces and read them
            forces = np.zeros((nat, 3))
            forces_lines = [l for l in lines if len(l) > 0 if l.split()[0] == "atom"] # All the lines that starts with atom will contain a force
            for i in range(nat):
                forces[i, :] = [float(x) for x in forces_lines[i].split()[-3:]] # Get the last three number from the line containing the force
            
            # Now we can take the stress tensor
            stress = np.zeros((3,3))
            # We pick the index of the line that starts with the words total stress
            index_before_stress = next(i for i, l in enumerate(lines) if len(l) > 0 if l.split()[0] == "total" and l.split()[1] == "stress")
            # The stress tensor is located just after it
            for i in range(3):
                index = i + index_before_stress + 1
                stress[i, :] = [float(x) for x in lines[index].split()[:3]]
        
            # We can save the forces_population1_X.dat and pressures_population1_X.dat files
            force_file = os.path.join( self.ENS_FOLDER + str(self.POPULATION), "forces_population"+ str(self.POPULATION) +"_{}.dat".format(id_number))
            stress_file = os.path.join( self.ENS_FOLDER + str(self.POPULATION), "pressures_population" + str(self.POPULATION) + "_{}.dat".format(id_number))
            np.savetxt(force_file, forces)
            np.savetxt(stress_file, stress)
        
        # Now we read all the configurations, we can save the energy file
        energy_file = os.path.join(self.ENS_FOLDER + str(self.POPULATION), "energies_supercell_population" + str(self.POPULATION) + ".dat")
        np.savetxt(energy_file, energies)
        return energies
    
    
    
  
  
  
  
    
    
    
    
    def read_OUTCAR(self,saved_ordering):
        directory = "scfin"
        output_filenames = [f for f in os.listdir(directory) if f.endswith("OUTCAR")] # We select only the output files
        output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly
        energies = np.zeros(self.N_RANDOM)
    
        for file in output_files:
            nat = 0
            # Get the number of the configuration.
            id_number = int(file.split(".")[0].split(self.PREFIX)[-1]) # The same as before, we need the to extract the configuration number from the filename
            # Load the file
            ff = open(file, "r")
            lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
            ff.close()
    
            # Lets look for the energy (in espresso the first line that starts with !)
            # next is used to find only the first occurrence
            energy_line = next(l for l in lines if len(l) > 0 if l.find("free  energy   TOTEN ") != -1)
    
            # Lets collect the energy (the actual number is the 5th item on the line, but python indexes start from 0)
            # note, also the id_number are saved starting from 1
            energies[id_number - 1] = float(energy_line.split()[4])/self.__Ry_to_eV__
    
            # Now we can collect the force
            # We need the number of atoms
            nat_line = next( l for l in lines if len(l) > 0 if l.split()[0] == "ions" and l.split()[2] == "type" )
    
    
            for i in range(4, len(nat_line.split())):
                nat = nat + int( nat_line.split()[i] )
    
            # Now allocate the forces and read them
            forces = np.zeros((nat, 3))
            flag = False
            force_lines = []
            for i in range(0,len(lines)):
                if lines[i].find("total drift:") != -1 and flag == True : flag = False
                if flag == True and len(lines[i].split()) > 1: force_lines.append(lines[i])
                if lines[i].find("TOTAL-FORCE (eV/Angst)") != -1 and flag == False : flag = True
                if len(lines[i].split())  > 2:
                    if lines[i].split()[0] == "in" and lines[i].split()[1] == "kB": index_before_stress = lines[i]
    
    
            for i in range(nat):
                forces[i, :] = [float(x)/self.__Ry_to_eV__*self.__Bohr_to_A__ for x in force_lines[saved_ordering[i]].split()[-3:]] # Get the last three number from the line containing the force
    
    
            # Now we can take the stress tensor
            stress = np.zeros((3,3))
            # We pick the index of the line that starts with the words total stress
            # The stress tensor is located just after it

            stress[0,0] = float(index_before_stress.split()[2])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
            stress[1,1] = float(index_before_stress.split()[3])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
            stress[2,2] = float(index_before_stress.split()[4])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
            stress[2,1] = float(index_before_stress.split()[5])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
            stress[1,2] = float(index_before_stress.split()[5])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
            stress[0,2] = float(index_before_stress.split()[6])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
            stress[2,0] = float(index_before_stress.split()[6])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
            stress[1,0] = float(index_before_stress.split()[7])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
            stress[0,1] = float(index_before_stress.split()[7])/10*math.pow(10,9)/math.pow(10,30)/self.__e_charge__*math.pow(10,19)/self.__Ry_to_eV__*math.pow(self.__Bohr_to_A__,3)
    
            # We can save the forces_population1_X.dat and pressures_population1_X.dat files
            force_file = os.path.join( self.ENS_FOLDER + str(self.POPULATION), "forces_population"+ str(self.POPULATION) +"_{}.dat".format(id_number))
            stress_file = os.path.join( self.ENS_FOLDER + str(self.POPULATION), "pressures_population" + str(self.POPULATION) + "_{}.dat".format(id_number))
            np.savetxt(force_file, forces)
            np.savetxt(stress_file, stress)
    
        # Now we read all the configurations, we can save the energy file
        energy_file = os.path.join(self.ENS_FOLDER + str(self.POPULATION), "energies_supercell_population" + str(self.POPULATION) + ".dat")
        np.savetxt(energy_file, energies)
        return energies

           
