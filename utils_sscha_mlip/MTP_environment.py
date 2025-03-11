import numpy as np
import sys, os
import time


class MTP_Environment:

    def __init__(self,POP=1,ENS_FOL="TMP_gen1",G=1,address="./",fol="./SSCHA",exe_path="./"):
        self.__Ry_to_eV__    = 13.6057039763
        self.__Bohr_to_A__   = 0.529177210903
        self.POPULATION      = POP
        self.ENS_FOLDER      = ENS_FOL
        self.atom_type_table = {}
        self.MLIP_PATH       = exe_path 
        self.PRETRAINED      = False
        self.GAMMA           = G
        self.conf_table_gamma = []
        self.ADDRESS         = address
        self.folder          = fol
        self.LOCAL           = True


    ###################################################    
    ###################################################
    ###################################################
    def Generate_CFG(self):
    # Generates the CFG files used by the MTP 
        cfg_file = open("MLIP_"+str(self.POPULATION) + ".cfg", "w")
        directory = self.ENS_FOLDER + str(self.POPULATION)
        FORCE_FILES = True
    
        energy_file = open(directory + "/energies_supercell_population" + str(self.POPULATION) + ".dat", "r")
        energy_lines = [l.strip() for l in energy_file.readlines()]
        energy_file.close()
    
        print(len(energy_lines),energy_lines)
    
        for j in range(1,len(energy_lines)+1):
            structure_file = open(directory + "/scf_population" + str(self.POPULATION) + "_" + str(j) + ".dat", "r")
            lines = [l.strip() for l in structure_file.readlines()]
            structure_file.close()
    
            try:
                force_file = open(directory + "/forces_population" + str(self.POPULATION) + "_" + str(j) + ".dat", "r")
                force_lines = [l.strip() for l in force_file.readlines()]
                force_file.close()
            except: FORCE_FILES = False
    
            try:
                pressure_file = open(directory + "/pressures_population" + str(self.POPULATION) + "_" + str(j) + ".dat", "r")
                pressure_lines = [l.strip() for l in pressure_file.readlines()]
                pressure_file.close()
            except:
                pressure_lines = [
                "0.0000 0.0000 0.0000",
                "0.0000 0.0000 0.0000",
                "0.0000 0.0000 0.0000",
                ]
    
            SIZE = len(lines)-6
            cell = np.zeros((3, 3))
            atoms = np.zeros((SIZE,3))
            forces = np.zeros((SIZE,3))
            atm_type = [None] * SIZE
            pressures = np.zeros((3, 3))
    
            for i in range(0,3):
                cell[i, :] = [float(x) for x in lines[i+1].split()[-3:]]
    
            for i in range(0,SIZE):
                atoms[i, :] = [float(x) for x in lines[i+6].split()[-3:]]
                atm_type[i] =  lines[i+6].split()[0]
                if FORCE_FILES == True:
                    forces[i,:] = [float(x) for x in force_lines[i].split()[-3:]]
                    forces[i,:] = forces[i, :] * self.__Ry_to_eV__ / self.__Bohr_to_A__
    
    
            for i in atm_type:
                try: self.atom_type_table[i]
                except: self.atom_type_table[i] = len(self.atom_type_table)
    
    
            Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
    
            for i in range(0,3):
    
                pressures[i, :] = [float(x) for x in pressure_lines[i].split()[-3:]]
                pressures[i, :] = pressures[i, :] * self.__Ry_to_eV__ / self.__Bohr_to_A__ / self.__Bohr_to_A__ / self.__Bohr_to_A__* Volume
    
            cfg_file.write("BEGIN_CFG\n")
            cfg_file.write(" Size\n")
            cfg_file.write("{: >5}".format(SIZE) +"\n")
            cfg_file.write(" Supercell\n")
            for row in cell:
                cfg_file.write("    {: >13f} {: >13f} {: >13f}\n".format(*row))
            cfg_file.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
            for i in range(0,len(atoms)):
                if FORCE_FILES == True:
                    cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(self.atom_type_table[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*forces[i,:]))
                else:
                    cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(self.atom_type_table[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*[0,0,0]))
            cfg_file.write(" Energy\n")
            cfg_file.write("     {: >13f}".format(float(energy_lines[j-1])*self.__Ry_to_eV__) +"\n")
            cfg_file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
            cfg_file.write("    {: >12f}".format(pressures[0,0]) + "{: >12f}".format(pressures[1,1]) +  "{: >12f}".format(pressures[2,2]))
            cfg_file.write("{: >12f}".format(pressures[1,2]) + "{: >12f}".format(pressures[0,2]) +  "{: >12f}\n".format(pressures[0,1]))
            cfg_file.write(" Feature atom_type_table " + str(self.atom_type_table) + "\n")
            cfg_file.write(" Feature conf_number " + str(j) + "\n")
            cfg_file.write(" Feature population " + str(self.POPULATION) + "\n")
            cfg_file.write("END_CFG\n\n")
        cfg_file.close()
        
        




    ###################################################    
    ###################################################
    ###################################################        
    def Calc_GRADE(self):
    #calculates the grade for the generation of the trainingset
    
        res = os.system(self.MLIP_PATH + ' mindist MLIP_' + str(self.POPULATION) + '.cfg')
        if self.POPULATION == 1 and self.PRETRAINED == False:
            res = os.system(self.MLIP_PATH + ' calc-grade potential.mtp MLIP_'+str(self.POPULATION)+'.cfg MLIP_'+str(self.POPULATION)+'.cfg MLIP_'+str(self.POPULATION)+'.gamma.cfg')
        else:
            res = os.system(self.MLIP_PATH + ' calc-grade potential.mtp trainset.cfg MLIP_'+str(self.POPULATION)+'.cfg MLIP_'+str(self.POPULATION)+'.gamma.cfg')        
        
        
    ###################################################    
    ###################################################
    ###################################################          
    def Fill_GAMMA_Table_0(self,N_RANDOM):
        self.conf_table_gamma = []
        for i in range(0,N_RANDOM):
            self.conf_table_gamma.append([i+1,1])        
        
 
        

    ###################################################    
    ###################################################
    ###################################################     
    def Fill_GAMMA_Table(self):

       #manages the gamma table for the analysis of the trainingset
        self.conf_table_gamma = []
        conf_output_file = open( "MLIP_"+str(self.POPULATION) + ".gamma.cfg", "r")
                
        for line in conf_output_file:
        #PRINT OUTPUTS CLOSES FILES 
            if line.find("END_CFG") != -1:
                if gamma_conf >= self.GAMMA or (self.POPULATION == 1 and self.PRETRAINED == False):# or abs(float(energy_lines[conf_number-1])) > 0.0001 :
                    self.conf_table_gamma.append([conf_number,gamma_conf])    
        #GETS THE GAMMA FOR ACTIVE LEARNING
            if line.find("Feature   MV_grade") != -1:
                gamma_conf = float(line.split()[-1])
         
            if line.find("Feature   conf_number") != -1:
                conf_number = int(line.split()[-1])   
            
        print(len(self.conf_table_gamma),self.conf_table_gamma) 
        conf_output_file.close()      
            
        

    ###################################################    
    ###################################################
    ###################################################
    def Compile_Trainingset(self):

    
        cfg_file = open("trainset.cfg", "a")
        directory = self.ENS_FOLDER + str(self.POPULATION)
        
        energy_file = open(directory + "/energies_supercell_population" + str(self.POPULATION) + ".dat", "r")
        energy_lines = [l.strip() for l in energy_file.readlines()]
        energy_file.close()
        
        print(len(energy_lines))
        
        for j in range(1,len(self.conf_table_gamma)+1):
            structure_file = open(directory + "/scf_population" + str(self.POPULATION) + "_" + str(self.conf_table_gamma[j-1][0]) + ".dat", "r")
            lines = [l.strip() for l in structure_file.readlines()]
            
            force_file = open(directory + "/forces_population" + str(self.POPULATION) + "_" + str(self.conf_table_gamma[j-1][0]) + ".dat", "r")
            force_lines = [l.strip() for l in force_file.readlines()]
            force_file.close()
        
            pressure_file = open(directory + "/pressures_population" + str(self.POPULATION) + "_" + str(self.conf_table_gamma[j-1][0]) + ".dat", "r")
            pressure_lines = [l.strip() for l in pressure_file.readlines()]
            pressure_file.close()
            
            SIZE = len(lines)-6
            cell = np.zeros((3, 3))
            atoms = np.zeros((SIZE,3))
            forces = np.zeros((SIZE,3))
            atm_type = [None] * SIZE
            pressures = np.zeros((3, 3))
            
            for i in range(0,3):
                cell[i, :] = [float(x) for x in lines[i+1].split()[-3:]]
                
            for i in range(0,SIZE):
                atoms[i, :] = [float(x) for x in lines[i+6].split()[-3:]]
                atm_type[i] =  lines[i+6].split()[0]
        
                forces[i,:] = [float(x) for x in force_lines[i].split()[-3:]]
                forces[i,:] = forces[i, :] * self.__Ry_to_eV__ / self.__Bohr_to_A__
                
                
            for i in atm_type:
                try: self.atom_type_table[i]
                except: self.atom_type_table[i] = len(self.atom_type_table)    
            
            
            Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
            
            
            for i in range(0,3):
                pressures[i, :] = [float(x) for x in pressure_lines[i].split()[-3:]]    
                pressures[i, :] = pressures[i, :] * self.__Ry_to_eV__ / self.__Bohr_to_A__ / self.__Bohr_to_A__ / self.__Bohr_to_A__* Volume
            
            cfg_file.write("BEGIN_CFG\n")
            cfg_file.write(" Size\n")
            cfg_file.write("{: >5}".format(SIZE) +"\n")
            cfg_file.write(" Supercell\n")
            for row in cell:
                cfg_file.write("    {: >13f} {: >13f} {: >13f}\n".format(*row))
            cfg_file.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
            for i in range(0,len(atoms)):
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(self.atom_type_table[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*forces[i,:]))
                
            cfg_file.write(" Energy\n")
            cfg_file.write("     {: >13f}".format(float(energy_lines[self.conf_table_gamma[j-1][0]-1])*self.__Ry_to_eV__) +"\n")
            cfg_file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
            cfg_file.write("    {: >12f}".format(pressures[0,0]) + "{: >12f}".format(pressures[1,1]) +  "{: >12f}".format(pressures[2,2]))
            cfg_file.write("{: >12f}".format(pressures[1,2]) + "{: >12f}".format(pressures[0,2]) +  "{: >12f}\n".format(pressures[0,1]))
            cfg_file.write(" Feature atom_type_table " + str(self.atom_type_table) + "\n")
            cfg_file.write(" Feature conf_number " + str(self.conf_table_gamma[j-1][0]) + "\n")
            cfg_file.write(" Feature population " + str(self.POPULATION) + "\n")
            cfg_file.write("END_CFG\n\n")
        cfg_file.close()            
        
        
        


        
    def Train_submit(self):
        current_dir = os.getcwd()
        if self.LOCAL == False:
            res = os.system('scp -q potential.mtp trainset.cfg MLIP_' + str(self.POPULATION) + '.cfg ' + self.ADDRESS + ':"' +   self.folder + '/"' )
        else:
            res = os.system('cp  potential.mtp trainset.cfg MLIP_' + str(self.POPULATION) + '.cfg ' +    self.folder + '/' )
            
        if len(self.conf_table_gamma) >= 1:
            print("training")
            if self.LOCAL == False:
                res = os.system('ssh -q ' + self.ADDRESS + ' "cd ' +  self.folder + '/; ./SEND_TRAINING.sh"')
            else:
                res = os.system('cd ' + current_dir) 
                res = os.system('cd ' +  self.folder + '/; ./SEND_TRAINING.sh') 
                res = os.system('cd ' + current_dir)  
            p=True  
            while(p==True):
                if self.LOCAL == False:
                    res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+  self.folder + '/;  ./queue_write.sh"')
                    res = os.system('scp -q ' + self.ADDRESS + ':"' +  self.folder + '/queue.tmp" ./')
                    res = os.system('ssh -q ' + self.ADDRESS + ' "cd '+  self.folder + '/; rm queue.tmp"')
                else:
                    res = os.system('cd ' + current_dir) 
                    res = os.system('cd ' +  self.folder + '/;  ./queue_write.sh')
                    res = os.system('cd ' + current_dir) 
                    res = os.system('cp ' +  self.folder + '/queue.tmp ./')
                    res = os.system('rm ' +  self.folder + '/queue.tmp')            
                queue_file = open("queue.tmp","r")
                queue = queue_file.read()
                still_in_queue = queue.find("test_ta")
                queue_file.close()
                if still_in_queue != -1:
                    print("Still in queue.")
                    res = os.system('rm queue.tmp')
                    ttt = time.localtime()
                    print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(600/60) + " minutes (" + str(100/3600) + " hours)")
                    time.sleep(100)                
                else:
                    p=False
        if self.LOCAL == False:
            res = os.system('scp -q ' + self.ADDRESS + ':"'+  self.folder + '/potential.mtp" ./')
        else:
            res = os.system('cd ' + current_dir)
            res = os.system('cp ' +  self.folder + '/potential.mtp ./')
        res = os.system(self.MLIP_PATH + ' calc-efs potential.mtp MLIP_'+str(self.POPULATION)+'.cfg MLIP_'+str(self.POPULATION)+'.out.cfg')          
            
        



    def read_CFG(self,energies):
        
        directory = self.ENS_FOLDER + str(self.POPULATION)
        cfg_file = "MLIP_" + str(self.POPULATION) + ".out.cfg"
        
        mlip_conf_list = []
        
        for i in range(0,len(self.conf_table_gamma)):
            mlip_conf_list.append(self.conf_table_gamma[i][0])
        
        
        conf_output_file = open( cfg_file, "r")
        energy_file = open(directory + "/energies_supercell_population" + str(self.POPULATION) + ".dat", "w")
        
        ITEM=0
        FORCE_FLAG = False
        ENERGY_FLAG = False
        PRESSURE_FLAG = False
        CELL_FLAG = False
        
        #Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
        Volume = 1
        
        print(energies)
        print(mlip_conf_list)
        
        TOTAL_NUMBER_OF_MLIP_CONF = 0
        for line in conf_output_file:
        #OPEN FILES ZEROS STORAGE ######################################################################################
            if line.find("BEGIN_CFG") != -1:
                ITEM += 1
                cell = np.zeros((3, 3))
                cell_tmp = []
                forces = []
                pressure = [[0,0,0],[0,0,0],[0,0,0]]
                energy = 0
        #PRINT OUTPUTS CLOSES FILES ######################################################################################
            if line.find("END_CFG") != -1:
                
                
                if conf_number in mlip_conf_list or abs(float(energies[conf_number-1])) > 0.0001:
                    
                    print("reading " + str(conf_number) + " from DFT" )
                    energy_file.write(str(energies[conf_number-1]) + "\n")
                    
                else:   
                    
                    print("reading " + str(ITEM) + " from MLIP" )
                    TOTAL_NUMBER_OF_MLIP_CONF = TOTAL_NUMBER_OF_MLIP_CONF +1
                    
                    pressure_file = open(directory + "/pressures_population" + str(self.POPULATION) + "_" + str(ITEM) + ".dat", "w")
                    force_file = open(directory + "/forces_population" + str(self.POPULATION) + "_" + str(ITEM) + ".dat", "w")
                    
                    energy_file.write(str(energy) + "\n")
                    
                    for data_to_print in forces:
                        force_file.write("{} {} {}".format(data_to_print[0],data_to_print[1],data_to_print[2]) + "\n")
        
                    cell[0][:] = cell_tmp[0][:]
                    cell[1][:] = cell_tmp[1][:]
                    cell[2][:] = cell_tmp[2][:]
                    Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
                
                
                    pressure_file.write("{} {} {}".format(pressure[0][0]/Volume,pressure[0][1]/Volume,pressure[0][2]/Volume) + "\n")
                    pressure_file.write("{} {} {}".format(pressure[1][0]/Volume,pressure[1][1]/Volume,pressure[1][2]/Volume) + "\n")
                    pressure_file.write("{} {} {}".format(pressure[2][0]/Volume,pressure[2][1]/Volume,pressure[2][2]/Volume) + "\n")
              
                    pressure_file.close()
                    force_file.close()
        #GETS CELL AND STORES ######################################################################################        
            if line.find("Supercell") != -1 or CELL_FLAG == True:
                CELL_FLAG = True
                fields_cell = line.split()
                if fields_cell[0] != "Supercell" and fields_cell[0] != "AtomData:":
                    cell_tmp.append([float(fields_cell[0]),float(fields_cell[1]),float(fields_cell[2])])
        #GETS FORCES AND STORES
            if line.find("AtomData:") != -1 or FORCE_FLAG == True:
                CELL_FLAG = False
                FORCE_FLAG = True
                fields_forces=line.split()
                if fields_forces[0] != "AtomData:" and fields_forces[0] != "Energy":
                    forces.append([float(fields_forces[5])/ self.__Ry_to_eV__ * self.__Bohr_to_A__,float(fields_forces[6])/ self.__Ry_to_eV__ * self.__Bohr_to_A__,float(fields_forces[7])/ self.__Ry_to_eV__ * self.__Bohr_to_A__])
        #GETS ENERGY AND STORES ######################################################################################
            if line.find("Energy") != -1 or ENERGY_FLAG == True:
                FORCE_FLAG = False
                ENERGY_FLAG = True
                fields_energy=line.split()
                if fields_energy[0] != "Energy" and fields_energy[0] != "PlusStress:":
                    energy = float(fields_energy[0])/self.__Ry_to_eV__
        #GETS PRESSURES AND STORES ######################################################################################
            if line.find("PlusStress:") != -1 or PRESSURE_FLAG == True:
                ENERGY_FLAG = False
                PRESSURE_FLAG = True
                fields_pressure = line.split()
                if fields_pressure[0] != "PlusStress:" and fields_pressure[0] != "Feature":
                    pressure[0][0] = float(fields_pressure[0]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
                    pressure[1][1] = float(fields_pressure[1]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
                    pressure[2][2] = float(fields_pressure[2]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
                    pressure[1][2] = float(fields_pressure[3]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
                    pressure[2][1] = float(fields_pressure[3]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
                    pressure[0][2] = float(fields_pressure[4]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
                    pressure[2][0] = float(fields_pressure[4]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
                    pressure[0][1] = float(fields_pressure[5]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
                    pressure[1][0] = float(fields_pressure[5]) / self.__Ry_to_eV__ * self.__Bohr_to_A__ * self.__Bohr_to_A__ * self.__Bohr_to_A__
        #STOPS LOOP  ######################################################################################
            if line.find("Feature") != -1:
                PRESSURE_FLAG = False
            
            if line.find("Feature   conf_number") != -1:
                conf_number = int(line.split()[-1])  
                                                    
        energy_file.close()        
        conf_output_file.close()
        return TOTAL_NUMBER_OF_MLIP_CONF        
        


