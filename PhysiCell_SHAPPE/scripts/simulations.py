#AUTHOR: Blanche Mongeon

import subprocess
import numpy as np
import os
import sys
import shutil #this will be used to copy files by using the command shutil.copyfile(src, dst)
              #in that case we are copying the content of the file src into the file dst
import xml.etree.ElementTree as ET #to parse through and modify a xml file
from datetime import datetime
import pandas as pd

# Get the arguments that were passed 
arg = sys.argv
if len(arg)>1:
    executable = sys.argv[1] #this will be the project called for all simulations
    projectName = sys.argv[2]
    alphaCD8 = int(sys.argv[3])
    init = int(sys.argv[4])
    mvt = int(sys.argv[5])
    use_voids = int(sys.argv[6])
    averageModel = int(sys.argv[7])
    jobNumber = int(sys.argv[8])
else:
    sys.exit('You must enter the name of the executable project for all simulations.')

Nconfigurations = 10
Nreplicates =  3 # Number of replicates per configuration
# Retrieve the sample's number and the configuration
if averageModel==1:
    patientNumber = 0
    nC = jobNumber
else:
    patientNumber = int(np.ceil(float(jobNumber)/Nconfigurations))
    nC = int(jobNumber-(patientNumber-1)*Nconfigurations)

# Set folders information and path
main_folder = os.getcwd() # ProjectABMODE folder
IMCinfo_folder = os.path.join(main_folder,'IMC_Data')
### Read patients information - will be useful to name folders
dfPatients = pd.read_csv(os.path.join(IMCinfo_folder,'fullPhysicellDataframe.csv'))
px_nbr = str(dfPatients.csv_ID[patientNumber])


RepulsiveRegionsFolder = os.path.join(IMCinfo_folder,"RepulsiveRegions")

TargetDistancesFolder = os.path.join(IMCinfo_folder,"TargetDistances",px_nbr)

# init=0 -> CSR | init=1 -> ECSR | init=2 -> SHAPPE
if(init==0 and use_voids==0):
    configurations_folder = os.path.join(IMCinfo_folder,'CSR_Reconstructions_ABMtypes',px_nbr)
    # Get the path to final outputMatlab folder
    final_output_folder = os.path.join(main_folder,'finalOutputNoLocalRecruitment_WellMixed')
    outfile_folder = os.path.join(main_folder,'outfilesNoLocalRecruitment_WellMixed')
    errorFileName = "simulationsWellMixed_"+px_nbr+"_config"+str(nC)+".txt"
elif(init==1 and use_voids==0):
    configurations_folder = os.path.join(IMCinfo_folder,'ECSR_Reconstructions_ABMtypes',px_nbr)
    # Get the path to final outputMatlab folder
    final_output_folder = os.path.join(main_folder,'finalOutputNoLocalRecruitment_ECSR')
    outfile_folder = os.path.join(main_folder,'outfilesNoLocalRecruitment_ECSR')
    errorFileName = "simulationsECSR_"+px_nbr+"_config"+str(nC)+".txt"
elif(init==1 and use_voids==1):
    configurations_folder = os.path.join(IMCinfo_folder,'ECSR_Reconstructions_ABMtypes',px_nbr)
    # Get the path to final outputMatlab folder
    final_output_folder = os.path.join(main_folder,'finalOutputNoLocalRecruitment_RandomInitialisation')
    outfile_folder = os.path.join(main_folder,'outfilesNoLocalRecruitment_RandomInitialisation')
    errorFileName = "simulationsRandomInitialisation_"+px_nbr+"_config"+str(nC)+".txt"
elif(init==2):
    configurations_folder = os.path.join(IMCinfo_folder,'Reconstructions_ABMtypes',px_nbr)
    final_output_folder = os.path.join(main_folder,'finalOutputNoLocalRecruitment')
    outfile_folder = os.path.join(main_folder,'outfilesNoLocalRecruitment')
    errorFileName = "simulations_"+px_nbr+"_config"+str(nC)+".txt"
else:
    sys.exit('4th inout can only take the value 1 (for random initialisation) or 2 (for SHAPPE initialisation).')
os.makedirs(final_output_folder, exist_ok=True)
os.makedirs(outfile_folder, exist_ok=True)

# Probabilities of repulsive regions
voidProbabilities = np.loadtxt(os.path.join(IMCinfo_folder,'voidProbs.csv'), delimiter=',')
voidCarryingCapacity = np.loadtxt(os.path.join(IMCinfo_folder,'voidProbsAreaUnit.csv'), delimiter=',')
if(patientNumber>0):
    voidProb = voidProbabilities[patientNumber-1]
    voidCC = voidCarryingCapacity[patientNumber-1]
else:
    voidProb = 0
    voidCC = 0

# Copy physicell folder into a new physicell folder on temporary storage ($SLURM_TMPDIR)
temporary_folder = os.path.join(os.environ['SLURM_TMPDIR'],'Folder_'+px_nbr+'_Config'+str(nC))
temporary_physicell_folder = os.path.join(temporary_folder,'PhysiCell'+px_nbr+'_Config'+str(nC))
# Create the directory if it does not already exist
os.makedirs(temporary_physicell_folder, exist_ok=True)
# Copy archived PhysiCell folder
template_physicell_folder = os.path.join(main_folder,'PhysiCellTemplate.tar.gz')
subprocess.run(["rsync", "--inplace", "--no-whole-file", template_physicell_folder, temporary_folder], check=True)
compressed_file = os.path.join(temporary_folder,'PhysiCellTemplate.tar.gz')
# Decompress file
subprocess.run(["tar", "-xzf", compressed_file, "-C", temporary_physicell_folder], check=True)
# Remove compressed file to save storage space
subprocess.run(["rm", compressed_file], check=True)

# Move to your temporary physicell folder
os.chdir(temporary_physicell_folder)

code_folder = os.getcwd()
output_dir_long = os.path.join(code_folder,'output')

errorFilePath = os.path.join(code_folder,errorFileName)
matlab_folder = os.path.join(code_folder,'matlab')

# Treatments we want to simulate
treatments = ["None"]
NTreatments = len(treatments) # Number of simulated treatments per patient
TMZsimul = ["false"]

xml_extension = ".xml"
template_file = "PhysiCell_settings"

errorFile = open(errorFilePath,"a")
now = datetime.now()
date_string = now.strftime("%d/%m/%Y %H:%M")
errorFile.write("\n"+date_string+": Starting "+px_nbr)

# Loop through your simulated treatments
### nT is the treatment index
for nT in range(NTreatments):

    now = datetime.now()
    date_string = now.strftime("%d/%m/%Y %H:%M")
    errorFile.write("\n"+date_string+": now running treatment "+treatments[nT]+" for "+px_nbr)

    # Make the project
    subprocess.run(["make","clean"],stdout=errorFile)
    errorFile.write('\n')
    subprocess.run(["make","reset"],stdout=errorFile)
    errorFile.write('\n')
    subprocess.run(["make",projectName],stdout=errorFile)
    errorFile.write('\n')
    subprocess.run(["make"],stdout=errorFile)
    errorFile.write("\n Project created.")

        
    # Keep track of what simulation you are starting and when
    now = datetime.now()
    date_string = now.strftime("%d/%m/%Y %H:%M")
    errorFile.write("\n"+date_string+": Starting the replicates for configurations "+
                    str(nC)+" out of "+str(Nconfigurations)+".")
    
    # Copy the initial conditions
    shutil.copyfile(os.path.join(configurations_folder,"rep"+str(nC)+".csv"),os.path.join(code_folder,'cells.csv'))
    shutil.copyfile(os.path.join(configurations_folder,"rep"+str(nC)+".csv"),os.path.join(code_folder,'config','cells.csv'))

    # Loop through your replicate
    for nR in range(Nreplicates):
        file_name = px_nbr+"_Treatment"+treatments[nT]+"_Config"+str(nC)+"_Replicate"+str(nR+1)

        # Copy the template file into the new project file
        shutil.copyfile(os.path.join('config',template_file+xml_extension),os.path.join('config',file_name+xml_extension))

        # Parse through and modify your xml file as needed
        tree = ET.parse(os.path.join('config',file_name+xml_extension))
        root = tree.getroot()
        rootTag = root.tag

        userparam_element = root.find('user_parameters')
        ### TMZ Treatment Boolean
        TMZ_tx = userparam_element.find("TMZ_tx")
        TMZ_tx.text = TMZsimul[nT]

        ### Alpha value for the recruitment of CD8+ T cells
        alphaCD8_XML = userparam_element.find("cd8_recruitment_alpha")
        alphaCD8_XML.text = str(float(alphaCD8/10.0))

        ### Tumour radius
        tumour_radius_initial = userparam_element.find("tumour_radius_initial")
        tumour_radius_initial.text = str(dfPatients.synthetic_radius_circle[patientNumber])
        tumour_radius = userparam_element.find("tumour_radius")
        tumour_radius.text = str(dfPatients.synthetic_radius_circle[patientNumber])

        ### PhysiCell Domain
        domainsize_elements = root.find('domain')
        x_min = domainsize_elements.find("x_min")
        x_min.text = str(-2*dfPatients.synthetic_radius_circle[patientNumber])
        x_max = domainsize_elements.find("x_max")
        x_max.text = str(2*dfPatients.synthetic_radius_circle[patientNumber])

        y_min = domainsize_elements.find("y_min")
        y_min.text = str(-2*dfPatients.synthetic_radius_circle[patientNumber])
        y_max = domainsize_elements.find("y_max")
        y_max.text = str(2*dfPatients.synthetic_radius_circle[patientNumber])

        ### URL for voids
        url_voids = userparam_element.find("url_voids")
        url_voids.text = os.path.join(RepulsiveRegionsFolder,px_nbr,"rep"+str(nC)+".csv")
        
        ### Void probability
        void_CCXML = userparam_element.find("void_carrying_capacity")
        void_CCXML.text = str(voidCC)
        errorFile.write("voidCC: "+str(voidCC))
        void_probXML = userparam_element.find("void_probability")
        void_probXML.text = str(voidProb)
        errorFile.write("voidProb: "+str(voidProb))

        ### INITIAL CONDITIONS
        ########## G(0)
        G0 = userparam_element.find("G0")
        N_G0 = dfPatients.N_G[patientNumber]
        G0.text = str(N_G0)
        ########## T(0)
        T0 = userparam_element.find("T_0")
        N_T0 = dfPatients.N_CD8[patientNumber]
        T0.text = str(N_T0)
        ########## S(0)
        S0 = userparam_element.find("S0")
        N_S0 = dfPatients.N_S[patientNumber]
        S0.text = str(N_S0)  
        ########## M0(0)
        M0_0 = userparam_element.find("M0_0")
        N_M00 = dfPatients.N_M0[patientNumber]
        M0_0.text = str(N_M00)
        ########## M1(0)
        M1_0 = userparam_element.find("M1_0")
        N_M10 = dfPatients.N_M1[patientNumber]
        M1_0.text = str(N_M10)
        ########## M2(0)
        M2_0 = userparam_element.find("M2_0")
        N_M20 = dfPatients.N_M2[patientNumber]
        M2_0.text = str(N_M20)

        ####### TME BIAS TOWARDS M2 PHENOTYPE
        TAM_recruitement_IC50 = userparam_element.find("TAM_recruitment_IC50")
        TAM_recruitement_IC50.text = str((1/99)*N_G0)
        TAM_recruitement_M0 = userparam_element.find("TAM_recruitment_M0")
        TAM_recruitement_M1 = userparam_element.find("TAM_recruitment_M1")
        TAM_recruitement_M2 = userparam_element.find("TAM_recruitment_M2")
        if((N_M00+N_M10+N_M20)>0):
            TAM_recruitement_M0.text = str(N_M00/(N_M00+N_M10+N_M20))
            TAM_recruitement_M1.text = str(N_M10/(N_M00+N_M10+N_M20))
            TAM_recruitement_M2.text = str(N_M20/(N_M00+N_M10+N_M20))
        else:
            TAM_recruitement_M0.text = str(0.0)
            TAM_recruitement_M1.text = str(0.0)
            TAM_recruitement_M2.text = str(0.0)
        M2proba = userparam_element.find("M2_probability")
        if((N_M10+N_M20)>0):
            M2proba.text = str(N_M20/(N_M10+N_M20))
        else:
            M2proba.text = str(0.0)

        ### Set random initialisation to false
        use_repulsive_regions_XML = userparam_element.find("use_repulsive_regions")
        if use_voids==0:
            use_repulsive_regions_XML.text = str(0.0)
        else:
            use_repulsive_regions_XML.text = str(1.0)

        ### Oriented movement
        oriented_movement_XML = userparam_element.find("Oriented_Movement")
        if mvt==1:
            oriented_movement_XML.text = "false"
        else:
            oriented_movement_XML.text = "true"

        ### Target distances
        target_distances_XML = userparam_element.find("target_distances_folder")
        target_distances_XML.text = TargetDistancesFolder

        # Save the resulting xml file
        tree.write(os.path.join('config',file_name+xml_extension))

        # Run the project
        outfile_name = 'outfile_'+file_name+'.txt'
        outfile_path = os.path.join(code_folder,outfile_name)
        with open(outfile_name, "w") as outfile:
            result = subprocess.run(["./"+executable,"config/"+file_name+xml_extension],stdout=outfile, stderr=outfile)

        if result.returncode != 0:
            errorFile.write(f"\n{file_name}: GBMTAM exited {result.returncode}, skipping MATLAB")
            errorFile.flush()
            continue

        # Process simulations on matlab
        if(init==1 and use_voids==0):
            px_nbr_matlab = px_nbr+'_random'
        elif(init==1):
            px_nbr_matlab = px_nbr+'_randomInitialisation'
        else:
            px_nbr_matlab = px_nbr
        matlab_cmd = (
                f"try, "
                f"processSimulation('{px_nbr_matlab}','{treatments[nT]}',{((nC-1)*Nreplicates)+nR+1},'{output_dir_long}','{final_output_folder}'); "
                f"catch ME, disp(getReport(ME,'extended')); exit(1); "
                f"end; exit(0);"
            )

        with open(outfile_path, "a") as outfile:
            try:
                subprocess.run(["matlab","-batch",matlab_cmd], check=True, cwd=matlab_folder,stdout=outfile, stderr=outfile)
            except subprocess.CalledProcessError as e:
                errorFile.write(f"\nMATLAB failed for {file_name}: {e}")
                errorFile.flush()

        subprocess.run(["make","data-cleanup"],stdout=errorFile)

        # Copy the outfile
        shutil.copyfile(outfile_path,os.path.join(outfile_folder,outfile_name))
        errorFile.flush()
        
errorFile.write("\n "+px_nbr+" completed!")
errorFile.close()

# Copy the error file
shutil.copyfile(errorFilePath,os.path.join(outfile_folder,errorFileName))