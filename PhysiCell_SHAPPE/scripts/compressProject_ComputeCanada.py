import subprocess
import os
import shutil

path_project = "/Users/morgancraig/Documents/Documents - MacBook Air de wksadm/GBMs/ProjectABMODE"
path_makefiles = "/Users/morgancraig/Documents/Documents - MacBook Air de wksadm/GBMs/ProjectABMODE/Makefiles/ComputeCanada"
path_projectPhysicell = os.path.join(path_project,"PhysiCell_TargetingTAMs")

# Firt input: from | Second: to
shutil.copyfile(os.path.join(path_makefiles,"Makefile"),os.path.join(path_projectPhysicell,"Makefile"))
shutil.copyfile(os.path.join(path_makefiles,"Makefile-backup"),os.path.join(path_projectPhysicell,"Makefile-backup"))
shutil.copyfile(os.path.join(path_makefiles,"Makefile-default"),os.path.join(path_projectPhysicell,"sample_projects","Makefile-backup"))
shutil.copyfile(os.path.join(path_makefiles,"glioblastomaTargetingTAMsProject","Makefile"),os.path.join(path_projectPhysicell,"sample_projects","GlioblastomaTargetingTAMs","Makefile"))

os.chdir(path_project)
subprocess.run(["tar","-czf","PhysiCellTemplateTest.tar.gz","-C","PhysiCell_TargetingTAMs","."])