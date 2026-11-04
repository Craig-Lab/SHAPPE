#!/bin/bash
#SBATCH --array=371-380,661-670,701-710
#SBATCH --time=50:00:00
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --output=logsPhysicellWellMixed_NoLocalRecruitment/%x_%A_%a.out
#SBATCH --error=logsPhysicellWellMixed_NoLocalRecruitment/%x_%A_%a.err

mkdir -p logsPhysicellWellMixed_NoLocalRecruitment

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

executable=GBMTAM
project=glioblastomaTargetingTAMs
alpha=10
init=0
mvt=1
use_voids=0
averageModel=0

module load intel
module load python/3.10
module load matlab/2024b.1
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirementsPhysiCellAllianceCan.txt
python3.10 simulations.py $executable $project $alpha $init $mvt $use_voids $averageModel $SLURM_ARRAY_TASK_ID