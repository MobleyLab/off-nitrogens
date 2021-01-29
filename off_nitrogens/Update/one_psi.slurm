#!/bin/bash

#SBATCH --job-name="perturb_valence_60"
#SBATCH --partition=mf_nes2.8
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10gb
# one hour
#SBATCH --time=2:00:00
#SBATCH --distribution=block:cyclic
#--------------

# Informational output
echo "=================================== SLURM JOB ==================================="
echo
echo "The job will be started on the following node(s):"
echo $SLURM_JOB_NODELIST
echo
echo "Slurm User:         $SLURM_JOB_USER"
echo "Run Directory:      $(pwd)"
echo "Job ID:             $SLURM_JOB_ID"
echo "Job Name:           $SLURM_JOB_NAME"
echo "Partition:          $SLURM_JOB_PARTITION"
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"
echo "Number of tasks:    $SLURM_NTASKS"
echo "Submitted From:     $SLURM_SUBMIT_HOST"
echo "Submit directory:   $SLURM_SUBMIT_DIR"
echo "=================================== SLURM JOB ==================================="
echo


cd $SLURM_SUBMIT_DIR
echo 'Working Directory:'
pwd

date

/beegfs/DATA/mobley/limvt/local/miniconda3/envs/oepython3/bin/psi4 -i input.dat -o output.dat

date
