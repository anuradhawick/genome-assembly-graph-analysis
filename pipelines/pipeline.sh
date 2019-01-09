#!/bin/bash
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l jobfs=20GB
#PBS -l walltime=10:00:00
#PBS -l wd

module load python3/3.6.2
module load python/2.7.11
module load java/jdk1.8.0_60

PATH=$PATH:/home/659/aw5153/canu/Linux-amd64/bin/:/home/659/aw5153/flye/bin/:/home/659/aw5153/PBSIM-PacBio-Simulator/src/
export model_clr_path=/home/659/aw5153/PBSIM-PacBio-Simulator/data/model_qc_clr
export genonepath=/home/659/aw5153/genome-assembly-graph-analysis/pipelines/input/yeast/YeastAll.fa

# export model_clr_path=/media/anuradhawick/data/Tools/PBSIM-PacBio-Simulator/data/model_qc_clr
# export genonepath=/media/anuradhawick/data/Experiments/Assembly_Graph/pipelines/input/yeast/YeastAll.fa



python3 assembly_improve.py > pipeline_run_$PBS_JOBID.log