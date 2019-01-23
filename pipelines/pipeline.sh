#!/bin/bash
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l jobfs=1GB
#PBS -l walltime=10:00:00
#PBS -l wd

module load python3/3.6.2
module load python/2.7.11
module load java/jdk1.8.0_60

PATH=$PATH:/home/659/aw5153/canu/Linux-amd64/bin/:/home/659/aw5153/flye/bin/:/home/659/aw5153/PBSIM-PacBio-Simulator/src/
java_path=$(which java)

export out_dir=/short/qr59/aw5153/pipeline_output_hap_5
# export out_dir=/media/anuradhawick/data/Experiments/Assembly_Graph/pipelines/output

mkdir -p $out_dir

export working_path=$(pwd)


sh generate_reads.sh > $out_dir/log1_$PBS_JOBID.log
sh canu.sh > $out_dir/log2_$PBS_JOBID.log
sh flye.sh > $out_dir/log3_$PBS_JOBID.log