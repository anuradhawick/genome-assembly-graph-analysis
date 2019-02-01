#!/bin/bash
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=128GB
#PBS -l jobfs=1GB
#PBS -l walltime=20:00:00
#PBS -l wd

export out_dir=/short/qr59/aw5153/pipeline_output_hap_5_vc

echo INFO::Running varian calls:variant discovery
bcftools call -c --threads 32 -O b -o $out_dir/vcf_files/initial_call.bcf  $out_dir/vcf_files/study.bcf