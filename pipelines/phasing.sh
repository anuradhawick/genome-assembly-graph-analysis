#!/bin/bash
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=128GB
#PBS -l jobfs=1GB
#PBS -l walltime=20:00:00
#PBS -l wd

export out_dir=/short/qr59/aw5153/pipeline_output_hap_5_vc
export graph_path=/short/qr59/aw5153/pipeline_output_hap_5/flye_yeast/assembly_graph.gfa
export reads_path=/short/qr59/aw5153/pipeline_output_hap_5/yeast_reads/all.fastq

echo INFO::Running phasing calls:preprocess BCF -> VCF
bcftools view -o $out_dir/vcf_files/initial_call.vcf $out_dir/vcf_files/initial_call.bcf
echo INFO::Running phasing calls:WhatsHap execution
whatshap phase -o $out_dir/vcf_files/phased.vcf $out_dir/vcf_files/initial_call.vcf $out_dir/mappings/mapping.sorted.bam
