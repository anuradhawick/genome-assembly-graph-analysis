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

awk '/^S/{print ">"$2"\n"$3}' $graph_path | fold > $out_dir/contigs.fasta

[ -e $out_dir ] && rm -r $out_dir
mkdir -p $out_dir/mappings

echo INFO::Running mapping functions

minimap2 --cs -a $out_dir/contigs.fasta $reads_path > $out_dir/mappings/mapping.sam
samtools view -Sb $out_dir/mappings/mapping.sam > $out_dir/mappings/mapping.bam
samtools sort $out_dir/mappings/mapping.bam > $out_dir/mappings/mapping.sorted.bam
samtools index $out_dir/mappings/mapping.sorted.bam

mkdir -p $out_dir/vcf_files

echo INFO::Running varian calls:pile up

bcftools mpileup --threads 32 -o $out_dir/vcf_files/study.bcf -Ou -f $contigs_path $out_dir/mappings/mapping.sorted.bam

echo INFO::Running varian calls:variant discovery

bcftools call --threads 32 -O b -o $out_dir/vcf_files/initial_call.bcf  $out_dir/vcf_files/study.bcf