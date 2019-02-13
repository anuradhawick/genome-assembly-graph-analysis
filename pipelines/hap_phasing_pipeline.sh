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

module load samtools/1.9

[ -e $out_dir ] && rm -r $out_dir
mkdir -p $out_dir/mappings

awk '/^S/{print ">"$2"\n"$3}' $graph_path | fold > $out_dir/contigs.fasta

echo INFO::Running mapping functions

minimap2 --cs -a $out_dir/contigs.fasta $reads_path > $out_dir/mappings/mapping.sam
samtools view -Sb $out_dir/mappings/mapping.sam > $out_dir/mappings/mapping.bam
samtools sort $out_dir/mappings/mapping.bam > $out_dir/mappings/mapping.sorted.bam
samtools index $out_dir/mappings/mapping.sorted.bam

mkdir -p $out_dir/vcf_files

sh ./hap_phasing_pipeline/pileup.sh
sh ./hap_phasing_pipeline/call.sh
sh ./hap_phasing_pipeline/phase.sh