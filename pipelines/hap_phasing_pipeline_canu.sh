#!/bin/bash
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=128GB
#PBS -l jobfs=1GB
#PBS -l walltime=10:00:00
#PBS -l wd

export out_dir=/short/qr59/aw5153/pipeline_output_hap_5_vc_canu
export contigs_path=/short/qr59/aw5153/pipeline_output_hap_5/canu_yeast/r1.contigs.fasta
export reads_path=/short/qr59/aw5153/pipeline_output_hap_5/yeast_reads/all.fastq

module use -a /short/qr59/aw5153/modules/

module load samtools/1.9
module load minimap2/2.15
module load bcftools/1.9

[ -e $out_dir ] && rm -r $out_dir
mkdir -p $out_dir/mappings

echo INFO::Running mapping functions

minimap2 --cs -a $contigs_path $reads_path > $out_dir/mappings/mapping.sam
samtools view -Sb $out_dir/mappings/mapping.sam > $out_dir/mappings/mapping.bam
samtools sort $out_dir/mappings/mapping.bam > $out_dir/mappings/mapping.sorted.bam
samtools index $out_dir/mappings/mapping.sorted.bam

mkdir -p $out_dir/vcf_files

sh ./hap_phasing_pipeline/pileup.sh
sh ./hap_phasing_pipeline/call.sh
sh ./hap_phasing_pipeline/phase.sh