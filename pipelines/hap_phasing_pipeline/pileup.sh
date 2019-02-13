module load bcftools/1.9

echo INFO::Running varian calls:pile up

bcftools mpileup --threads 32 -o $out_dir/vcf_files/study.bcf -Ou -f $contigs_path $out_dir/mappings/mapping.sorted.bam