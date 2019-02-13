module load bcftools/1.9

echo INFO::Running phasing calls:preprocess BCF -> VCF
bcftools view -o $out_dir/vcf_files/initial_call.vcf $out_dir/vcf_files/initial_call.bcf
echo INFO::Running phasing calls:WhatsHap execution
whatshap phase --ignore-read-groups -r $out_dir/contigs.fasta -o $out_dir/vcf_files/phased.vcf $out_dir/vcf_files/initial_call.vcf $out_dir/mappings/mapping.sorted.bam
