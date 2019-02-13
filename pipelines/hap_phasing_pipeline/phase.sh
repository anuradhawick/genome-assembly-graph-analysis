echo INFO::Running phasing calls:preprocess BCF -> VCF
bcftools view -o $out_dir/vcf_files/initial_call.vcf $out_dir/vcf_files/initial_call.bcf
echo INFO::Running phasing calls:WhatsHap execution
whatshap phase --ignore-read-groups -r $contigs_path -o $out_dir/vcf_files/phased.vcf $out_dir/vcf_files/initial_call.vcf $out_dir/mappings/mapping.sorted.bam
