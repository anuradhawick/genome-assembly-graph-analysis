echo INFO::Running varian calls:variant discovery
bcftools call -c --threads 32 -O b -o $out_dir/vcf_files/initial_call.bcf  $out_dir/vcf_files/study.bcf