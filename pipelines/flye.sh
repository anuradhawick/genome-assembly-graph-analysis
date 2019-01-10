[ -e $out_dir/flye_yeast ] && rm -r $out_dir/flye_yeast && mkdir $out_dir/flye_yeast
flye --pacbio-raw $out_dir/yeast_reads/all.fastq --out-dir $out_dir/flye_yeast --threads 32 --genome-size 12m
    