[ -e $out_dir/canu_yeast ] && rm -r $out_dir/canu_yeast && mkdir $out_dir/canu_yeast
canu -p r1 -pacbio-raw $out_dir/yeast_reads/all.fastq -d $out_dir/canu_yeast genomeSize=12m maxMemory=60 maxThreads=32 useGrid=false # java=$java_path
