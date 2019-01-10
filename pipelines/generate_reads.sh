export model_clr_path=/home/659/aw5153/PBSIM-PacBio-Simulator/data/model_qc_clr
export genomepath=/short/qr59/aw5153/genome-assembly-graph-analysis/pipelines/input/yeast/all.fa

# export model_clr_path=/media/anuradhawick/data/Tools/PBSIM-PacBio-Simulator/data/model_qc_clr
# export genomepath=/media/anuradhawick/data/Experiments/Assembly_Graph/pipelines/input/yeast/all.fa

# gather all reads to one file
cat input/yeast/chr*.fa > input/yeast/all.fa

# generate reads
mkdir -p $out_dir/yeast_reads
cd $out_dir/yeast_reads
pbsim --data-type CLR --depth 50 --model_qc $model_clr_path --length-min 20000 --length-max 30000 --difference-ratio 0:0:0 $genomepath
cat sd*.fastq > all.fastq
cd $working_path