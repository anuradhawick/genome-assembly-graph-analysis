module load python3/3.6.2

PATH=$PATH:/home/659/aw5153/canu/Linux-amd64/bin/:/home/659/aw5153/flye/bin/:/home/659/aw5153/PBSIM-PacBio-Simulator/src/
export model_clr_path=/home/659/aw5153/PBSIM-PacBio-Simulator/data/model_qc_clr
export genonepath=/home/659/aw5153/genome-assembly-graph-analysis/pipelines/input/yeast/YeastAll.fa

# export model_clr_path=/media/anuradhawick/data/Tools/PBSIM-PacBio-Simulator/data/model_qc_clr
# export genonepath=/media/anuradhawick/data/Experiments/Assembly_Graph/pipelines/input/yeast/YeastAll.fa



python3 assembly_improve.py