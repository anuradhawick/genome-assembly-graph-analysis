from read_generation.generate_reads import generate_reads 
from assembly.canu import run_canu
from assembly.flye import run_flye

generate_reads()
run_canu()
run_flye()
