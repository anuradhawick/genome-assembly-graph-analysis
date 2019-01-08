from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio import SeqIO
import glob

from utils.command_liner import cmd_exit

def generate_reads():
    cmd_exit("[ -e output/Sampled_Reads_Yeast ] && rm -r output/Sampled_Reads_Yeast", "cleaning path")
    cmd_exit("mkdir 'output/Sampled_Reads_Yeast'", "creating the folder")




    file_list = list(glob.glob("./input/yeast/chr*.fa"))
    records = []

    for f in file_list:
        strings = open(f).read().split("\n")
        record = SeqRecord(Seq("".join(strings[1::]), DNAAlphabet), id=strings[0][1::], description="")
        records.append(record)

    SeqIO.write(records, "./input/yeast/YeastAll.fa", "fasta")

    command = """cd 'output/Sampled_Reads_Yeast' && pbsim --data-type CLR --depth 50 --model_qc $model_clr_path --length-min 20000 --length-max 30000 --difference-ratio 0:0:0 '/media/anuradhawick/data/Experiments/Assembly_Graph/pipelines/input/yeast/YeastAll.fa'"""

    cmd_exit(command, "generating reads from yeast")

    file_list = list(glob.glob("./output/Sampled_Reads_Yeast/*.fastq"))
    file_list = map(str, file_list)

    command = "cat " + " ".join(file_list) + " > ./output/Sampled_Reads_Yeast/all.fastq"

    cmd_exit(command, "gathering all reads into one file")