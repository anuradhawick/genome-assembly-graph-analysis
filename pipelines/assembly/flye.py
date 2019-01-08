from utils.command_liner import cmd_exit

def run_flye():
    command1 = "[ -e 'output/Flye_yeast' ] && rm -r 'output/Flye_yeast' && mkdir 'output/Flye_yeast'"
    command2 = "flye --pacbio-raw './output/Sampled_Reads_Yeast/all.fastq' --out-dir 'output/Flye_yeast' --threads 16 --genome-size 12m"
    
    cmd_exit(command1, "creating assembly output path for Flye")
    cmd_exit(command2, "running Flye")

    return