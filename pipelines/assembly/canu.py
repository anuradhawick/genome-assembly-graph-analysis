from utils.command_liner import cmd_exit

def run_canu():
    command1 = "[ -e 'output/Canu_yeast' ] && rm -r 'output/Canu_yeast' && mkdir 'output/Canu_yeast'"
    command2 = """canu -p r1 -pacbio-raw './output/Sampled_Reads_Yeast/all.fastq' -d 'output/Canu_yeast' genomeSize=12m cnsMemory=32 cnsConcurrency=32 gridEngineThreadsOption='-l ncpus=32'"""
    
    cmd_exit(command1, "creating assembly output path for canu")
    cmd_exit(command2, "running canu")

    return