##############################################################################################################
# Produce de novo assembly of RNA-Seq reads + differential analysis + annotation of the assembled transcripts#
##############################################################################################################
import os

# configuration file for fastq file directory and other specific parameters
configfile:"config.json" 

# TOOLS
TRIMMOMATIC = "/zfs/datastore0/software/src/Trimmomatic-0.30/trimmomatic-0.30.jar"


# Parameters taken from config file
TRIMMOMATIC_PARAMS = config["trimmomatic_params"]
TRINITY_PARAMS = config["trinity_params"]

DIRS = ["unzip","./trim/"]



rule all:
    input: 
        expand("trim/{data}_{r}.fastq.gz",data=config["data"],r=["FP","FU","RP","RU"]),
        expand("trinity/{data}.Trinity.fasta",data=config["data"])


########### Rules ###################



####################################
rule denovo:
    input:
        FP = "trim/{data}_FP.fastq.gz",
        RP = "trim/{data}_RP.fastq.gz"
    output:"trinity/{data}.Trinity.fasta"
    message:"de novo assembly of {wildcards.data} transcripts using Trinity"
    shell:
        "Trinity --seqType fq --left {input.FP} --right {input.RP}"
        "{TRINITY_PARAMS} --output {output}"

#####################################
rule trimmomatic:
    input:
        forward = lambda wildcards: config["basedir"] + config["data"][wildcards.data]["forward"],
        reverse = lambda wildcards: config["basedir"] + config["data"][wildcards.data]["reverse"]
    output:
        FP = "trim/{data}_FP.fastq.gz", # forward paired
        RP = "trim/{data}_RP.fastq.gz", # reverse paired
        FU = "trim/{data}_FU.fastq.gz", # forward unpaired
        RU = "trim/{data}_RU.fastq.gz"  # reverse unpaired
    message:"Trimming {wildcards.data} file"
    shell:"""java -jar {TRIMMOMATIC} PE {input.forward} {input.reverse} {output.FP} {output.FU} {output.RP} {output.RU} {TRIMMOMATIC_PARAMS}"""
