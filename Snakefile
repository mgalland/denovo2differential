##############################################################################################################
# Produce de novo assembly of RNA-Seq reads + differential analysis + annotation of the assembled transcripts#
##############################################################################################################
import os


FASTQDIR = "/home/mgalland/data/01_sequencing_runs/20151013_72523_InHolland/"
#SAMPLES = list(set([s[0:3] for s in os.listdir(FASTQDIR)])) # set function allows to filter duplicate samples


DIRS = ["unzip","./trim/"]

configfile:"config.json"

rule all:
    input: 
        expand("gunzip/{data}_{r}.fastq",data=config["data"],r=["R1","R2"])




########### Rules ###################
rule gunzip:
    input:
        forward = lambda wildcards: config["basedir"] + config["data"][wildcards.data]["forward"],
        reverse = lambda wildcards: config["basedir"] + config["data"][wildcards.data]["reverse"]
    output:
        forward = "gunzip/{data}_R1.fastq",
        reverse = "gunzip/{data}_R2.fastq"
    message:"unzipping original files"
    shell:"""
          gunzip -c {input.forward} > {output.forward}
          gunzip -c {input.reverse} > {output.reverse}
	  """

