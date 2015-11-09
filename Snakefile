##############################################################################################################
# Produce de novo assembly of RNA-Seq reads + differential analysis + annotation of the assembled transcripts#
##############################################################################################################
import os

# configuration file for fastq file directory and other specific parameters
configfile:"config.json" 

# Trimmomatic (read trimming)
TRIMMOMATIC_PARAMS = config["trimmomatic_params"]
TRIMMOMATIC = "/zfs/datastore0/software/src/Trimmomatic-0.30/trimmomatic-0.30.jar"

# Trinity parameters
TRINITY_ASSEMBLY_PARAMS = config["trinity_assembly_params"]
TRINITY_ESTIMATE_ABUNDANCE = config["trinity_estimate_abundance"]
TRINITY_ABUNDANCE_PARAMS = config["trinity_abundance_params"]

rule all:
    input: 
        expand("trim/{data}_{r}.fastq.gz",data=config["data"],r=["FP","FU","RP","RU"]),
        "trinity_out_dir.Trinity.fasta",
        expand("express/{data}/results.xprs",data=config["data"])
    message:"all done"



########### Rules ###################



###################################################
# rule annotation of de novo assembled transcripts#
###################################################



####################################
# rule quantification with eXpress
###################################
rule transcript_abundance:
    input:
        FP = "trim/{data}_FP.fastq.gz",
        RP = "trim/{data}_RP.fastq.gz",
        assembly = "trinity_out_dir.Trinity.fasta"
    output:
        res = "express/{data}/results.xprs"        
    message:"estimating transcript abundance for {wildcards.data}"
    params:"express/{data}/"
    shell:
        "{TRINITY_ESTIMATE_ABUNDANCE} --transcripts {input.assembly} "
        "--seqType fq "
        "--left {input.FP} --right {input.RP} "
        "--output_dir {params} "
        "{TRINITY_ABUNDANCE_PARAMS}"
    

rule denovo:
    input:
        left = "trim/all.left.fastq.gz",
        right = "trim/all.right.fastq.gz"
    output:"trinity_out_dir.Trinity.fasta"
    message:"de novo assembly of all reads"
    log:"assembly.log.txt"
    shell:
        "Trinity --seqType fq --left {input.left} --right {input.right} "
        "{TRINITY_ASSEMBLY_PARAMS} 2>{log}"

rule gzip:
    input:
        left = "trim/all.left.fastq",
        right = "trim/all.right.fastq"
    output:
        left = "trim/all.left.fastq.gz",
        right = "trim/all.right.fastq.gz"
    message:"gunzip concatenated left/right read files"
    shell:"""
        gzip --force {input.left}
        gzip --force {input.right}
	  """

rule concatenate_reads:
    input:
        left = expand("trim/{data}_FP.fastq.gz",data=config["data"]),
        right = expand("trim/{data}_RP.fastq.gz",data=config["data"])
    output:
        left = "trim/all.left.fastq",
        right = "trim/all.right.fastq"
    message:"concatenating all reads in one file" 
    shell:"""
        touch {output.left}
        touch {output.right}
	gunzip -d -c {input.left}| cat - >> {output.left}
	gunzip -d -c {input.right}| cat - >> {output.right}
	"""

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
    log:"trim/{data}_log.txt"
    shell:
        "java -jar {TRIMMOMATIC} PE -phred33 {input.forward} {input.reverse} "
        "{output.FP} {output.FU} {output.RP} {output.RU} "
        "{TRIMMOMATIC_PARAMS} "
        "2>{log}"
