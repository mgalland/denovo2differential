##############################################################################################################
# Produce de novo assembly of RNA-Seq reads + differential analysis + annotation of the assembled transcripts#
##############################################################################################################
import os

# configuration file for fastq file directory and other specific parameters
configfile:"config_test.json" 

# TOOLS
TRIMMOMATIC = "/zfs/datastore0/software/src/Trimmomatic-0.30/trimmomatic-0.30.jar"


# Parameters taken from config file
TRIMMOMATIC_PARAMS = config["trimmomatic_params"]
TRINITY_PARAMS = config["trinity_params"]

DIRS = ["unzip","./trim/"]



rule all:
    input: 
        expand("trim/{data}_{r}.fastq.gz",data=config["data"],r=["FP","FU","RP","RU"]),
        #expand("trinity/{data}.Trinity.fasta",data=config["data"]),
        "trinity_out_dir.Trinity.fasta"


########### Rules ###################


####################################
# rule annotation of de novo assembled transcripts


####################################
# rule quantification with eXpress

####################################

# trinity_out_dir.Trinity.fasta


rule denovo:
    input:
        left = "trim/all.left.fastq.gz",
        right = "trim/all.right.fastq.gz"
    output:"trinity_out_dir.Trinity.fasta"
    message:"de novo assembly of all reads"
    log:"assembly.log.txt"
    shell:
        "Trinity --seqType fq --left {input.left} --right {input.right} "
        "{TRINITY_PARAMS} 2>{log}"


rule zip_left:
    input:
        left = "trim/all.left.fastq"
    output:
        left = "trim/all.left.fastq.gz"
    message:"gunzip concatenated left read files"
    shell:"""
        gzip --force {input.left}
	  """

rule zip_right:
    input:
         right = "trim/all.right.fastq"
    output:
         right = "trim/all.right.fastq.gz"
    message:"gunzip concatenated right read files"
    shell:"""
        gzip --force {input.right}
	  """
rule concatenate_right_reads:
    input:expand("trim/{data}_RP.fastq.gz",data=config["data"])
    output:"trim/all.right.fastq"
    message:"concatenating all right reads in one file" 
    shell:"""
        touch {output}
	gunzip -d -c {input}| cat - >> {output}
	"""

rule concatenate_left_reads:
    input:expand("trim/{data}_FP.fastq.gz",data=config["data"])
    output:"trim/all.left.fastq"
    message:"concatenating all left reads in one file" 
    shell:"""
        touch {output}
	gunzip -d -c {input}| cat - >> {output}
	"""

####################################
#rule single_denovo:
 #   input:
  #      FP = "trim/{data}_FP.fastq.gz",
   #     RP = "trim/{data}_RP.fastq.gz"
  #  output:"trinity/{data}/Trinity.fasta"
  #  message:"de novo assembly of {wildcards.data} transcripts using Trinity - single assembly"
  #  log:"trinity/{data}_trinity_log.txt"
  #  shell:
   #     "mkdir -p ./assembly/{wildcards.data}/
    #    "Trinity --seqType fq --left {input.FP} --right {input.RP} "
     #   "{TRINITY_PARAMS} 2>{log}|"
      #  "mv trinity_out_dir.Trinity.fasta /trinity/{ 

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
