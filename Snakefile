##############################################################################################################
# Produce de novo assembly of RNA-Seq reads + differential analysis + annotation of the assembled transcripts#
##############################################################################################################
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# configuration file for fastq file directory and other specific parameters
configfile:"config.json" 

# Trimmomatic (read trimming)
TRIMMOMATIC_PARAMS = config["trimmomatic_params"]
TRIMMOMATIC = "/zfs/datastore0/software/src/Trimmomatic-0.30/trimmomatic-0.30.jar"

# Trinity parameters
TRINITY_ASSEMBLY_PARAMS = config["trinity_assembly_params"]
TRINITY_ESTIMATE_ABUNDANCE = config["trinity_estimate_abundance"]
TRINITY_ABUNDANCE_PARAMS = config["trinity_abundance_params"]

####### Outputs ######
TRIMMED_READS = expand("trim/{data}_{r}.fastq.gz",data=config["data"],r=["FP","FU","RP","RU"])
ASSEMBLY = "trinity/trinity_out_dir.Trinity.fasta"
EXPRESS_XPRS = expand("express/{data}/results.xprs",data=config["data"])
ASSEMBLY_SHORTNAMES_NR = "trinity/trinity_out_dir.Trinity.shortnames.nr.fasta"
BLASTN = "blastn/trinity.outfmt6"


rule all:
    input: 
        TRIMMED_READS,
        ASSEMBLY,
        EXPRESS_XPRS,
        ASSEMBLY_SHORTNAMES_NR,
        BLASTN
    message:"all done"


###################################################
# rule annotation of de novo assembled transcripts#
###################################################
rule blastn:
    input:
        "trinity/trinity_out_dir.Trinity.shortnames.fasta"
    output:
        "blastn/trinity.outfmt6"
    message:"blastn of the Trinity assembled transcripts"
    params:
        db = lambda wildcards: config["blastn"]["db"],
        targets = lambda wildcards: config["blastn"]["max_target_seqs"],
        threads = lambda wildcards: config["blastn"]["num_threads"],
        outfmt = lambda wildcards: config["blastn"]["outfmt"],
        hsps = lambda wildcards: config["blastn"]["max_hsps"]
    shell:
        "blastn -db {params.db} -query {input} -out {output} "
        "-max_target_seqs {params.targets} "
        "-num_threads {params.threads} "
        "-outfmt {params.outfmt} "
        "-max_hsps {params.hsps}"

#ule update_blastdb:
  #  input:
   #     lambda wildcards: config["blastdbdir"]
   # output:
    #message:
    #params:lambda wildcards: config["blast"]["dir"]
    #shell:
    #    "gunzip /home/mgall"
    #    "cat *.tar | tar -xvf - -i "
        

rule shorten_seq_names:
    input:
        "trinity/trinity_out_dir.Trinity.fasta"
    output:
        "trinity/trinity_out_dir.Trinity.shortnames.fasta"
    message:"shorten sequence names in Trinity assembly fasta file"
    run:
        with open(input[0],"r") as filin, open(output[0],"w") as fileout:
            records = []
            for record in SeqIO.parse(filin,"fasta"):
                short_name = record.id.split(" ")[0]
                records.append(SeqRecord(record.seq,id=short_name))
            SeqIO.write(records,fileout,"fasta")        

####################################
# rule quantification with eXpress
###################################
rule transcript_abundance:
    input:
        FP = "trim/{data}_FP.fastq.gz",
        RP = "trim/{data}_RP.fastq.gz",
        assembly = "trinity/trinity_out_dir.Trinity.fasta"
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
    output:protected("trinity/trinity_out_dir.Trinity.fasta")
    message:"de novo assembly of all reads"
    log:"assembly.log.txt"
    shell:
        "Trinity --seqType fq --left {input.left} --right {input.right} "
        "{TRINITY_ASSEMBLY_PARAMS} 2>{log};"
        "mv trinity_out* trinity/"

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
