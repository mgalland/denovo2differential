##############################################################################################################
# Produce de novo assembly of RNA-Seq reads + differential analysis + annotation of the assembled transcripts#
##############################################################################################################
import os
import io
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
from denovotools import stdout_redirector

# configuration file for fastq file directory and other specific parameters
configfile:"config.json" 

# Trimmomatic (read trimming)
TRIMMOMATIC = config["trimmomatic"]["jar"]
TRIMMOMATIC_PARAMS = config["trimmomatic"]["params"]

# Trinity & related tools parameters 
TRINITY_ASSEMBLY_PARAMS = config["trinity"]["assembly_params"]
TRINITY_ESTIMATE_ABUNDANCE = config["trinity"]["estimate_abundance"]
TRINITY_ABUNDANCE_PARAMS = config["trinity"]["abundance_params"]
TRINITY_UTILS = config["trinity"]["utils"]
TRANSDECODER = config["transdecoder"]["dir"]
TRANSDECODER_MIN_PROT_LENGTH = config["transdecoder"]["min_prot_length"]
TRANSDECODER_MIN_ORF_NUCLEOTIDE_LENGTH = config["transdecoder"]["min_orf_nucleotide_length"]

# BLAST
BLAST_HEADER = config["blast"]["blastn"]["params"]["outfmt"]
BLASTN_PARAMS = list(config["blast"]["blastn"]["params"].values())
BLASTP_PARAMS = list(config["blast"]["blastp"]["params"].values())

# THREADS
THREADS = 8

####### Outputs ######
QUALS = expand("plots/{data}.quality_boxplot.png",data=config["454"]["data"])
ILLUMINA_TRIMMED_READS = expand("trim/{illumina}_{r}.fastq.gz",illumina=config["illumina"],r=["FP","FU","RP","RU"])
ASSEMBLY = "trinity/trinity_out_dir.Trinity.fasta"
QC_MAPPING = "qc/mapping.txt"
SAMS = expand("qc/{illumina}.mapping.sam",illumina=config["illumina"])
TRINITY_XPRS = expand("transcript_abundance/{illumina}/results.xprs",illumina=config["illumina"])
TRINITY_PEP = "trinity/transcripts.fasta.transdecoder.pep"
BLASTN = "blast/trinity_vs_nt.outfmt6"
BLASTP = ["blast/trinitry_vs_swissprot.outfmt6","blast/trinity_vs_Rgenes.outfmt6","blast/trinity_vs_carotenoids.outfmt6","blast/trinity_vs_anthocyanins.outfmt6"]
VERSIONS = "report/software_versions.txt"

rule all:
    input: 
        ILLUMINA_TRIMMED_READS,
	QUALS,
        ASSEMBLY,
	QC_MAPPING,
	SAMS,
        TRINITY_XPRS,
        BLASTN,
	BLASTP,
        VERSIONS
    message:"all done"


###################################################
# Report
###################################################
rule software_versions:
    output:
        "report/software_versions.txt"
    message:"compiling tool versions"
    run:
        shell("mkdir -p report/")
        trimmomatic_version = os.path.basename(TRIMMOMATIC)
        shell("Trinity --version |cat >> {output}")
        shell("bowtie2 --version |cat >> {output}")
        shell("blastn -h |cat >> {output}")
        shell("express --help |cat >> {output}")
        with open(output[0],"w") as fileout:
            fileout.write(trimmomatic_version + "\n")


###################################################
## Filter blastp output
########################

###################################################
# rule annotation of de novo assembled transcripts#
###################################################
rule blastp_against_anthocyanins:
    input:"trinity/trinity_out_dir.Trinity.shortnames.fasta.transdecoder.pep"
    output:"blast/trinity_vs_anthocyanins.outfmt6"
    message:"blasting against anthocyanins genes"
    shell:
        "blastp -db /home/mgalland/data/02_refs/tulip/anthocyanins.fasta "
        "{BLASTP_PARAMS} "
        "-query {input} -out {output} "
        "-num_threads {THREADS}"


rule blastp_against_carotenoids:
    input:"trinity/trinity_out_dir.Trinity.shortnames.fasta.transdecoder.pep"
    output:"blast/trinity_vs_carotenoids.outfmt6"
    message:"blasting against carotenoid genes"
    shell:
        "blastp -db /home/mgalland/data/02_refs/tulip/carotenoids.fasta "
        "{BLASTP_PARAMS} "
        "-query {input} -out {output} "
        "-num_threads {THREADS}"

rule blastp_against_Rgenes:
    input:"trinity/trinity_out_dir.Trinity.shortnames.fasta.transdecoder.pep"
    output:"blast/trinity_vs_Rgenes.outfmt6"
    message:"blasting against R genes"
    shell:
        "blastp -db /home/mgalland/data/02_refs/tulip/R_genes.fasta "
        "{BLASTP_PARAMS} "
        "-query {input} -out {output} "
        "-num_threads {THREADS}"

rule blastp_against_swissprot:
    input:"trinity/trinity_out_dir.Trinity.shortnames.fasta.transdecoder.pep"
    output:"blast/trinitry_vs_swissprot.outfmt6"
    message:"blastp {input} against swissprot"
    params:
        db = config["blast"]["dbdir"] + config["blast"]["db"]["swissprot"]   
    shell:
        "blastp -db {params.db} "
        "{BLASTP_PARAMS} "
        "-query {input} -out {output} "
        "-num_threads {THREADS}"
 

rule trinity_transdecoder:
    input:"trinity/trinity_out_dir.Trinity.shortnames.fasta"
    output:"trinity/trinity_out_dir.Trinity.shortnames.fasta.transdecoder.pep"
    message:"predicting ORF within transcripts"
    shell:
        "{TRANSDECODER} -t {input} {TRANSDECODER_MIN_PROT_LENGTH} "
        "{TRANSDECODER_MIN_ORF_NUCLEOTIDE_LENGTH} "
        "--CPU {THREADS};" 
        "mv trinity_out_dir* trinity/;"
        "rm -r transdecoder.tmp*"

rule blastn:
    input:
        "trinity/trinity_out_dir.Trinity.shortnames.fasta"
    output:
        "blast/trinity_vs_nt.outfmt6"
    message:"blastn of the Trinity assembled transcripts"
    params: 
        db = config["blast"]["dbdir"] + config["blast"]["db"]["nt"]   
    shell:
        "blastn -db {params.db} "
        "{BLASTN_PARAMS} "
        "-query {input} -out {output} "
        "-num_threads {THREADS}"
    
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
rule abundance_estimates_to_matrix:
    input: "transcript_abundance/{illumina}/results.xprs" 
    output:"transcript_abundance/matrix.TMM.EXPR.matrix"
    message:"TMM normalization of eXpress values"
    params:"transcript_abundance/*/results.xprs"
    run:
        shell(TRINITY_UTILS + "abundance_estimates_to_matrix.pl --est_method eXpress --name_sample_by_basedir {params}")

rule estimate_transcript_abundance:
    input:
        FP = "trim/{illumina}_FP.fastq.gz",
        RP = "trim/{illumina}_RP.fastq.gz",
        assembly = "trinity/trinity_out_dir.Trinity.fasta"
    output:
        res = "transcript_abundance/{illumina}/results.xprs"        
    message:"estimating transcript abundance for {wildcards.illumina}"
    params:"express/{illumina}/"
    shell:
        "{TRINITY_ESTIMATE_ABUNDANCE} --transcripts {input.assembly} "
        "--seqType fq "
        "--left {input.FP} --right {input.RP} "
        "--output_dir {params} "
        "{TRINITY_ABUNDANCE_PARAMS}"

###################################
# QC of the assembly
###################################
rule map_all2contigs:
    input:
        left = "trim/all.left.fastq.gz",
        right = "trim/all.right.fastq.gz",
        assembly = "trinity/trinity_out_dir.Trinity.fasta"
    output:
        log = "qc/mapping.txt",
        sam = "qc/mapping.sam"
    message:"mapping reads back to Trinity assembly to estimate quality of assembly"
    shell:
        "bowtie2-build {input.assembly} trinity/assembly ;"
        "bowtie2 -x trinity/assembly -1 {input.left} -2 {input.right} -S {output.sam} 2>{output.log}"

rule map_sample2contigs:
    input:
        left = expand("trim/{illumina}_FP.fastq.gz",illumina=config["illumina"]),
        right = expand("trim/{illumina}_RP.fastq.gz",illumina=config["illumina"]),
        assembly = "trinity/trinity_out_dir.Trinity.fasta"
    output:
       sam = "qc/{illumina}.mapping.sam"
    message:"mapping sample reads back to Trinity assembly"
    log: "qc/{illumina}.mapping_log.txt"
    shell:
        "bowtie2-build {input.assembly} trinity/assembly ;"
        "bowtie2 -x trinity/assembly -1 {input.left} -2 {input.right} -S {output.sam} 2>{log}"

#####################################    
# de novo assembly of the reads
######################################
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
        left = expand("trim/{illumina}_FP.fastq.gz",illumina=config["illumina"]),
        right = expand("trim/{illumina}_RP.fastq.gz",illumina=config["illumina"])
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
# Trimming and QC of 454 reads
#####################################
rule fastx_boxplot:
    input:
        "trim/{data}.trimmed.stats.txt"
    output:
        "plots/{data}.quality_boxplot.png"
    message:"generating boxplot for quality scores for {wildcards.data}"
    shell:"fastq_quality_boxplot_graph.sh -i {input} -o {output} -t {wildcards.data}"

rule fastx_quality_stats:
    input:
        "trim/{data}.trimmed.fastq"
    output:
        "trim/{data}.trimmed.stats.txt"
    message:"generating fastx quality stats for {wildcards.data}"
    shell:"fastx_quality_stats -i {input} -o {output}"
    
rule fastx_quality_filter:
    input: 
        reads = lambda wildcards: config["454"]["dir"] + config["454"]["data"][wildcards.data]      
    output: "trim/{data}.trimmed.fastq"
    message:"trimming {input.reads} file"
    log:"fastq/{wildcards.data}.log.txt"
    shell:"fastq_quality_filter -v -q 10 -p 80 -i {input.reads} -o {output} 2>{log}"

##############################################
# Trimming and QC of Illumina paired-end reads
##############################################
rule trimmomatic_for_IlluminaPE:
    input:
        forward = lambda wildcards: config["basedir"] + config["illumina"][wildcards.illumina]["forward"],
        reverse = lambda wildcards: config["basedir"] + config["illumina"][wildcards.illumina]["reverse"]
    output:
        FP = "trim/{illumina}_FP.fastq.gz", # forward paired
        RP = "trim/{illumina}_RP.fastq.gz", # reverse paired
        FU = "trim/{illumina}_FU.fastq.gz", # forward unpaired
        RU = "trim/{illumina}_RU.fastq.gz"  # reverse unpaired
    message:"Trimming {wildcards.illumina} file"
    log:"trim/{illumina}_log.txt"
    shell:
        "java -jar {TRIMMOMATIC} PE -phred33 {input.forward} {input.reverse} "
        "{output.FP} {output.FU} {output.RP} {output.RU} "
        "{TRIMMOMATIC_PARAMS} "
        "2>{log}"

