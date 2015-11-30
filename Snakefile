##############################################################################################################
# Produce de novo assembly of RNA-Seq reads + differential analysis + annotation of the assembled transcripts#
##############################################################################################################
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# configuration file for fastq file directory and other specific parameters
configfile:"config.json" 

# input sequencing files
DIR454 = config["454"]["dir"]
DICT_454 = config["454"]["data"]
FILES_454 = [os.path.join(DIR454,f) for f in list(DICT_454.values())]
NAMES_454 = list(DICT_454.keys())

# Trimmomatic (read trimming)
TRIMMOMATIC = config["trimmomatic"]["jar"]
TRIMMOMATIC_PARAMS = config["trimmomatic"]["params"]

# Trinity & related tools parameters 
TRINITY_ASSEMBLY_PARAMS = config["trinity"]["assembly_params"]
TRINITY_ESTIMATE_ABUNDANCE = config["trinity"]["estimate_abundance"]
TRINITY_ABUNDANCE_PARAMS = config["trinity"]["abundance_params"]
TRANSDECODER = config["transdecoder"]["dir"]
TRANSDECODER_MIN_PROT_LENGTH = config["transdecoder"]["min_prot_length"]
TRANSDECODER_MIN_ORF_NUCLEOTIDE_LENGTH = config["transdecoder"]["min_orf_nucleotide_length"]

# BLAST
BLAST_HEADER = config["blastn"]["outfmt"]

####### Outputs ######
QUALS = expand("trim/{data}.png",data=config["454"]["data"])
ILLUMINA_TRIMMED_READS = expand("trim/{illumina}_{r}.fastq.gz",illumina=config["illumina"],r=["FP","FU","RP","RU"])
ASSEMBLY = "trinity/trinity_out_dir.Trinity.fasta"
QC_MAPPING = "qc/mapping.txt"
XPRS = expand("transcript_abundance/{illumina}/results.xprs",illumina=config["illumina"])
ASSEMBLY_SHORTNAMES_NR = "trinity/trinity_out_dir.Trinity.shortnames.nr.fasta"
BLASTN = "blastn/trinity.outfmt6"
VERSIONS = "report/software_versions.txt"

rule all:
    input: 
        ILLUMINA_TRIMMED_READS,
	QUALS,
        ASSEMBLY,
	QC_MAPPING,
        XPRS,
        ASSEMBLY_SHORTNAMES_NR,
        BLASTN,
        VERSIONS
    message:"all done"


###################################################
# Report
###################################################
rule software_verions:
    output:
        "report/software_versions.txt"
    message:"compiling tool versions"
    run:
        shell("mkdir -p report/")
        shell("Trinity --version > report/software_versions.txt")
        trimmomatic_version = os.path.basename(TRIMMOMATIC)
        shell("cat {trimmomatic_version} >> report/software_versions.txt")
        shell("bowtie2 --version |head -1 >> report/software_versions.txt")
        shell("blastn -h |grep -A1 'DESCRIPTION' >> report/software_versions.txt")
        shell("express --help 2>test.txt; cat test.txt |head -2 >> report/software_versions.txt;rm test.txt")
        

###################################################
# rule annotation of de novo assembled transcripts#
###################################################


# rule blastp search + PFAM search

rule trinity_transdecoder:
    input:"trinity/trinity_out_dir.Trinity.fasta"
    output:"transcripts.fasta.transdecoder.pep"
    message:"predicting ORF within transcripts"
    shell:
        "{TRANSDECODER} -t {input} -m" 
        "mv transcripts.fasta.transdecoder.pep trinity/transcripts.fasta.transdecoder.pep"

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
rule map2contigs:
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
rule fastx_boxplot:
    input:
        "trim/{data}.trimmed.stats.txt"
    output:
        "trim/{data}.png"
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

#########################################
rule trimmomatic_for_Illumina:
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
