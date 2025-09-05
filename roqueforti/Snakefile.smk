shell.executable("/bin/bash")
from os.path import join

# read config info into this namespace
#configfile: "config.yaml"
PATH_TO_DATA="/home/thibault/Penicillium_roqueforti/data"
RESULTS_PATH="/home/thibault/Penicillium_roqueforti/results" # Path and name of your results folder
REFERENCE_BOWTIE2="/home/thibault/Penicillium_roqueforti/reference_genomes/GCA_001599855.1_JCM_22842_assembly_v001_genomic.fna" # Path and file
DICT_SAMTOOLS="/home/thibault/Penicillium_roqueforti/reference_genomes/GCA_001599855.1_JCM_22842_assembly_v001_genomic.dict" ### TO DO 
ADAPTERS_FILE="/home/thibault/Penicilium_roqueforti/data/adapter_long_list.fa"

(STRAINS,SAMPLES,STRANDS) = glob_wildcards(join(PATH_TO_DATA,"{strain}/{sample}_{strand}.fastq.gz"))
#STRANDS=["R1","R2"]

rule all:
	input:
		expand("{REFERENCE_BOWTIE2}.1.bt2",REFERENCE_BOWTIE2=REFERENCE_BOWTIE2),
		expand("{DICT_SAMTOOLS}",DICT_SAMTOOLS=DICT_SAMTOOLS),
		#expand(join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_cleaned_{strand}.fastq.gz"),zip,strain=STRAINS,sample=SAMPLES,strand=STRANDS),
		#expand(join(RESULTS_PATH,"mapping/{strain}/{sample}.bam"),zip,strain=STRAINS,sample=SAMPLES),
		expand("{REFERENCE_BOWTIE2}.fai",REFERENCE_BOWTIE2=REFERENCE_BOWTIE2),
		#expand(join(RESULTS_PATH,"fastqc_cleaned/{strain}/{sample}"),zip,strain=STRAINS,sample=SAMPLES),
		expand(join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score_pos_markdup.bam"),zip,strain=STRAINS,sample=SAMPLES),
		expand(join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score_pos_markdup.bam.bai"),zip,strain=STRAINS,sample=SAMPLES),
		expand(join(RESULTS_PATH,"variant_calling/{strain}/{sample}.vcf"),zip,strain=STRAINS,sample=SAMPLES)
# indexation bowtie2
 
rule bowtie2_indexation:
	input:  
 		REFERENCE_BOWTIE2
	output:
 		"{REFERENCE_BOWTIE2}.1.bt2",
 		"{REFERENCE_BOWTIE2}.2.bt2",
 		"{REFERENCE_BOWTIE2}.3.bt2",
 		"{REFERENCE_BOWTIE2}.4.bt2",
 		"{REFERENCE_BOWTIE2}.rev.1.bt2",
 		"{REFERENCE_BOWTIE2}.rev.2.bt2"
	priority:20
	log:    
		"logs/bowtie2/index.err"
	shell:
		"bowtie2-build {input} {input} > {log} 2>&1"

# Clean reads with Trimmomatic:
rule trimmomatic:
	input:
		fv=join(PATH_TO_DATA,"{strain}/{sample}_R1.fastq.gz"),
		rv=join(PATH_TO_DATA,"{strain}/{sample}_R2.fastq.gz")

	output:
		out1=join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_cleaned_R1.fastq.gz"),
		out2=join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_R1_unpaired.fastq.gz"),
		out3=join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_cleaned_R2.fastq.gz"),
		out4=join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_R2_unpaired.fastq.gz")
	log:
		"logs/trimmomatic/{sample}_trimmomatic.err"
	threads: 4
	shell:
		"TrimmomaticPE -threads {threads} -phred33 {input.fv} {input.rv} {output.out1} {output.out2} {output.out3} {output.out4} ILLUMINACLIP:{ADAPTERS_FILE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 > {log} 2>&1"

# Mapping on reference with Bowtie2:

rule Bowtie2:
	input:
		fv=join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_cleaned_R1.fastq.gz"),
		rv=join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_cleaned_R2.fastq.gz"),
		ind=REFERENCE_BOWTIE2
	output:
		temp(join(RESULTS_PATH,"mapping/{strain}/{sample}.bam"))
	threads:4
	log:
		"logs/bowtie2/{sample}_bowtie2.err"
	shell:
		"bowtie2 --threads {threads} --very-sensitive-local --local --phred33 --rg-id {wildcards.strain} --rg SM:{wildcards.strain} --rg LB:{wildcards.strain} --rg PL:ILLUMINA -x {input.ind} -1 {input.fv} -2 {input.rv} -X 1000 | samtools view -bS - > {output}"


# fastqc after cleaning reads:

rule fastqc_after_cleaning:
	input:
		fv=join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_cleaned_R1.fastq.gz"),
		rv=join(RESULTS_PATH,"cleaned_reads/{strain}/{sample}_cleaned_R2.fastq.gz")
	output:
		join(RESULTS_PATH,"fastqc_cleaned/{strain}/{sample}")
	shell:
		"mkdir -p {output} ; fastqc {input.fv} {input.rv} -o {output}"

# indexing reference with samtools:

rule index_samtools:
	input:
		{REFERENCE_BOWTIE2}
	output:
		"{REFERENCE_BOWTIE2}.fai"
	priority:20
	shell:
		"samtools faidx {input}"

# making dict for reference with samtools:

rule dict_samtools:
	input:
		{REFERENCE_BOWTIE2}
	output:
		{DICT_SAMTOOLS}
	priority:20
	shell:
		"samtools dict -o {output} {input}"

# Sorting bam with samtools:

rule sort_samtools:
	input:
		join(RESULTS_PATH,"mapping/{strain}/{sample}.bam")
	output:
		temp(join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted.bam"))
	priority:15
	shell:
		"samtools sort -n -o {output} {input}"

# Adding Mate Score (MS) and Mate Cigar (MC) for markdup

rule add_score_samtools:
	input:
		join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted.bam")
	output:
		temp(join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score.bam"))
	priority:10
	shell:	
		"samtools fixmate -m {input} {output}"

# Position order with samtools:

rule position_samtools:
	input:
		join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score.bam")
	output:
		temp(join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score_pos.bam"))
	priority:5
	shell:
		"samtools sort -o {output} {input}"

# Marking dupe with samtools:

rule markdup_samtools:
	input:
		join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score_pos.bam")
	output:
		join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score_pos_markdup.bam")
	priority:3
	shell:
		"samtools markdup {input} {output}"

# Indexing bam with samtools:

rule index_bam_samtools:
	input:
		join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score_pos_markdup.bam")
	output:
		join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score_pos_markdup.bam.bai")
	priority:1
	shell:
		"samtools index {input}"

# Call germline SNPs and indels via local re-assembly of haplotypes (doc GATK):

rule haplotypecaller_gatk:
	input:
		entri1=join(RESULTS_PATH,"mapping/{strain}/{sample}_sorted_score_pos_markdup.bam"),
		ref=REFERENCE_BOWTIE2
	output:
		join(RESULTS_PATH,"variant_calling/{strain}/{sample}.vcf")
	shell:
		"gatk HaplotypeCaller -I {input.entri1} -O {output} -R {input.ref} -ERC GVCF --sample-ploidy 1"
 

