### This pipeline is for general QC and alignment of miRNA-seq data

## vim: ft=python
import sys
import os
import glob
import itertools

shell.prefix("set -eo pipefail; ")
localrules: all

# define wildcards
def parse_sampleID(fname):
    return fname.split('/')[-1].split('_')[0]

file = sorted(glob.glob('merged_fastq/*.fastq.gz'), key=parse_sampleID)

d = {}
for key, value in itertools.groupby(file, parse_sampleID):
    d[key] = list(value)
       
rule all:
    input:
          expand("star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=d.keys()),
          "pretrim_qc/preQC_multiqc_report.html",
          "posttrim_qc/postQC_multiqc_report.html",
          "star_align/log/star_align_multiqc_report.html",
          "reads_count/reads_count.csv"         

rule cutadapt:
    input: 
          "merged_fastq/{sample}.fastq.gz"
    output:
          "trimmed/{sample}.trim.fastq.gz",
          "trimmed/{sample}_too_short.fastq.gz",
          "trimmed/{sample}_too_long.fastq.gz"
    params: 
          adapter = config["ADAPTER_FA"]
    threads: 10
    shell:
          """
          cutadapt -b {params.adapter} -m 15 -M 31 --too-short-output={output[1]} --too-long-output={output[2]} -q 10,10 -o {output[0]} {input} 2>log/{wildcards.sample}_cutadapt.err
          """

rule pretrim_qc:
    input: 
          "merged_fastq/{sample}.fastq.gz" 
    output: 
          "pretrim_qc/{sample}_fastqc.zip",
          "pretrim_qc/{sample}_fastqc.html"
    threads: 10
    shell:
          """
          fastqc {input} -o pretrim_qc -f fastq --noextract 2>log/{wildcards.sample}_preqc.err
          """

rule posttrim_qc:
    input: 
          "trimmed/{sample}.trim.fastq.gz"
    output: 
          "posttrim_qc/{sample}.trim_fastqc.zip",
          "posttrim_qc/{sample}.trim_fastqc.html"
    threads: 10
    shell:
          """
          fastqc {input} -o posttrim_qc -f fastq --noextract 2>log/{wildcards.sample}_postqc.err
          """

rule star_index:
    input:
          "/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa",
          "star_index/gencode.v24.primary_assembly.annotation-tRNAs-ERCC_phiX.gtf"
    output:
          "star_index/complete.txt"
    threads: 24
    shell:
          """
          STAR --runThreadN 24 --runMode genomeGenerate --genomeDir star_index --sjdbGTFfile {input[1]} --sjdbOverhang 1 --genomeFastaFiles {input[0]} 2>log/star_index.err 
          touch {output}
          """

rule star_align:
    input:
          "star_index/ENCFF628BVT.gtf",
          "trimmed/{sample}.trim.fastq.gz"
          "star_index/complete.txt"
    output:
          "star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",
          "star_align/{sample}/{sample}Log.final.out",
          "star_align/{sample}/{sample}ReadsPerGene.out.tab"
    threads: 24
    params:
          index="star_index"
    shell:
          """
          STAR --runThreadN 24 --genomeDir {params.index} --readFilesIn {input[1]} --outFileNamePrefix star_align/{wildcards.sample}/{wildcards.sample} --readFilesCommand zcat --sjdbGTFfile {input[0]} --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapScoreRange 0 --quantMode TranscriptomeSAM GeneCounts --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 --outSAMunmapped Within --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16 --alignSJDBoverhangMin 1000 --alignIntronMax 1 --outWigType wiggle --outWigStrand Stranded --outWigNorm RPM 2>log/{wildcards.sample}_star_align.err 
          """

rule multiqc:
    input:
          expand("pretrim_qc/{sample}_fastqc.html",sample=d.keys()),
          expand("posttrim_qc/{sample}.trim_fastqc.html",sample=d.keys()),
          expand("star_align/{sample}/{sample}Log.final.out",sample=d.keys())
    output:
          "pretrim_qc/preQC_multiqc_report.html",
          "posttrim_qc/postQC_multiqc_report.html",
          "star_align/log/star_align_multiqc_report.html"
    threads: 8
    shell:
          """
          multiqc pretrim_qc/. --title preQC -o pretrim_qc 2>log/multiqc_preqc.err
          multiqc posttrim_qc/. --title postQC -o posttrim_qc/ 2>log/multiqc_postqc.err
          mkdir star_align/log
          cp star_align/*/*Log.final.out star_align/log
          multiqc star_align/log/. --title star_align -o star_align/log 2>log/multiqc_star.err
          """

## library is sense stranded
rule merge:
    input:  
          expand("star_align/{sample}/{sample}ReadsPerGene.out.tab",sample=d.keys())
    output:
          "reads_count/reads_count.csv"
    shell:
          """
          Rscript merge.R 2>log/merge_count.err
          """

