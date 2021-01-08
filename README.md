# Chernobyl Thyroid Cancer - miRNAseq
## I. Description
This workflow was used for general QC and alignment of microRNA-seq data in the Chernbobyl thyroid cancer study. The sequencing reads were processed according to the ENCODE microRNA-seq pipeline (https://www.encodeproject.org/microrna/microrna-seq/#references). 

Major steps in the workflow are:
1) Trimming of adapters and low-quality reads using cutadapt and only reads with lengths of 15-31 nt were kept
2) Generating QC reports using FASTQC and aggregating results using multiQC
3) Aligning trimmed reads to GRCh38 human reference genome (illumine iGenomes NCBI GRCh38) using STAR, and microRNA reads count was quantified according to GENCODE V24 genome annotation file which was the microRNA subset from comprehensive GENCODE annotations.
4) Merging reads-count tables of all samples
## II. Dependencies
* [Python](https://www.python.org)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Multiqc](https://multiqc.info)
* [Star](https://github.com/alexdobin/STAR)
* [R](https://www.r-project.org)
## III. Input requirements
* [config.yaml](https://github.com/NCI-CGR/ChernobylThyroidCancer-miRNAseq/blob/main/config.yaml)
* [sample_names.txt](https://github.com/NCI-CGR/ChernobylThyroidCancer-miRNAseq/blob/main/sample_names.txt)
* [merged fastq files stored in directory: merged_fastq/](https://github.com/NCI-CGR/ChernobylThyroidCancer-miRNAseq/tree/main/merged_fastq)
* [adapters.fa with adapter sequences](https://github.com/NCI-CGR/ChernobylThyroidCancer-miRNAseq/blob/main/adapters.fa)
* reference genome sequence and annotation files
* [microRNA annotation file](https://github.com/NCI-CGR/ChernobylThyroidCancer-miRNAseq/tree/main/star_index)
## IV. Working directory structure
```bash
.
├── adapters.fa
├── config.yaml
├── log
├── merged_fastq
│   └── {sample}.fastq.gz
├── merge.R
├── posttrim_qc
│   ├── posttrim_qc_multiqc_report.html
│   ├── {sample}.trim_fastqc.html
│   └── {sample}.trim_fastqc.zip
├── pretrim_qc
│   ├── pretrim_qc_multiqc_report.html
│   ├── {sample}_fastqc.html
│   └── {sample}_fastqc.zip
├── reads_count
│   └── reads_count.csv
├── run.sh
├── sample_names.txt
├── Snakefile
├── snakemake.batch
├── star_align
│   ├── log
│   │   └── star_align_multiqc_report.html
│   └── {sample}
│       ├── {sample}Aligned.sortedByCoord.out.bam
│       ├── {sample}Log.final.out
│       └── {sample}ReadsPerGene.out.tab
├── star_index
│   ├── chrLength.txt
│   ├── chrNameLength.txt
│   ├── chrName.txt
│   ├── chrStart.txt
│   ├── ENCFF628BVT.gtf
│   ├── exonGeTrInfo.tab
│   ├── exonInfo.tab
│   ├── gencode.v24.primary_assembly.annotation-tRNAs-ERCC_phiX.gtf
│   ├── geneInfo.tab
│   ├── Genome
│   ├── genomeParameters.txt
│   ├── SA
│   ├── SAindex
│   ├── sjdbInfo.txt
│   ├── sjdbList.fromGTF.out.tab
│   ├── sjdbList.out.tab
│   └── transcriptInfo.tab
└── trimmed
    ├── {sample}_too_long.fastq.gz
    ├── {sample}_too_short.fastq.gz
    ├── {sample}.trim.fastq.gz
    └── {sample}.trim.log
```
