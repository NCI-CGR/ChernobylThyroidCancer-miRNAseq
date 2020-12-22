# Chernobyl Thyroid Cancer - miRNAseq
## I. Description
This workflow was used for general QC and alignment of microRNA-seq data in the Chernbobyl thyroid cancer study. The sequenced reads were processed according to the ENCODE microRNA-seq pipeline (https://www.encodeproject.org/microrna/microrna-seq/#references). 
Major steps in the workflow are:
1) Trimming of adapters and low-quality reads using cutadapt and only reads with lengths of 15-31 nt were kept
2) Generating QC reports using FASTQC and aggregating results using multiQC
3) Aligning trimmed reads to GRCh38 human reference genome (illumine iGenomes NCBI GRCh38) using STAR, and microRNA reads count was quantified according to GENCODE V24 genome annotation file which was the microRNA subset from comprehensive GENCODE annotations.
4) Merging reads-count tables of all samples
