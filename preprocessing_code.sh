#!/bin/env bash
# AI usage note: this analysis code was written by the author. Claude Code
# (Anthropic) was used afterward, retrospectively, to organize the
# repository (file layout, README) — not to write or modify the analysis
# logic below. See AI_DISCLOSURE.md for the full disclosure.
## WHOLE-TISSUE SEQ DATA PRE-PROCESSING USING STAR-ALIGNER
module use ~/ModuleFiles/
module load STAR2.7

#Genome Generate
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ~/STAR_index --genomeFastaFiles ~/Drosophila_melanogaster.BDGP6.28.dna.toplevel_GAL4_EGFP.fa --sjdbGTFfile ~/Drosophila_melanogaster.BDGP6.28.100.chr_filtered.GAL4.EGFP.gtf --sjdbOverhang 49

#Align 24h-tjTS>GFP-repl1 (Control)
STAR --runThreadN 40 --outFileNamePrefix WT_24_1 --genomeDir ~/STAR_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --chimSegmentMin 10 --readFilesCommand gunzip -c --readFilesIn ~/A-1-1_ATCACG_L001_R1_001.fastq.gz,~/A-1-1_ATCACG_L001_R1_002.fastq.gz,~/A-1-1_ATCACG_L002_R1_001.fastq.gz,~/A-1-1_ATCACG_L002_R1_002.fastq.gz
#Align 96h-tjTS>GFP-repl1 (Control)
STAR --runThreadN 40 --outFileNamePrefix WT_96_1 --genomeDir ~/STAR_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --chimSegmentMin 10 --readFilesCommand gunzip -c --readFilesIn ~/A-1-2_CGATGT_L001_R1_001.fastq.gz,~/A-1-2_CGATGT_L001_R1_002.fastq.gz,~/A-1-2_CGATGT_L002_R1_001.fastq.gz,~/A-1-2_CGATGT_L002_R1_002.fastq.gz
#Align 24h-tjTS>GFP-repl2 (Control)
STAR --runThreadN 40 --outFileNamePrefix WT_24_2 --genomeDir ~/STAR_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --chimSegmentMin 10 --readFilesCommand gunzip -c --readFilesIn ~/A-2-1_TTAGGC_L001_R1_001.fastq.gz,~/A-2-1_TTAGGC_L001_R1_002.fastq.gz,~/A-2-1_TTAGGC_L002_R1_001.fastq.gz,~/A-2-1_TTAGGC_L002_R1_002.fastq.gz
#Align 96h-tjTS>GFP-repl2 (Control)
STAR --runThreadN 40 --outFileNamePrefix WT_96_2 --genomeDir ~/STAR_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --chimSegmentMin 10 --readFilesCommand gunzip -c --readFilesIn ~/A-2-2_TGACCA_L001_R1_001.fastq.gz,~/A-2-2_TGACCA_L001_R1_002.fastq.gz,~/A-2-2_TGACCA_L002_R1_001.fastq.gz,~/A-2-2_TGACCA_L002_R1_002.fastq.gz
#Align 24h-tjTS>lglRNAi-repl1 (Sample)
STAR --runThreadN 40 --outFileNamePrefix LglIR_24_1 --genomeDir ~/STAR_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --chimSegmentMin 10 --readFilesCommand gunzip -c --readFilesIn ~/B-1-1_ACAGTG_L001_R1_001.fastq.gz,~/B-1-1_ACAGTG_L001_R1_002.fastq.gz,~/B-1-1_ACAGTG_L002_R1_001.fastq.gz,~/B-1-1_ACAGTG_L002_R1_002.fastq.gz
#Align 96h-tjTS>lglRNAi-repl1 (Sample)
STAR --runThreadN 40 --outFileNamePrefix LglIR_96_1 --genomeDir ~/STAR_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --chimSegmentMin 10 --readFilesCommand gunzip -c --readFilesIn ~/B-1-2_GCCAAT_L001_R1_001.fastq.gz,~/B-1-2_GCCAAT_L001_R1_002.fastq.gz,~/B-1-2_GCCAAT_L002_R1_001.fastq.gz,~/B-1-2_GCCAAT_L002_R1_002.fastq.gz
#Align 24h-tjTS>lglRNAi-repl2 (Sample)
STAR --runThreadN 40 --outFileNamePrefix LglIR_24_2 --genomeDir ~/STAR_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --chimSegmentMin 10 --readFilesCommand gunzip -c --readFilesIn ~/B-2-1_CAGATC_L001_R1_001.fastq.gz,~/B-2-1_CAGATC_L001_R1_002.fastq.gz,~/B-2-1_CAGATC_L002_R1_001.fastq.gz,~/B-2-1_CAGATC_L002_R1_002.fastq.gz
#Align 96h-tjTS>lglRNAi-repl2 (Sample)
STAR --runThreadN 40 --outFileNamePrefix LglIR_96_2 --genomeDir ~/STAR_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --chimSegmentMin 10 --readFilesCommand gunzip -c --readFilesIn ~/B-2-2_ACTTGA_L001_R1_001.fastq.gz,~/B-2-2_ACTTGA_L001_R1_002.fastq.gz,~/B-2-2_ACTTGA_L002_R1_001.fastq.gz,~/B-2-2_ACTTGA_L002_R1_002.fastq.gz

#featureCounts:
featureCounts -f -t gene -g gene_name -a ~/Drosophila_melanogaster.BDGP6.28.100.chr_filtered.GAL4.EGFP.gtf -O -o featurecounts_bulkSeq.txt ~/*.bam

## SINGLE-CELL SEQ DATA PRE-PROCESSING
module load cellranger/3.0.1 Velocyto Subread-2.0.0

#Cellranger:
#w1118
cellranger count --transcriptome=~/refgenome4cellranger/ --id=Ctrl_FC --chemistry=SC3Pv2 --fastqs=~/fastq --sample=SI-GA-B4_1,SI-GA-B4_2,SI-GA-B4_3,SI-GA-B4_4 --indices=ACTTCATA,GAGATGAC,TGCCGTGG,CTAGACCT #SRA dtataset: SRX7814226
#cellranger count --transcriptome=~/refgenome4cellranger/ --id=Ctrl_FC --fastqs=~/fastq/ --sample=B_1,B_2,B_3,B_4 --chemistry=SC3Pv2 --indices=CAGTACTG,AGTAGTCT,GCAGTAGA,TTCCCGAC #Replicate1 #SRA dataset: SRP321151
#lglRNAi
cellranger count --transcriptome=~/refgenome4cellranger/ --id=LglKD_FC --fastqs=~/fastq/ --sample=B4_1,B4_2,B4_3,B4_4 --chemistry=SC3Pv2 --indices=GTCCGGTC,AAGATCAT,CCTGAAGG,TGATCTCA

#Velocyto (for each sample):
velocyto run -b ~/filtered_feature_bc_matrix/barcodes.tsv.gz -o ~/velocyto/ -m ~/dm6_rmsk.gtf ~/possorted_genome_bam.bam ~/Drosophila_melanogaster.BDGP6.28.100.chr_filtered.GAL4.EGFP.gtf
