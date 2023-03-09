# NGS82

## Description

- Team : [Pharmacologie des dystrophies musculaires](https://www.istem.eu/randd/pharmacologie-des-dystrophies-musculaires/) 
- Sequencing type : Quantseq
- Samples : 24
- Organism : Human GRCh37.87

## Objectives  

- Raw & normalized counts, 
- List of filtered DE genes (p-adj < 0.05 & |log2FC| >= 0.4 & BM >= 20). 

## Tools 

- FastQC ( v0.11.9 ) 
- cutadapt ( 4.1 ) 
- PRINSEQ-lite ( 0.20.4 ) 
- STAR ( 2.7.6a ) 
- samtools ( 1.13 )  
- HTSeq-count ( 1.99.2 ) 
- R Packages ( 4.1.2 ) 
( RColorBrewer, gplots, DESeq2, pheatmap , tidyverse, gridExtra, reshape, ggforce, cowplot, reshape2, viridis)

## DEG results 


|  | Down  | Up | Total |
| :------------------------- | :----:  | :----: | :----: |
| GSDIII vs WT | 3 | 6 | 9 |
| 0576c12 vs WT | 101 | 70 | 171 |
| 3390c15 vs WT | 155 | 132 | 287 |
| 2523c3 vs WT | 83 | 27 | 110 |
| 0303c3 vs WT | 90 | 73 | 163 |
| AGL_KO vs PCi1426_WT | 3 | 1 | 4 |



## Details samples 

- Size : 76 bp  
- Encoding : Sanger / Illumina 1.9 
- Date : 01/03/2022 

| Names | #Reads  | #Reads post-filters | STAR : % mapped reads  | # mapped reads | # Filtered mapped reads | STAR : % Filtered mapped reads |
| :----: | :--------:  | :--------: | :----:  | :--------: | :--------: | :----: |
| 1_S1 | xx | xx | xx | xx | xx | xx |

