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
( clusterProfiler, enrichplot, org.Hs.eg.db, ggnewscale, pathview, DOSE, plyr, R.utils)

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
|MyoT_starv2D_WT_PCi1426_B2 | 19 037 718 | 18 739 410 | 97,74 | 18 606 686 | 11 960 119 | 62,82 | 
|MyoT_starv2D_WT_PCi1426_B3 | 20 685 887 | 20 150 636 | 96,49 | 19 958 820 | 10 303 221 | 49,81 | 
|MyoT_starv2D_WT_PCi1426_B4 | 19 649 172 | 19 229 002 | 96,99 | 19 057 821 | 9 758 237 | 49,66 | 
|MyoT_starv2D_WT_i90_B1 | 13 839 354 | 13 594 602 | 97,57 | 13 503 131 | 7 570 791 | 54,70 | 
|MyoT_starv2D_WT_i90_B3 | 20 926 037 | 20 546 139 | 96,84 | 20 264 054 | 10 460 614 | 49,99 | 
|MyoT_starv2D_WT_i90_B4 | 18 124 459 | 17 672 598 | 96,68 | 17 523 180 | 8 049 908 | 44,41 | 
|MyoT_starv2D_WT_iCas9_B1 | 18 643 894 | 18 245 146 | 96,89 | 18 064 432 | 9 971 946 | 53,49 | 
|MyoT_starv2D_WT_iCas9_B2 | 19 568 711 | 18 975 160 | 95,47 | 18 682 222 | 10 373 533 | 53,01 | 
|MyoT_starv2D_WT_iCas9_B3 | 17 904 744 | 17 471 592 | 96,35 | 17 250 421 | 10 662 038 | 59,55 | 
|MyoT_starv2D_GSDIII_AGL_KO_B1 | 19 563 667 | 19 092 760 | 96,63 | 18 904 644 | 12 035 856 | 61,52 | 
|MyoT_starv2D_GSDIII_AGL_KO_B3 | 19 483 435 | 19 101 814 | 97,14 | 18 927 070 | 12 203 586 | 62,64 | 
|MyoT_starv2D_GSDIII_AGL_KO_B4 | 20 173 895 | 19 943 780 | 98,05 | 19 780 000 | 12 380 047 | 61,37 | 
|MyoT_starv2D_GSDIII_0576c12_B1 | 18 098 984 | 17 788 945 | 97,46 | 17 639 370 | 10 445 084 | 57,71 | 
|MyoT_starv2D_GSDIII_0576c12_B2 | 17 216 572 | 16 951 254 | 97,52 | 16 789 342 | 10 137 569 | 58,88 | 
|MyoT_starv2D_GSDIII_0576c12_B3 | 19 492 542 | 19 205 297 | 97,45 | 18 995 864 | 11 358 042 | 58,27 | 
|MyoT_starv2D_GSDIII_3390c15_B1 | 20 838 302 | 20 593 953 | 98,13 | 20 449 068 | 13 509 946 | 64,83 | 
|MyoT_starv2D_GSDIII_3390c15_B2 | 19 506 588 | 19 213 283 | 97,70 | 19 057 839 | 12 165 331 | 62,37 | 
|MyoT_starv2D_GSDIII_3390c15_B4 | 19 949 287 | 19 708 899 | 98,15 | 19 580 354 | 12 002 370 | 60,16 | 
|MyoT_starv2D_GSDIII_0303c3_B1 | 21 259 436 | 20 991 268 | 97,97 | 20 827 905 | 12 039 543 | 56,63 | 
|MyoT_starv2D_GSDIII_0303c3_B3 | 19 299 978 | 18 978 957 | 97,56 | 18 829 832 | 9 796 632 | 50,76 | 
|MyoT_starv2D_GSDIII_0303c3_B4 | 21 210 397 | 20 916 959 | 97,91 | 20 768 008 | 11 399 817 | 53,75 | 
|MyoT_starv2D_GSDIII_2523c3_B1 | 20 946 365 | 20 678 390 | 98,06 | 20 541 022 | 11 100 370 | 52,99 | 
|MyoT_starv2D_GSDIII_2523c3_B2 | 21 056 505 | 20 771 638 | 97,74 | 20 580 864 | 11 293 400 | 53,63 | 
|MyoT_starv2D_GSDIII_2523c3_B3 | 20 426 475 | 20 154 835 | 97,79 | 19 974 088 | 11 749 369 | 57,52 | 

