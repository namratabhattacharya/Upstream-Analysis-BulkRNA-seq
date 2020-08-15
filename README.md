# Upstream-Analysis

## Pipeline followed -  
Differential Expression Analysis using Trim-galore, STAR, RSEM, R(EdgeR, Limma)

## Run Command -


```bash
bash DEpipeline1_v2.sh
```
DEpipelinev1.sh will produce the following result-
1. Directory */Result_ReadQuant* will store RSEM output ( genes.results, isoforms.results)
2. Directory */Result_AlignPct* will store Trim-galore output 
3. *filelabel.txt* stored in */Result_ReadQuant* is a tab-seperated file which contain sample information. 

```bash
Rscript DEpipeline2_v2.R
```
