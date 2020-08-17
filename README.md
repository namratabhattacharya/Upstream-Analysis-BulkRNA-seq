# Upstream-Analysis

## Pipeline followed   
Differential Expression Analysis using Trim-galore, STAR, RSEM, R(EdgeR, Limma)


## Run Command 
Run-time input to DEpipeline1_v2.sh are as follows- 
1. Genome type (human/mouse)
2. Number of samples
4. Read types (PE/SE)
5. Accession number e.g. SRR11771595
6. Label (subtype1/subtype2)

```bash
bash DEpipeline1_v2.sh
```
DEpipeline1_v2.sh will produce the following result-
1. Directory */Result_ReadQuant* will store RSEM output ( genes.results, isoforms.results)
2. Directory */Result_AlignPct* will store Trim-galore output 
3. *filelabel.txt* stored in */Result_ReadQuant* is a tab-seperated file which contain sample information

```bash
setwd( "/path/to/Result_ReadQuant" )
```
Before executing DEpipeline2_v2.R change the working directory to  */Result_ReadQuant*. 
```bash
Rscript DEpipeline2_v2.R
```
DEpipeline2_v2.R will produce the following result-
1. *DE_tumor_v_nontumor.txt* gives the list of differentially expressed genes
2.  Differential expression levels result will be shown in the terminal
3.  Rplots.pdf gives the MDS and voom plot
