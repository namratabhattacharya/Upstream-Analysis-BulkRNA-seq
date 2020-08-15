# Upstream-Analysis


Input - GSE75688
Study - Single-cell RNA sequencing of primary breast cancer
Here I have considered four samples, two from tumor and 2 from non-tumor as follows-
BC01_53 = tumor, SRR2973289
BC01_50 = non_tumor, SRR5023384
BC01_72 = tumor, SRR2973295
BC01_74 = non_tumor, SRR2973296
Pipeline followed -  
Trim-galore, STAR, RSEM, R(EdgeR, Limma)

Run Command -
bash DEpipeline1.sh
Rscript DEpipeline.R
The parameters in the script are static. I will change it after the discussion.

Output- Here I have tried to do differential gene expression analysis between tumor and non-tumor samples with DE levels.
Output Files are attached in the mail -  Result_DE_tumor_v_nontumor.txt, Rplots.pdf
Voom plot is not giving the correct result. Further filtration is required.
