library(limma)
library(edgeR)

#Read Files and label
input <- read.csv(file="filelabel.txt",head=TRUE,sep="\t")
file <- c(input$filename)
#temp = list.files(pattern="*.genes.results")
labels <- c(input$label)

#Create DGElist
dlist <- readDGE(file, columns=c(1,5))

#Display DGElist
cat("\nDisplay DGElist:\n")
dlist

#remove .genes from samplename
samplenames <- substring(colnames(dlist), 1, nchar(colnames(dlist))-6)
samplenames
colnames(dlist) <- samplenames

#Organise sample information
group <- as.factor(labels)
dlist$samples$group <- group
dlist$counts <- trunc(dlist$counts)

#Filter low-expressed genes 
#cutoff <- 1
#drop <- which(apply(cpm(dlist), 1, max) < cutoff)
#dlist <- dlist[-drop,] 
keep.exprs <- filterByExpr(dlist, group=group)
dlist <- dlist[keep.exprs,, keep.lib.sizes=FALSE]
cat("\nFiltered low-expressed genes:\n")
dim(dlist)

#Calculate normalization factors using TMM
dlist <- calcNormFactors(dlist, method = "TMM")
cat("\nNormalised gene expression:\n")
dlist$samples$norm.factors

#boxplot
x2 <- calcNormFactors(dlist)  
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm)

#MDS plot
plotMDS(dlist, col = as.numeric(group))

#Creating a design matrix and contrasts
##model matrix
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))
##Comparison between tumor and nontumor
contr.matrix <- makeContrasts(tumorvsnontumor = tumor-nontumor, levels = colnames(design))
cat("\nContrast matrix:\n")
contr.matrix


#Voom transformation
v <- voom(dlist, design, plot=TRUE)
#Fitting linear models in limma
vfit <- lmFit(v, design)
#Estimate contrast for each gene
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
#Empirical Bayes smoothing of standard errors
efit <- eBayes(vfit)

#list of DE genes
top.genes <- topTable(efit, sort.by = "P", n = Inf)
cat("\nTop DE genes:\n")
head(top.genes, 20)
top.genes$Gene <- rownames(top.genes)
top.genes <- top.genes[,c("Gene", names(top.genes)[1:6])]
write.table(top.genes, file = "DE_tumor_v_nontumor.txt", row.names = F, sep = "\t", quote = F)

#number of DE genes
cat("\nNumber of DE genes:\n")
length(which(top.genes$adj.P.Val < 0.05))
#differential expression levels
cat("\nDE levels:\n")
summary(decideTests(efit))
