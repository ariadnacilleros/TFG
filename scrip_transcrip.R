#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("digest")

library(GEOquery)
library(digest)
library(oligo)

setwd("C:/Users/Miquel/OneDrive/TFG/GSE64998")
gseid<-"GSE64998"

list.files()

celFiles = list.files( pattern = "CEL.gz")

GEOFS <- read.celfiles(celFiles)

#exploration
exprs(GEOFS)

#Quality control
image(GEOFS)
colos<-rainbow(21)
hist(GEOFS,target="core",col=colos)
boxplot(GEOFS,target="core",col=colos,las=3,cex.axis=0.5)

#NORMALITZACIÓ
GEOFS.rma<-rma(GEOFS)
dim(exprs(GEOFS.rma)) #21 mostres i 33297 sondes
GEOFS.rma
class(GEOFS.rma)

boxplot(GEOFS.rma, col=colos, target="core")
hist(GEOFS.rma, col=colos, target="core")

#canviem nom mostres 

sampleNames(GEOFS.rma)
samplenames<-c("ND-obese1", "ND-obese2", "ND-obese3", "ND-obese4", "ND-obese5", "ND-obese6", "ND-obese7", "ND-obese8", "T2D-obese1", "T2D-obese2", "T2D-obese3", "T2D-obese4", "T2D-obese5", "T2D-obese6", "T2D-obese7", "non-obese1", "non-obese2", "non-obese3", "non-obese4", "non-obese5", "non-obese6")
sampleNames(GEOFS.rma)<-samplenames

#DATA AGREGATION
library(factoextra)

x<-exprs(GEOFS.rma)

summary(pca.filt <- prcomp(t(x), cor = TRUE ))

cond<-as.factor(c(rep("ND-O",8),rep("D-O",7), rep("ND-NO", 6)))

fviz_pca_ind(pca.filt, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = cond, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Llegenda") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

covar<-as.factor(c(rep("obese", "15"), rep("non-obese", 6)))


#DEG

library(limma)


variable.diab<-as.factor(c(rep("ND",8), rep("D", 7), rep("ND", 6)))
variable.obess<-as.factor(c(rep("O", 15), rep("NO",6)))
design<-model.matrix(~0+variable.diab+variable.obess)
colnames(design)<-c("D", "ND", "O")
rownames(design)<-samplenames
contrast.matrix<-makeContrasts(D-ND,levels=design)


fit<-lmFit(GEOFS.rma,design) #model could also be fitted to the batch corrected data
contrast.matrix<-makeContrasts(D-ND,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fite<-eBayes(fit2)
top.table<-topTable(fite,coef=1,number=Inf,adjust="BH")
results<-decideTests(fite)
table(results) 

#distribució p-valor 
hist(top.table$P.Value,breaks=100,main="results P")

results.p0.05<-top.table[top.table$P.Value<0.05,]
dim(results.p0.05) #1294 sondes

results.p0.05.logFC1<-top.table[top.table$P.Value<0.05 & abs(top.table$logFC)>1,]
dim(results.p0.05.logFC1) #4 sondes

###volcano 
volcanoplot(fite,coef=1,highlight=nrow(results.p0.05.logFC1),names=rownames(results.p0.05.logFC1),main="Diabetes vs normal")


#ANNOTATION
library(annotate)
library(hugene10sttranscriptcluster.db)
hugene10sttranscriptcluster() #conté tota la info sobre els identificadors dels gens

dat <- exprs(GEOFS.rma)[rownames(results.p0.05.logFC1),] #extreure la info del meu objecte expressionSet
logFC <- results.p0.05.logFC1$logFC
pval <- results.p0.05.logFC1$P.Value
adj.pval<-results.p0.05.logFC1$adj.P.Val

sym<-mget(rownames(results.p0.05.logFC1), env=hugene10sttranscriptclusterSYMBOL) #extrect el nom el símbol del gen
name<-mget(rownames(results.p0.05.logFC1), env=hugene10sttranscriptclusterGENENAME) #extrect el nom del gen
chr<-mget(rownames(results.p0.05.logFC1), env=hugene10sttranscriptclusterCHR) #extrect el cromosome del gen

affyids<-rownames(results.p0.05.logFC1)
genelist <- list(affyids)
filename <- "Results.html"
title <- "Differentially expressed genes at p<0.05 and abs logFC1"
othernames <- list(sym,name,chr,round(logFC, 1), round(pval, 4), round(adj.pval, 4), round(dat, 2)) 
head <- c("Probe ID", "Symbol", "Gene Name", "Chr","logFC", "p-value","adj.p-value",sampleNames(GEOFS.rma))
repository <- list("affy")
htmlpage(genelist, filename, title, othernames, head, repository = repository)

#volcano 
volcanoplot(fite,coef=1,highlight=nrow(topcpg),names=rownames(betas),main="Diabetes vs normal")


#METLYOME GENES
dat <- exprs(GEOFS.rma)[rownames(top.table),] #extreure la info del meu objecte expressionSet
logFC <- top.table$logFC
pval <- top.table$P.Value
adj.pval<-top.table$adj.P.Val

sym<-mget(rownames(top.table), env=hugene10sttranscriptclusterSYMBOL) #extrect el nom el símbol del gen
name<-mget(rownames(top.table), env=hugene10sttranscriptclusterGENENAME) #extrect el nom del gen
chr<-mget(rownames(top.table), env=hugene10sttranscriptclusterCHR) #extrect el cromosome del gen

affyids<-rownames(top.table)
genelist <- list(affyids)
filename <- "ResultsTOTAL.html"
title <- "Differentially expressed genes at p<0.05 and abs logFC1"
othernames <- list(sym,name,chr,round(logFC, 1), round(pval, 4), round(adj.pval, 4), round(dat, 2)) 
head <- c("Probe ID", "Symbol", "Gene Name", "Chr","logFC", "p-value","adj.p-value",sampleNames(GEOFS.rma))
repository <- list("affy")
htmlpage(genelist, filename, title, othernames, head, repository = repository)
head(othernames)
