library(GEOquery) 
library(rlang)
library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(doParallel)
library(shiny)
library(shinyMethyl)

registerDoParallel(cores = 4)

idat.folder <- "C:/Users/Miquel/OneDrive/TFG/GSE65057"
targets <- read.metharray.sheet(base=idat.folder)

#loading data
rgset <- read.metharray.exp(targets = targets)
phenoData <- rgset$Sample_Group
manifest <- getManifest(rgset)

MSet <- preprocessRaw(rgset) 

#passem de M values i U values a M values i Beta values
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
beta <- getBeta(RSet)
M <- getM(MSet)

#GenomicRatioSet: annotem les sondes
GRset <- mapToGenome(RSet)

beta <- getBeta(GRset)
M <- getM(GRset)
CN <- getCN(GRset)

sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)

gr <- granges(GRset)
head(gr, n= 3)

##Full annotation
annotation <- getAnnotation(GRset)
names(annotation)

###Normalization
gRatioSet.quantile <- preprocessQuantile(rgset)##SQN
summary <- shinySummarize(rgset)
summary.norm <- shinySummarize(gRatioSet.quantile)
runShinyMethyl(summary,summary.norm)
#####anotation
ann450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

keep <- !(featureNames(gRatioSet.quantile) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)

gRatioSet.quantile <- gRatioSet.quantile[keep,]
## remove probes with SNPs at CpG or SBE site
gRatioSet.quantile<- dropLociWithSnps(gRatioSet.quantile)

betas<-getBeta(gRatioSet.quantile)
M <- getM(gRatioSet.quantile)
#save(betas,file="betas.rda")
#load("betas.rda")
ann450k<-as.data.frame(ann450k)

## We will take care only of CpGs that are in promoter and gene body
ann450k <-ann450k [grep("TSS1500|TSS200|5'UTR|1stExon|Body",ann450k$UCSC_RefGene_Group),]
betas<-betas[rownames(betas) %in% rownames(ann450k) ,]
ann450k <-ann450k[rownames(betas) %in% rownames(ann450k) ,]

annotation2 <- data.frame(row = 1:length(ann450k$UCSC_RefGene_Name),
                          pos = ann450k$pos,
                          site=rownames(ann450k),
                          chr=ann450k$chr,
                          gene = ann450k$UCSC_RefGene_Name,
                          stringsAsFactors = F)


#DMA
#library(MEAL)
variable.diab<-as.factor(c("ND", rep("D", 2), rep("ND", 4), "D", "ND", "D", rep("ND", 2), rep("D", 5), rep("ND", 7)))
variable.obess<-as.factor(c(rep("O", 17), rep("NO",7)))
design<-model.matrix(~0+variable.diab+variable.obess)
colnames(design)<-c("D", "ND", "O")
rownames(design)<-sampleNames
contrast.matrix<-makeContrasts(D-ND,levels=design)

# limma 

#DMP<-runPipeline(set=gRatioSet.quantile, variable_names = "Sample_Group", covariable_names = "Sample_covar", betas = TRUE, model=~0+variable.diab+variable.obess)

eset <- ExpressionSet(assayData=betas)
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topgenes<-topTable(fit2, coef=1, adjust="BH") 
topgenes


logFCthreshold<-0.01
pvaluethreshold<-1.1e-04
topcpg<-topgenes[(abs(topgenes$logFC)>logFCthreshold)&(topgenes$P.Val<pvaluethreshold),]
topcpg

volcanoplot(fit2,coef=1,highlight=nrow(topcpg),names=rownames(betas),main="Diabetes vs normal")

tapply(betas["cg25653204",], variable.diab, summary)

boxplot(betas["cg25653204",]~variable.diab)

annotation[rownames(topcpg),"chr"]
annotation[rownames(topcpg),"UCSC_RefGene_Name"]
ann450k[rownames(topcpg),"UCSC_RefGene_Accession"]
