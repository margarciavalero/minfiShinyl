### minfi pipeline for minfiShinyl web App ###


##Loading libraries
library(shiny)
library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggfortify)
library(Hmisc)
library(ComplexHeatmap)
library(circlize) 
library(DT)
library(dplyr)
library("IlluminaHumanMethylationEPICmanifest")
library(pvclust)


##Loading data
baseDir <- "C:/Users/alumne/Desktop/TFM/data/shinyApp"

input <- read.csv(file.path(baseDir, "SampleSheet_breast.csv"), sep=",")

targets <- read.metharray.sheet(as.character(input$appPath[1]))
RGset <- read.metharray.exp(targets=targets, force=TRUE)

pd <- pData(RGset)


## fist samples quality control
qc <- getQC(MSet)
rownames(qc) <- pd$Sample_Name
plotQC(qc)


## minfi objects normalization
detP <-detectionP(RGset)

keep <- colMeans(detP)<0.05

RGset <- RGset[,keep]

targets <- targets[keep,]

MSet.QNorm <- preprocessQuantile(RGset)

#probes in the same order
detP <- detP[match(featureNames(MSet.QNorm), rownames(detP)),]

#removing probes failed in one or more samples
keep <- rowSums(detP<0.01) == ncol(MSet.QNorm)

MSet.QNorm.filt <- MSet.QNorm[keep,]

#removing probes on the sex chromosome

annotEPIC <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

keep <- !(featureNames(MSet.QNorm.filt) %in% annotEPIC$Name[annotEPIC$chr %in% c("chrX", "chrY")])

MSet.QNorm.filt <- MSet.QNorm.filt[keep,]

#removing probes with SNPs at CpG site:

MSet.QNorm.filt <- dropLociWithSnps(MSet.QNorm.filt)


#exclude cross-reactive probes and polymorphic probes-> El ususario debe tener este documento! o debe estar subido a algún sitio para que el pipeline lo lea automáticamente desde dónde sea

reactive.probes <- read.table(file="C:/Users/alumne/Desktop/TFM/data/shinyApp/reactive_probes.txt")

keep <- !(featureNames(MSet.QNorm.filt) %in% reactive.probes$V1)

MSet.QNorm.filt <- MSet.QNorm.filt[keep,]

#GRset Normalization
GRset.funnorm <- preprocessFunnorm(RGset)



## Beta and M values calculation
mVals <- getM(MSet.QNorm.filt)
betaVals <- getBeta(MSet.QNorm.filt)
colnames(betaVals) <- pd$Sample_Name
colnames(mVals) <- pd$Sample_Name



## Probes Density Plot
densityPlot(getBeta(MSet.QNorm.filt), sampGroups = pd$status)



## PCA
pcs <- prcomp(t(getM(MSet.QNorm.filt)))
rownames(pcs$x) <- pd$Sample_Name
autoplot(pcs, label=T)



## The Heatmap: hierarchical analysis
Heatmap(t(betaVals[1:500, 1:6]),
        name = "CpGs beta values",
        column_title_gp = gpar(fontsize=9),
        col=colorRamp2(c(0, 0.5, 1),c("blue", "white", "red")),
        show_column_names = F,
        cluster_rows = T,
        cluster_columns = T,
        column_dend_side = "top",
        #column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 4),
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        #top_annotation = ha,
        column_title = "Heatmap of CpGs beta values in BR cohort"
)



## The Cluster: hierarchical analysis
fit <- pvclust(betaVals[1:500, 1:6], method.hclust="ward",
               method.dist="euclidean"
)
plot(fit) # dendogram with p values



## dmp calculation
dmp <- dmpFinder(betaVals, pheno = pd$status, type = "categorical")

dmpX <- dmp %>% mutate(
  X = rownames(dmp)
)

cpg.dmp.sign05 <- dmpX %>% filter(
  dmp$pval < 0.05
)
cpg.dmp.sign05.sort <- cpg.dmp.sign05[
  order(-cpg.dmp.sign05$pval), , drop=F
  ]


## dmp genes annotation
annotEPIC <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotEPIC.df <- as.data.frame(annotEPIC)

dmps.annot <- merge(
  cpg.dmp.sign05.sort,
  annotEPIC.df,
  by.x = "X",
  by.y = "Name"
)

genes.list.dmp.annot <- dmps.annot$UCSC_RefGene_Name



## DMR calculation
designMatrix <- model.matrix(~ pd$status)

dmrs <- bumphunter(GRset.funnorm, design = designMatrix, 
                   coef=1, 
                   cutoff=0.9,
                   type="Beta")

head(dmrs$table)
head(dmrs$coef)
