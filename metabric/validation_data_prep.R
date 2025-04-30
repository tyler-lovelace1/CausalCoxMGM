library(Biobase)
library(GEOquery)
library(limma)
library(huge)
library(survival)
library(survminer)
library(affycoretools)
library(hgu133plus2.db)
library(oligo)
library(sva)
library(umap)
library(parallel)
library(tidyverse)

readNormGeneAveCEL <- function(gseid) {

    gset <- getGEO(gseid, destdir = paste0('external_validation/', gseid), GSEMatrix=TRUE)

    geneData.list <- list()

    if (length(gset)==1) {
        celfiles <- c(list.files(paste0('external_validation/',gseid,'/rawCEL'),
                                 pattern='*.CEL.gz', full.names=T),
                      list.files(paste0('external_validation/',gseid,'/rawCEL'),
                                 pattern='*.cel.gz', full.names=T))
        
        celData <- read.celfiles(filenames=celfiles)

        oligo::hist(celData, main='Raw Probe Level')

        rmaData <- oligo::rma(celData)

        if (celData@annotation=='pd.hg.u133.plus.2') {
            require("hgu133plus2.db")
            rmaData <- annotateEset(rmaData, hgu133plus2.db, columns = c("PROBEID", "SYMBOL", "GENENAME"), multivals = "first")
        } else if (celData@annotation=='pd.hg.u133a') {
            require("hgu133a.db")
            rmaData <- annotateEset(rmaData, hgu133a.db, columns = c("PROBEID", "SYMBOL", "GENENAME"), multivals = "first")
        } else if (celData@annotation=='pd.hg.u133b') {
            require("hgu133b.db")
            rmaData <- annotateEset(rmaData, hgu133b.db, columns = c("PROBEID", "SYMBOL", "GENENAME"), multivals = "first")
        }

        geoid <- gset[[1]]@phenoData@data$geo_accession

        for (jdx in 1:length(geoid)) {
            colnames(rmaData)[grepl(geoid[jdx], colnames(rmaData), ignore.case=TRUE)] <- geoid[jdx]
        }
        
        rmaData <- rmaData[!is.na(rmaData@featureData@data$SYMBOL),]

        oligo::hist(rmaData, main='Normalized Probe Level')

        ex.ave <- avereps(exprs(rmaData), ID=rmaData@featureData@data[['SYMBOL']])

        featData <- AnnotatedDataFrame(rmaData@featureData@data %>% dplyr::group_by(SYMBOL) %>% dplyr::filter(row_number()==1) %>% dplyr::select(-c('PROBEID')) %>% data.frame(row.names=rownames(ex.ave)))

        geneData <- ExpressionSet(assayData=ex.ave, phenoData=phenoData(rmaData), featureData=featData)

        oligo::hist(geneData, main='Gene Level')

        annotation(geneData) <- annotation(gset[[1]])

        geneData.list[[annotation(geneData)]] <- geneData
        
    } else {

        platforms <- sapply(gset, annotation)

        celfiles <- c(list.files(paste0('external_validation/',gseid,'/rawCEL'),
                                 pattern='*.CEL.gz', full.names=T),
                      list.files(paste0('external_validation/',gseid,'/rawCEL'),
                                 pattern='*.cel.gz', full.names=T))

        ids <- gsub('[.]cel[.]gz', '',
                    gsub('[.]CEL[.]gz', '',
                         gsub(paste0('external_validation/', gseid, '/rawCEL/'), '', celfiles)))

        ids <- toupper(ids)

        ## idx <- 1

        for (pl in platforms[platforms!='GPL97']) {
        
            celData <- read.celfiles(filenames=celfiles[ids %in% colnames(gset[[which(platforms==pl)]])])
            
            oligo::hist(celData, main='Raw Probe Level')

            rmaData <- oligo::rma(celData)

            if (celData@annotation=='pd.hg.u133.plus.2') {
                require("hgu133plus2.db")
                rmaData <- annotateEset(rmaData, hgu133plus2.db, columns = c("PROBEID", "SYMBOL", "GENENAME"), multivals = "first")
            } else if (celData@annotation=='pd.hg.u133a') {
                require("hgu133a.db")
                rmaData <- annotateEset(rmaData, hgu133a.db, columns = c("PROBEID", "SYMBOL", "GENENAME"), multivals = "first")
            } else if (celData@annotation=='pd.hg.u133b') {
                require("hgu133b.db")
                rmaData <- annotateEset(rmaData, hgu133b.db, columns = c("PROBEID", "SYMBOL", "GENENAME"), multivals = "first")
            }

            geoid <- gset[[which(platforms==pl)]]@phenoData@data$geo_accession

            for (jdx in 1:length(geoid)) {
                colnames(rmaData)[grepl(geoid[jdx], colnames(rmaData), ignore.case=TRUE)] <- geoid[jdx]
            }
                        
            rmaData <- rmaData[!is.na(rmaData@featureData@data$SYMBOL),]

            oligo::hist(rmaData, main='Normalized Probe Level')

            ex.ave <- avereps(exprs(rmaData), ID=rmaData@featureData@data[['SYMBOL']])

            featData <- AnnotatedDataFrame(rmaData@featureData@data %>% dplyr::group_by(SYMBOL) %>% dplyr::filter(row_number()==1) %>% dplyr::select(-c('PROBEID')) %>% data.frame(row.names=rownames(ex.ave)))

            geneData <- ExpressionSet(assayData=ex.ave, phenoData=phenoData(rmaData), featureData=featData)

            oligo::hist(geneData, main='Gene Level')

            annotation(geneData) <- pl

            geneData.list[[annotation(geneData)]] <- geneData

            ## idx <- idx + 1
        }
    }

    return(geneData.list)
}

#### Exclude GSE11121: No value for age

ma.list <- c('GSE19615', 'GSE3494', 'GSE42568', 'GSE45255', 'GSE6532', 'GSE7390', 'GSE9195')

geneData.list <- lapply(ma.list, readNormGeneAveCEL)
names(geneData.list) <- ma.list

lapply(geneData.list, colnames)

metabric <- read.table('data_expression_median.txt', sep='\t', row.names=1, header=T)

metabric.sds <- apply(metabric, 1, sd, na.rm=T)

dim(metabric)

scanb <- read.csv('external_validation/scan-b/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv', row.names=1)

scanb <- scanb[,!grepl('repl', colnames(scanb))]

scanb.sds <- apply(scanb, 1, sd, na.rm=T)

names(scanb.sds[scanb.sds>0])

common.genes <- intersect(names(metabric.sds[metabric.sds!=0]), names(scanb.sds[scanb.sds!=0]))

for (gseid in ma.list) {
    common.genes <- intersect(common.genes, rownames(geneData.list[[gseid]]))
}

length(common.genes)

metabric <- metabric[common.genes,]
metabric <- metabric[,colSums(is.na(metabric))==0]

scanb <- scanb[common.genes,]

for (gseid in ma.list) {
    geneData.list[[gseid]] <- geneData.list[[gseid]][common.genes,]
}

clin.pat <- read.delim('../data_clinical_patient.txt', sep='\t', row.names=1, header=T, skip=4, na.strings=c('NA', ''))

rownames(clin.pat) <- make.names(rownames(clin.pat))

clin.pat$LYMPH_NODE_STATUS <- ifelse(clin.pat$LYMPH_NODES_EXAMINED_POSITIVE>0, 'Positive', 'Negative')

clin.samp <- read.delim('../data_clinical_sample.txt', sep='\t', row.names=1, header=T, skip=4, na.strings=c('NA', ''))

rownames(clin.samp) <- make.names(rownames(clin.samp))

clin.samp

clin.metabric <- cbind(clin.samp[colnames(metabric),], clin.pat[colnames(metabric),])

dim(clin.metabric)

clinList <- list()
clinList[['METABRIC']] <- na.omit(clin.metabric[,c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE', 'ER_STATUS', 'LYMPH_NODE_STATUS')])


metabric <- metabric[common.genes,-1]
metabric <- metabric[,colSums(is.na(metabric))==0]

joint.expr <- metabric[common.genes,-1]
joint.expr <- joint.expr[,colSums(is.na(joint.expr))==0]
dim(joint.expr)

batchVec <- rep('METABRIC', ncol(metabric[,colSums(is.na(metabric[common.genes,]))==0]))

pc.plot <- prcomp(t(scale(t(joint.expr))))$rotation[,1:30]

pc.plot <- as.data.frame(pc.plot)

pdf('individual_dataset_pca_plots.pdf', width=5, height=5)

ggplot(pc.plot, aes(x=PC1, y=PC2)) +
    geom_point() +
    geom_hline(aes(yintercept=mean(PC2)-6*sd(PC2))) +
    geom_hline(aes(yintercept=mean(PC2)+6*sd(PC2))) +
    geom_vline(aes(xintercept=mean(PC1)-6*sd(PC1))) +
    geom_vline(aes(xintercept=mean(PC1)+6*sd(PC1))) + 
    theme_bw() +
    ggtitle('METABRIC')

joint.expr <- cbind.data.frame(joint.expr, scanb[common.genes,])
## joint.expr <- joint.expr[,colSums(is.na(joint.expr))==0]
dim(joint.expr)

batchVec <- c(batchVec, rep('SCAN-B', ncol(scanb)))

pc.plot <- prcomp(t(scale(t(scanb[common.genes,]))))$rotation[,1:30]

pc.plot <- as.data.frame(pc.plot)

ggplot(pc.plot, aes(x=PC1, y=PC2)) +
    geom_point() +
    geom_hline(aes(yintercept=mean(PC2)-6*sd(PC2))) +
    geom_hline(aes(yintercept=mean(PC2)+6*sd(PC2))) +
    geom_vline(aes(xintercept=mean(PC1)-6*sd(PC1))) +
    geom_vline(aes(xintercept=mean(PC1)+6*sd(PC1))) + 
    theme_bw() +
    ggtitle('SCAN-B')

for (gseid in ma.list) {
    ## ex <- exprs(geneData.list[[gseid]])[common.genes,]

    ## pc.plot <- prcomp(t(scale(t(ex))))$rotation[,1:30]
    
    ## pc.plot <- as.data.frame(pc.plot)

    ## print(
    ##     ggplot(pc.plot, aes(x=PC1, y=PC2)) +
    ##     geom_point() +
    ##     geom_hline(aes(yintercept=mean(PC2)-6*sd(PC2))) +
    ##     geom_hline(aes(yintercept=mean(PC2)+6*sd(PC2))) +
    ##     geom_vline(aes(xintercept=mean(PC1)-6*sd(PC1))) +
    ##     geom_vline(aes(xintercept=mean(PC1)+6*sd(PC1))) + 
    ##     theme_bw() +
    ##     ggtitle(gseid)
    ## )
    
    ## joint.expr <- cbind.data.frame(joint.expr, ex)
    batchVec <- c(batchVec, rep(gseid, ncol(geneData.list[[gseid]])))
}

dev.off()

dim(joint.expr)

pc.plot <- prcomp(t(scale(t(joint.expr))))$rotation[,1:30]

pc.plot <- as.data.frame(pc.plot)

pc.plot$Batch <- batchVec

pdf('all_datasets_pca_uncorrected.pdf', width=5.5, height=5)
ggplot(pc.plot, aes(x=PC1, y=PC2, color=Batch)) +
    geom_point() +
    theme_bw()
dev.off()

plotDensities(metabric[,-1])

combat.expr <- ComBat(joint.expr, batchVec, ref.batch='METABRIC')


pc.plot <- prcomp(t(scale(t(combat.expr))))$rotation[,1:30]

pc.plot <- as.data.frame(pc.plot)

pc.plot$Batch <- batchVec

pdf('all_datasets_pca_combat.pdf', width=5.5, height=5)
ggplot(pc.plot, aes(x=PC1, y=PC2, color=Batch)) +
    geom_point() +
    theme_bw()
dev.off()


ump <- umap(pc.plot[,1:30], n_neighbors = 30, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)


ggplot(cbind.data.frame(UMAP=ump$layout, Batch=pc.plot$Batch)[apply(ump$layout, 1, max)<10,],
       aes(x=UMAP.1, y=UMAP.2, color=Batch)) +
    geom_point() +
    scale_color_brewer(palette='Set3') +
    theme_bw()


plotDensities(metabric[common.genes,-1])

oligo::boxplot(joint.expr)
oligo::boxplot(combat.expr)

plotDensities(joint.expr)
plotDensities(combat.expr)


## ERsplit.metabric <- read.csv('data/ERsplit_metabric_survival.csv', header=T)

## ERsplit.genes <- colnames(metabric)[1:(ncol(metabric)-17)]

#### GSE19615 DFHCC cohort ####
##
##  Li Y, Zou L, Li Q, Haibe-Kains B et al. Amplification of LAPTM4B
##  and YWHAZ contributes to chemotherapy resistance and recurrence of
##  breast cancer. Nat Med 2010 Feb;16(2):214-8. PMID: 20098429

geneData <- readNormGeneAveCEL('GSE19615')

celData <- read.celfiles(filenames=list.files('external_validation/GSE19615/rawCEL', pattern='*.CEL.gz', full.names=T))

celData

oligo::boxplot(celData)

oligo::hist(celData)

rmaData <- oligo::rma(celData)

rmaData

oligo::boxplot(rmaData)

oligo::hist(rmaData)

rmaData <- annotateEset(rmaData, hgu133plus2.db, columns = c("PROBEID", "SYMBOL", "GENENAME"), multivals = "first")

dim(rmaData)

rmaData <- rmaData[!is.na(rmaData@featureData@data$SYMBOL),]

dim(rmaData)

## gset <- getGEO("GSE19615", destdir = 'external_validation/GSE19615', GSEMatrix =TRUE)

## gset <- gset[[1]]

ex <- exprs(rmaData)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste("GSE19615", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE19615")

## ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(prcomp(t(scale(t(ex))))$rotation[,1:30], n_neighbors = 5, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

ex.ave <- avereps(ex, ID=rmaData@featureData@data[['SYMBOL']])

dim(ex.ave)

## ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

dim(ex.ave)

title <- paste("GSE19615", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE19165")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

common.genes <- intersect(common.genes, rownames(ex.ave))


#### GSE42568 Dublin cohort ####
##
##  Clarke C, Madden SF, Doolan P, Aherne ST et al. Correlating
##  transcriptional networks to breast cancer survival: a large-scale
##  coexpression analysis. Carcinogenesis 2013
##  Oct;34(10):2300-8. PMID: 23740839

gset <- getGEO("GSE42568", destdir = 'external_validation/GSE42568', GSEMatrix =TRUE)

gset <- gset[[1]]

geneData <- readNormGeneAveCEL('GSE42568')

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE42568", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE42568")

## ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(prcomp(t(scale(t(ex))))$rotation[,1:30], n_neighbors = 10, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

ex.ave <- avereps(ex, ID=rmaData@featureData@data[['SYMBOL']])

dim(ex.ave)

## ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

dim(ex.ave)

title <- paste("GSE42568", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE42568")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

common.genes <- intersect(common.genes,
                          rownames(ex.ave))



#### GSE9195 Guyt2 cohort ####
##
##  Loi S, Haibe-Kains B, Desmedt C, Wirapati P et al. Predicting
##  prognosis using molecular profiling in estrogen receptor-positive
##  breast cancer treated with tamoxifen. BMC Genomics 2008 May
##  22;9:239. PMID: 18498629
##
##  Loi S, Haibe-Kains B, Majjaj S, Lallemand F et al. PIK3CA
##  mutations associated with gene signature of low mTORC1 signaling
##  and better outcomes in estrogen receptor-positive breast
##  cancer. Proc Natl Acad Sci U S A 2010 Jun
##  1;107(22):10208-13. PMID: 20479250

gset <- getGEO("GSE9195", destdir = 'external_validation/GSE9195', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste("GSE9195", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE9195")

## ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

dim(ex.ave)
ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]
dim(ex.ave)

title <- paste("GSE45255", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE45255")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

common.genes <- intersect(common.genes,
                          unlist(strsplit(rownames(ex.ave), ' /// ')))
length(common.genes)

#### GSE11121 Maintz cohort ####
##
## Schmidt M, Böhm D, von Törne C, Steiner E et al. The humoral immune
## system has a key prognostic impact in node-negative breast
## cancer. Cancer Res 2008 Jul 1;68(13):5405-13. PMID: 18593943
##
## Cadenas C, van de Sandt L, Edlund K, Lohr M et al. Loss of
## circadian clock gene expression is associated with tumor
## progression in breast cancer. Cell Cycle 2014;13(20):3282-91. PMID:
## 25485508
##
## Hellwig B, Hengstler JG, Schmidt M, Gehrmann MC et al. Comparison
## of scores for bimodality of gene expression distributions and
## genome-wide evaluation of the prognostic relevance of high-scoring
## genes. BMC Bioinformatics 2010 May 25;11:276. PMID: 20500820
##
## Heimes AS, Härtner F, Almstedt K, Krajnak S et al. Prognostic
## Significance of Interferon-γ and Its Signaling Pathway in Early
## Breast Cancer Depends on the Molecular Subtypes. Int J Mol Sci 2020
## Sep 29;21(19). PMID: 33003293

gset <- getGEO("GSE11121", destdir = 'external_validation/GSE11121', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE11121", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE11121")

## ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

dim(ex.ave)
ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]
dim(ex.ave)

title <- paste("GSE11121", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE11121")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

common.genes <- intersect(common.genes,
                          unlist(strsplit(rownames(ex.ave), ' /// ')))
length(common.genes)


#### GSE6532 TAM cohort ####
##
## Loi S, Haibe-Kains B, Desmedt C, Lallemand F et al. Definition of
## clinically distinct molecular subtypes in estrogen
## receptor-positive breast carcinomas through genomic grade. J Clin
## Oncol 2007 Apr 1;25(10):1239-46. PMID: 17401012
##
## Loi S, Haibe-Kains B, Desmedt C, Wirapati P et al. Predicting
## prognosis using molecular profiling in estrogen receptor-positive
## breast cancer treated with tamoxifen. BMC Genomics 2008 May
## 22;9:239. PMID: 18498629
##
## Loi S, Haibe-Kains B, Majjaj S, Lallemand F et al. PIK3CA mutations
## associated with gene signature of low mTORC1 signaling and better
## outcomes in estrogen receptor-positive breast cancer. Proc Natl
## Acad Sci U S A 2010 Jun 1;107(22):10208-13. PMID: 20479250
##
## Loi S, Sotiriou C, Haibe-Kains B, Lallemand F et al. Gene
## expression profiling identifies activated growth factor signaling
## in poor prognosis (Luminal-B) estrogen receptor positive breast
## cancer. BMC Med Genomics 2009 Jun 24;2:37. PMID: 19552798

gset <- getGEO("GSE6532", destdir = 'external_validation/GSE6532', GSEMatrix =TRUE)

## gset <- gset[-1]

gset[[1]]$`samplename:ch1`[is.na(gset[[1]]$`samplename:ch1`)] <- gset[[1]]$`sample name:ch1`[is.na(gset[[1]]$`samplename:ch1`)]

ex96 <- exprs(gset[[1]])
ex97 <- exprs(gset[[2]])

colnames(ex96) <- gset[[1]]$`samplename:ch1`
colnames(ex97) <- gset[[2]]$`samplename:ch1`

ex <- rbind(ex96, ex97[,colnames(ex96)])


# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE6532", "/", annotation(gset[[1]]), "+", annotation(gset[[2]]), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

## ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE6532")

## ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

ex.ave <- avereps(ex, ID=c(gset[[1]]@featureData@data[['Gene Symbol']],
                           gset[[2]]@featureData@data[['Gene Symbol']]))

dim(ex.ave)
ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]
dim(ex.ave)

title <- paste("GSE6532", "/", annotation(gset[[1]]), "+", annotation(gset[[2]]), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE6532")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

common.genes <- intersect(common.genes,
                          unlist(strsplit(rownames(ex.ave), ' /// ')))
length(common.genes)

#### GSE7390 Transbig cohort ####
##
##  Loi S, Haibe-Kains B, Desmedt C, Wirapati P et al. Predicting
##  prognosis using molecular profiling in estrogen receptor-positive
##  breast cancer treated with tamoxifen. BMC Genomics 2008 May
##  22;9:239. PMID: 18498629
##
##  Loi S, Haibe-Kains B, Majjaj S, Lallemand F et al. PIK3CA
##  mutations associated with gene signature of low mTORC1 signaling
##  and better outcomes in estrogen receptor-positive breast
##  cancer. Proc Natl Acad Sci U S A 2010 Jun
##  1;107(22):10208-13. PMID: 20479250

gset <- getGEO("GSE7390", destdir = 'external_validation/GSE7390', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE7390", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE7390")

## ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

dim(ex.ave)
ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]
dim(ex.ave)

title <- paste("GSE7390", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE45255")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

common.genes <- intersect(common.genes,
                          unlist(strsplit(rownames(ex.ave), ' /// ')))
length(common.genes)


#### GSE45255 IRB/JNR/NUH cohort ####
##
##  Nagalla, S. et al. Interactions between immunity, proliferation
##  and molecular subtype in breast cancer prognosis. Genome Biol. 14,
##  R34 (2013).

gset <- getGEO("GSE45255", destdir = 'external_validation/GSE45255', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE45255", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE45255")

## ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

dim(ex.ave)
ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]
dim(ex.ave)

title <- paste("GSE45255", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE45255")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

common.genes <- intersect(common.genes,
                          unlist(strsplit(rownames(ex.ave), ' /// ')))
length(common.genes)

clin.GSE45255 <- gset@phenoData@data

clin.GSE45255 <- clin.GSE45255 %>% select(c('dss event (defined as death from breast cancer):ch1',
                                            'DSS time:ch1', 'er status:ch1',
                                            'ln status:ch1', 'patient age:ch1', 'size (mm):ch1'))
colnames(clin.GSE45255) <- c('event', 'survival', 'ER.Status',
                             'Lymph.Node.Status', 'Age.at.Diagnosis', 'Tumor.Size')

GSE45255 <- cbind(t(ex.ave),
                  clin.GSE45255)

tail(GSE45255)
row.has.na <- apply(GSE45255, 1, function(x){any(x=="NA")})
GSE45255 <- GSE45255[!row.has.na,]
tail(GSE45255)

dim(GSE45255)

GSE45255$Age.at.Diagnosis <- as.numeric(GSE45255$Age.at.Diagnosis)
GSE45255$Tumor.Size <- as.numeric(GSE45255$Tumor.Size)

## GSE45255 <- GSE45255 %>% mutate_if(function(x) { is.numeric(x) },
##                                    function(x) { (x-mean(x))/sd(x) })

GSE45255 <- GSE45255 %>% mutate_if(function(x) { is.numeric(x) },
                                   function(x) { huge.npn((matrix(x)-mean(x))/sd(x), npn.func='truncation') })

GSE45255$survival <- as.numeric(GSE45255$survival)
GSE45255$event <- as.numeric(GSE45255$event)

GSE45255$ER.Status[GSE45255$ER.Status=='ER+'] = 'Positive'
GSE45255$ER.Status[GSE45255$ER.Status=='ER-'] = 'Negative'

GSE45255$HER2.Status[GSE45255$HER2.Status=='He+'] = 'Positive'
GSE45255$HER2.Status[GSE45255$HER2.Status=='He-'] = 'Negative'

GSE45255$Lymph.Node.Status[GSE45255$Lymph.Node.Status=='LN+'] = 'NodePositive'
GSE45255$Lymph.Node.Status[GSE45255$Lymph.Node.Status=='LN-'] = 'NodeNegative'

ggsurvplot(survfit(Surv(survival, event) ~ ER.Status + Lymph.Node.Status, GSE45255))

write.table(GSE45255, 'external_validation/GSE45255/data/ERsplit_GSE45255_survival.csv', sep=',', row.names=F)



#### GSE3494 Upp Cohort 


gset <- getGEO("GSE3494", destdir = 'external_validation/GSE3494', GSEMatrix=TRUE)

gset96 <- gset[[1]]
annotation(gset96)

ex96 <- exprs(gset96)
# log2 transform
qx96 <- as.numeric(quantile(ex96, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx96[5] > 100) || (qx96[6]-qx96[1] > 50 && qx96[2] > 0)
if (LogC) {
    ex96[which(ex96 <= 0)] <- NaN
    ex96 <- log2(ex96)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE3494", "/", annotation(gset96), " value distribution", sep ="")
plotDensities(ex96, main=title, legend=F)

ex96 <- na.omit(ex96) # eliminate rows with NAs
plotSA(lmFit(ex96), main="Mean variance trend, GSE3494")


gset97 <- gset[[2]]
annotation(gset97)

ex97 <- exprs(gset97)
# log2 transform
qx97 <- as.numeric(quantile(ex97, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx97[5] > 100) || (qx97[6]-qx97[1] > 50 && qx97[2] > 0)
LogC
if (LogC) {
    ex97[which(ex97 <= 0)] <- NaN
    ex97 <- log2(ex97)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE3494", "/", annotation(gset97), " value distribution", sep ="")
plotDensities(ex97, main=title, legend=F)

ex97 <- na.omit(ex97) # eliminate rows with NAs
plotSA(lmFit(ex97), main="Mean variance trend, GSE3494")

sampMap <- read.csv('external_validation/GSE3494/samplePatientMap.csv', header=T, sep='\t')
rownames(sampMap) <- sampMap$GEO.Sample.Accession..
head(sampMap)

colnames(ex96) <- sampMap[colnames(ex96), 'Patient.ID']
colnames(ex97) <- sampMap[colnames(ex97), 'Patient.ID']

ex.joint <- rbind(ex96, ex97[,colnames(ex96)])

dim(ex.joint)

dim(ex.joint)

title <- paste ("GSE3494", "/", annotation(gset96), "+", annotation(gset97), " value distribution", sep ="")
plotDensities(ex.joint, main=title, legend=F)

ex.ave <- avereps(ex.joint, ID=c(gset96@featureData@data[['Gene Symbol']],
                                 gset97@featureData@data[['Gene Symbol']]))

dim(ex.ave)
ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]
dim(ex.ave)

title <- paste ("GSE3494", "/", annotation(gset96), "+", annotation(gset97), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)

plotSA(lmFit(ex.ave), main="Mean variance trend, GSE3494")

common.genes <- intersect(common.genes,
                          unlist(strsplit(rownames(ex.ave), ' /// ')))
length(common.genes)

clin.GSE3494 <- read.csv('external_validation/GSE3494/clinData.csv', header=T, sep='\t')
rownames(clin.GSE3494) <- clin.GSE3494$INDEX..ID.
head(clin.GSE3494)

clin.GSE3494 <- clin.GSE3494 %>% select(c('DSS.EVENT..Disease.Specific.Survival.EVENT..1.death.from.breast.cancer..0.alive.or.censored..',
                                          'DSS.TIME..Disease.Specific.Survival.Time.in.years.',
                                          'ER.status', 'PgR.status',
                                          'Lymph.node.status', 'age.at.diagnosis', 'tumor.size..mm.'))

colnames(clin.GSE3494) <- c('event', 'survival', 'ER.Status', 'PR.Status',
                            'Lymph.Node.Status', 'Age.at.Diagnosis', 'Tumor.Size')

GSE3494 <- cbind(t(ex.ave[unique(rownames(ex.ave)[rownames(ex.ave) %in% ERsplit.genes]),]), clin.GSE3494[colnames(ex.ave),])
head(GSE3494)


GSE3494$survival <- as.character(GSE3494$survival)
GSE3494$event <- as.character(GSE3494$event)
GSE3494$Age.at.Diagnosis <- as.character(GSE3494$Age.at.Diagnosis)
GSE3494$Tumor.Size <- as.character(GSE3494$Tumor.Size)

GSE3494[GSE3494$ER.Status == 'ER?','ESR1']

table(GSE3494$ER.Status) / nrow(GSE3494)

quantile(GSE3494$ESR1, c(0, 0.25, 0.5, 0.75, 1.0))

GSE3494[GSE3494$ER.Status == 'ER?','ER.Status'] = 'ER+'

GSE3494$ER.Status[GSE3494$ER.Status=='ER+'] = 'Positive'
GSE3494$ER.Status[GSE3494$ER.Status=='ER-'] = 'Negative'

GSE3494$PR.Status[GSE3494$PR.Status=='PgR+'] = 'Positive'
GSE3494$PR.Status[GSE3494$PR.Status=='PgR-'] = 'Negative'

GSE3494$Lymph.Node.Status[GSE3494$Lymph.Node.Status=='LN+'] = 'NodePositive'
GSE3494$Lymph.Node.Status[GSE3494$Lymph.Node.Status=='LN-'] = 'NodeNegative'
GSE3494$Lymph.Node.Status[GSE3494$Lymph.Node.Status=='LN?'] = "NA"

row.has.na <- apply(GSE3494, 1, function(x){any(x=="NA") || any(is.na(x))})
row.has.na
GSE3494 <- GSE3494[!row.has.na,]

head(GSE3494)
dim(GSE3494)

GSE3494$Age.at.Diagnosis <- as.numeric(GSE3494$Age.at.Diagnosis)
GSE3494$Tumor.Size <- as.numeric(GSE3494$Tumor.Size)

GSE3494 <- GSE3494 %>% mutate_if(function(x) { is.numeric(x) },
                                 function(x) { huge.npn((matrix(x)-mean(x))/sd(x), npn.func='truncation') })

head(GSE3494)

GSE3494$survival <- as.numeric(GSE3494$survival)
GSE3494$event <- as.numeric(GSE3494$event)

ggsurvplot(survfit(Surv(survival, event) ~ strata(ER.Status) + (ENC1>0) + (CCT6B>0) + Lymph.Node.Status, GSE3494))

write.table(GSE3494, 'external_validation/GSE3494/data/ERsplit_GSE3494_survival.csv', sep=',', row.names=F)



#### GSE96058 SCAN-B Cohort


gset <- getGEO("GSE96058", destdir = 'external_validation/scan-b', GSEMatrix=TRUE)

## gset <- gset[[1]]
## rownames(gset[[1]]@phenoData@data)
## rownames(gset[[2]]@phenoData@data) %in% rownames(gset[[1]]@phenoData@data)

## samps <- c(gset[[1]]@phenoData@data$geo_accession,
##            gset[[2]]@phenoData@data$geo_accession)

clin.GSE96058 <- rbind(gset[[1]]@phenoData@data, gset[[2]]@phenoData@data)

rownames(clin.GSE96058) <- clin.GSE96058$title

## clin.GSE96058 <- read.csv('external_validation/scan-b/file3ea7511f82fb.tsv', sep='\t', skip=19)

## head(clin.GSE96058$SAMPLE)
## rownames(clin.GSE96058) <- clin.GSE96058$SAMPLE

clin.GSE96058 <- clin.GSE96058 %>% dplyr::select(c('overall survival event:ch1', 'overall survival days:ch1',
                                            'er status:ch1', 'lymph node status:ch1',
                                            'age at diagnosis:ch1', 'tumor size:ch1',
                                            'pam50 subtype:ch1'))

colnames(clin.GSE96058) <- c('event', 'survival', 'ER.Status',
                             'Lymph.Node.Status', 'Age.at.Diagnosis', 'Tumor.Size',
                             'Pam50.Subtype')

ex <- read.table('external_validation/scan-b/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv', sep=',', row.names=1, header=T, colClasses=c('character', rep('numeric',3409)))[,rownames(clin.GSE96058)]

dim(ex)

## ex <- fread('external_validation/scan-b/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv', sep=',')

## ex <- vroom('external_validation/scan-b/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv', delim=',', col_select=rownames(clin.GSE96058))

## ex <- ex[ERsplit.genes,]
head(ex)

## ex <- ex[,colnames(ex) %in% rownames(clin.GSE96058)]
## ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE96058", " value distribution", sep ="")
plotDensities(ex[common.genes,], main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex[common.genes,]), main="Mean variance trend, GSE96058")

## ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(scale(t(ex)), n_neighbors = 20, random_state = 123)
plotdf <- data.frame(UMAP1 = ump$layout[,1],
                     UMAP2 = ump$layout[,2],
                     PAM50 = factor(clin.GSE96058$Pam50.Subtype))
ggplot(plotdf, aes(x=UMAP1, y = UMAP2, color=PAM50)) +
    geom_point() +
    theme_bw()

plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

## unique(gset@featureData@data[['Gene Symbol']][gset@featureData@data[['Gene Symbol']] %in% ERsplit.genes])

## ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

## dim(ex.ave)
## ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]
## dim(ex.ave)

## title <- paste ("GSE96058", "/", annotation(gset), " value distribution", sep ="")
## plotDensities(ex.ave, main=title, legend=F)

## head(gset$GSE96058_series_matrix.txt.gz@assayData$exprs)
## head(ex.ave[unique(gset@featureData@data[['Gene Symbol']][gset@featureData@data[['Gene Symbol']] %in% ERsplit.genes]),])



common.genes <- intersect(common.genes,
                          rownames(ex)[apply(ex, 1, sd)!=0])
length(common.genes)

saveRDS(common.genes, 'meta_cohort_common_genes.rds')


common.genes <- readRDS('meta_cohort_common_genes.rds')
length(common.genes)

hvar.genes <- rownames(read.csv('data/hvar_genes_500_hub_summarized.csv', row.names=1))

hvar.genes <- gsub('_grp', '', hvar.genes)

length(hvar.genes)

#### GSE19615 DFHCC cohort ####
##
##  Li Y, Zou L, Li Q, Haibe-Kains B et al. Amplification of LAPTM4B
##  and YWHAZ contributes to chemotherapy resistance and recurrence of
##  breast cancer. Nat Med 2010 Feb;16(2):214-8. PMID: 20098429

gset <- getGEO("GSE19615", destdir = 'external_validation/GSE19615', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

rownames(ex.ave) <- sapply(strsplit(rownames(ex.ave), ' /// '), function(x) unique(x))

hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex.ave <- ex.ave[hvar.genes,]

title <- paste("GSE19615", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE19165")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

## gene.subset <- rownames(ex.ave)[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) any(x %in% common.genes))]

## ex.ave <- ex.ave[gene.subset,]

gset <- getGEO("GSE19615", destdir = 'external_validation/GSE19615', GSEMatrix =TRUE)

clin <- gset[[1]]@phenoData@data %>% dplyr::select(`er:ch1`, `age:ch1`, `her.2:ch1`, `lymph nodes:ch1`, `tumor size (cm):ch1`, `distant recur (yn):ch1`, `distant recurrence free survival (mo):ch1`, `time of followup (mo):ch1`) %>% mutate_at(c('age:ch1', 'tumor size (cm):ch1', 'distant recurrence free survival (mo):ch1', 'time of followup (mo):ch1'), as.numeric)

clin$`distant recurrence free survival (mo):ch1` <- ifelse(!is.na(clin$`distant recurrence free survival (mo):ch1`), clin$`distant recurrence free survival (mo):ch1`, clin$`time of followup (mo):ch1`)

clin <- clin[,-8]

colnames(clin) <- c('ER_STATUS', 'AGE_AT_DIAGNOSIS', 'HER2_STATUS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'DRFS.status', 'DRFS.time')

clin$DRFS.time <- 365.25/12 * clin$DRFS.time

clin$ER_STATUS <- ifelse(grepl('pos', clin$ER_STATUS), 'Positive', 'Negative')
clin$HER2_STATUS <- ifelse(grepl('pos', clin$HER2_STATUS), 'Positive', 'Negative')
clin$LYMPH_NODE_STATUS <- ifelse(grepl('pos', clin$LYMPH_NODE_STATUS), 'Positive', 'Negative')
clin$DRFS.status <- as.integer(clin$DRFS.status=='Y')

clin$TUMOR_SIZE <- 10*clin$TUMOR_SIZE

ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ HER2_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

dim(na.omit(clin))

clinList[['GSE19615']] <- na.omit(clin)

GSE19615.hvar <- cbind(t(combat.expr[,rownames(clin)]), clin)

dim(GSE19615.hvar)
dim(na.omit(GSE19615.hvar))

write.csv(GSE19615.hvar, 'external_validation/GSE19615/GSE19615_corrected.csv')


#### GSE42568 Dublin cohort ####
##
##  Clarke C, Madden SF, Doolan P, Aherne ST et al. Correlating
##  transcriptional networks to breast cancer survival: a large-scale
##  coexpression analysis. Carcinogenesis 2013
##  Oct;34(10):2300-8. PMID: 23740839

gset <- getGEO("GSE42568", destdir = 'external_validation/GSE42568', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

ex <- na.omit(ex) # eliminate rows with NAs
## sapply(strsplit(gset@featureData@data[['Gene Symbol']], ' /// ')
ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

rownames(ex.ave) <- sapply(strsplit(rownames(ex.ave), ' /// '), function(x) unique(x))

hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex.ave <- ex.ave[hvar.genes,]


title <- paste("GSE42568", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE42568")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

## gene.subset <- rownames(ex.ave)[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) any(x %in% common.genes))]

## ex.ave <- ex.ave[gene.subset,]

clin <- gset@phenoData@data %>% dplyr::select(`er_status:ch1`, `age:ch1`, `lymph node status:ch1`, `size:ch1`, `overall survival event:ch1`, `overall survival time_days:ch1`, `relapse free survival event:ch1`, `relapse free survival time_days:ch1`) %>% mutate_at(c('er_status:ch1', 'age:ch1', 'size:ch1', 'overall survival time_days:ch1', 'relapse free survival time_days:ch1'), as.numeric)

clin <- na.omit(clin)

colnames(clin) <- c('ER_STATUS', 'AGE_AT_DIAGNOSIS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'OS.status', 'OS.time', 'DFS.status', 'DFS.time')

clin$ER_STATUS <- ifelse(grepl(1, clin$ER_STATUS), 'Positive', 'Negative')
clin$LYMPH_NODE_STATUS <- ifelse(grepl('1', clin$LYMPH_NODE_STATUS), 'Positive', 'Negative')
clin$DFS.status <- as.integer(clin$DFS.status=='1')
clin$OS.status <- as.integer(clin$OS.status=='1')
clin$TUMOR_SIZE <- 10 * clin$TUMOR_SIZE

ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

ggsurvplot(survfit(Surv(OS.time, OS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(OS.time, OS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

dim(na.omit(clin))

clinList[['GSE42568']] <- na.omit(clin)

GSE42568.hvar <- cbind(t(ex.ave), clin[colnames(ex.ave),])

dim(GSE42568.hvar)

write.csv(GSE42568.hvar, 'external_validation/GSE42568/GSE42568_hvar.csv')

GSE42568.hvar <- cbind(t(combat.expr[,rownames(clin)]), clin)

dim(GSE19615.hvar)
dim(na.omit(GSE19615.hvar))

write.csv(GSE19615.hvar, 'external_validation/GSE19615/GSE19615_corrected.csv')



#### GSE9195 Guyt2 cohort ####
##
##  Loi S, Haibe-Kains B, Desmedt C, Wirapati P et al. Predicting
##  prognosis using molecular profiling in estrogen receptor-positive
##  breast cancer treated with tamoxifen. BMC Genomics 2008 May
##  22;9:239. PMID: 18498629
##
##  Loi S, Haibe-Kains B, Majjaj S, Lallemand F et al. PIK3CA
##  mutations associated with gene signature of low mTORC1 signaling
##  and better outcomes in estrogen receptor-positive breast
##  cancer. Proc Natl Acad Sci U S A 2010 Jun
##  1;107(22):10208-13. PMID: 20479250

gset <- getGEO("GSE9195", destdir = 'external_validation/GSE9195', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE9195", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
## sapply(strsplit(gset@featureData@data[['Gene Symbol']], ' /// ')
ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

rownames(ex.ave) <- sapply(strsplit(rownames(ex.ave), ' /// '), function(x) unique(x))

hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex.ave <- ex.ave[hvar.genes,]


title <- paste("GSE9195", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE9195")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

clin <- gset@phenoData@data %>% dplyr::select(`er:ch1`, `pgr:ch1`, `age:ch1`, `node:ch1`, `size:ch1`, `e.dmfs:ch1`, `t.dmfs:ch1`, `e.rfs:ch1`, `t.rfs:ch1`) %>% mutate_at(c('er:ch1', 'pgr:ch1', 'age:ch1', 'node:ch1', 'size:ch1', 'e.dmfs:ch1', 't.dmfs:ch1', 'e.rfs:ch1', 't.rfs:ch1'), as.numeric)

clin <- na.omit(clin)

## er.gmm <- GMM(as.matrix(ex.ave['ESR1',rownames(clin)]), 2)
## er.gmm

## er.res <- glm(as.numeric(`er_status:ch1`) ~  ESR1, cbind(clin[clin$`er_status:ch1`!='NA',], data.frame(ESR1=ex.ave['ESR1',rownames(clin[clin$`er_status:ch1`!='NA',])])), family=binomial())

## er.pred <- round(predict(er.res, data.frame(ESR1=ex.ave['ESR1',rownames(clin)]), type='response'))

## table(er.pred, clin$`er_status:ch1`)

colnames(clin) <- c('ER_STATUS', 'PR_STATUS', 'AGE_AT_DIAGNOSIS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'DRFS.status', 'DRFS.time', 'DFS.status', 'DFS.time')

clin$ER_STATUS <- ifelse(grepl(1, clin$ER_STATUS), 'Positive', 'Negative')
clin$LYMPH_NODE_STATUS <- ifelse(grepl('1', clin$LYMPH_NODE_STATUS), 'Positive', 'Negative')
clin$DFS.status <- as.integer(clin$DFS.status=='1')
clin$DRFS.status <- as.integer(clin$DRFS.status=='1')

clin$TUMOR_SIZE <- 10 * clin$TUMOR_SIZE

ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

clinList[['GSE9195']] <- na.omit(clin)

GSE9195.hvar <- cbind(t(ex.ave), clin[colnames(ex.ave),])

dim(GSE9195.hvar)

write.csv(GSE9195.hvar, 'external_validation/GSE9195/GSE9195_hvar.csv')


#### GSE11121 Maintz cohort ####
##
## Schmidt M, Böhm D, von Törne C, Steiner E et al. The humoral immune
## system has a key prognostic impact in node-negative breast
## cancer. Cancer Res 2008 Jul 1;68(13):5405-13. PMID: 18593943
##
## Cadenas C, van de Sandt L, Edlund K, Lohr M et al. Loss of
## circadian clock gene expression is associated with tumor
## progression in breast cancer. Cell Cycle 2014;13(20):3282-91. PMID:
## 25485508
##
## Hellwig B, Hengstler JG, Schmidt M, Gehrmann MC et al. Comparison
## of scores for bimodality of gene expression distributions and
## genome-wide evaluation of the prognostic relevance of high-scoring
## genes. BMC Bioinformatics 2010 May 25;11:276. PMID: 20500820
##
## Heimes AS, Härtner F, Almstedt K, Krajnak S et al. Prognostic
## Significance of Interferon-γ and Its Signaling Pathway in Early
## Breast Cancer Depends on the Molecular Subtypes. Int J Mol Sci 2020
## Sep 29;21(19). PMID: 33003293

gset <- getGEO("GSE11121", destdir = 'external_validation/GSE11121', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE11121", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)


ex <- na.omit(ex) # eliminate rows with NAs
## sapply(strsplit(gset@featureData@data[['Gene Symbol']], ' /// ')
ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

rownames(ex.ave) <- sapply(strsplit(rownames(ex.ave), ' /// '), function(x) unique(x))

hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex.ave <- ex.ave[hvar.genes,]


title <- paste("GSE11121", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE11121")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

clin <- gset@phenoData@data %>% select(`node:ch1`, `size_in_cm:ch1`, `e.dmfs:ch1`, `t.dmfs:ch1`) %>% mutate_at(c('node:ch1', 'size_in_cm:ch1', 'e.dmfs:ch1', 't.dmfs:ch1'), as.numeric)

clin <- na.omit(clin)

## er.gmm <- GMM(as.matrix(ex.ave['ESR1',rownames(clin)]), 2)
## er.gmm

## er.res <- glm(as.numeric(`er_status:ch1`) ~  ESR1, cbind(clin[clin$`er_status:ch1`!='NA',], data.frame(ESR1=ex.ave['ESR1',rownames(clin[clin$`er_status:ch1`!='NA',])])), family=binomial())

## er.pred <- round(predict(er.res, data.frame(ESR1=ex.ave['ESR1',rownames(clin)]), type='response'))

## table(er.pred, clin$`er_status:ch1`)

colnames(clin) <- c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'DRFS.status', 'DRFS.time')

## clin$ER_STATUS <- ifelse(grepl(1, clin$ER_STATUS), 'Positive', 'Negative')
clin$LYMPH_NODE_STATUS <- ifelse(grepl('1', clin$LYMPH_NODE_STATUS), 'Positive', 'Negative')
clin$DRFS.status <- as.integer(clin$DRFS.status=='1')

ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

GSE11121.hvar <- cbind(t(ex.ave), clin[colnames(ex.ave),])

library(ClusterR)

hist(GSE11121.hvar$ESR1)

#### Thresholed to give the same number of ER- patients (44) as
#### reported in the original study
GSE11121.hvar$ER_STATUS <- ifelse(GSE11121.hvar$ESR1 > quantile(GSE11121.hvar$ESR1, probs=0.22), 'Positive', 'Negative')

ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ ER_STATUS, GSE11121.hvar), pval=T)

dim(GSE11121.hvar)

write.csv(GSE11121.hvar, 'external_validation/GSE11121/GSE11121_hvar.csv')


#### GSE6532 TAM cohort ####
##
## Loi S, Haibe-Kains B, Desmedt C, Lallemand F et al. Definition of
## clinically distinct molecular subtypes in estrogen
## receptor-positive breast carcinomas through genomic grade. J Clin
## Oncol 2007 Apr 1;25(10):1239-46. PMID: 17401012
##
## Loi S, Haibe-Kains B, Desmedt C, Wirapati P et al. Predicting
## prognosis using molecular profiling in estrogen receptor-positive
## breast cancer treated with tamoxifen. BMC Genomics 2008 May
## 22;9:239. PMID: 18498629
##
## Loi S, Haibe-Kains B, Majjaj S, Lallemand F et al. PIK3CA mutations
## associated with gene signature of low mTORC1 signaling and better
## outcomes in estrogen receptor-positive breast cancer. Proc Natl
## Acad Sci U S A 2010 Jun 1;107(22):10208-13. PMID: 20479250
##
## Loi S, Sotiriou C, Haibe-Kains B, Lallemand F et al. Gene
## expression profiling identifies activated growth factor signaling
## in poor prognosis (Luminal-B) estrogen receptor positive breast
## cancer. BMC Med Genomics 2009 Jun 24;2:37. PMID: 19552798

gset <- getGEO("GSE6532", destdir = 'external_validation/GSE6532', GSEMatrix =TRUE)

gset <- gset[-1]

gset[[1]]

gset[[1]]$`samplename:ch1`[is.na(gset[[1]]$`samplename:ch1`)] <- gset[[1]]$`sample name:ch1`[is.na(gset[[1]]$`samplename:ch1`)]

ex96 <- exprs(gset[[1]])
ex97 <- exprs(gset[[2]])

clin

colnames(ex96) <- gset[[1]]$`samplename:ch1`
colnames(ex97) <- gset[[2]]$`samplename:ch1`

ex <- rbind(ex96, ex97[,colnames(ex96)])


# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE6532", "/", annotation(gset[[1]]), "+", annotation(gset[[2]]), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex.ave <- avereps(ex, ID=c(gset[[1]]@featureData@data[['Gene Symbol']],
                           gset[[2]]@featureData@data[['Gene Symbol']]))

ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

rownames(ex.ave) <- sapply(strsplit(rownames(ex.ave), ' /// '), function(x) unique(x))

hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex.ave <- ex.ave[hvar.genes,]


title <- paste("GSE6532", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE6532")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

## gene.subset <- rownames(ex.ave)[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) any(x %in% common.genes))]

## ex.ave <- ex.ave[gene.subset,]

clin570 <- gset[[1]]@phenoData@data
rownames(clin570) <- clin570$`samplename:ch1`
clin96 <- gset[[2]]@phenoData@data
samp.ids96 <- ifelse(is.na(clin96$`samplename:ch1`), clin96$`sample name:ch1`, clin96$`samplename:ch1`)
rownames(clin96) <- samp.ids96
clin97 <- gset[[3]]@phenoData@data
rownames(clin97) <- clin97$`samplename:ch1`
## clin97 <- clin97[rownames(clin96),]

clin570 <- clin570 %>% dplyr::select(geo_accession, `er:ch1`, `age:ch1`, `node:ch1`, `size:ch1`, `e.dmfs:ch1`, `t.dmfs:ch1`, `e.rfs:ch1`, `t.rfs:ch1`) %>% mutate_at(c('er:ch1', 'age:ch1', 'node:ch1', 'size:ch1', 'e.dmfs:ch1', 't.dmfs:ch1', 'e.rfs:ch1', 't.rfs:ch1'), as.numeric)

clin96 <- clin96 %>% dplyr::select(geo_accession, `er:ch1`, `age:ch1`, `node:ch1`, `size:ch1`, `e.dmfs:ch1`, `t.dmfs:ch1`, `e.rfs:ch1`, `t.rfs:ch1`) %>% mutate_at(c('er:ch1', 'age:ch1', 'node:ch1', 'size:ch1', 'e.dmfs:ch1', 't.dmfs:ch1', 'e.rfs:ch1', 't.rfs:ch1'), as.numeric)

clin97 <- clin97 %>% dplyr::select(geo_accession, `er:ch1`, `age:ch1`, `node:ch1`, `size:ch1`, `e.dmfs:ch1`, `t.dmfs:ch1`, `e.rfs:ch1`, `t.rfs:ch1`) %>% mutate_at(c('er:ch1', 'age:ch1', 'node:ch1', 'size:ch1', 'e.dmfs:ch1', 't.dmfs:ch1', 'e.rfs:ch1', 't.rfs:ch1'), as.numeric)

intersect(rownames(clin570), rownames(clin96))

samps.shared <- intersect(rownames(clin96), rownames(clin97))

clin96[samps.shared,][is.na(clin96[samps.shared,])] <- clin97[samps.shared,][is.na(clin96[samps.shared,])]

## clin[is.na(clin97)] <- clin96[is.na(clin97)]

rownames(clin96) <- clin96$geo_accession

clin96 <- clin96[,-1]

clin96 <- na.omit(clin96)

colnames(clin96) <- c('ER_STATUS', 'AGE_AT_DIAGNOSIS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'DRFS.status', 'DRFS.time', 'DFS.status', 'DFS.time')

clin96 <- clin96 %>% mutate_all(as.numeric)

clin96$ER_STATUS <- ifelse(grepl(1, clin96$ER_STATUS), 'Positive', 'Negative')
clin96$LYMPH_NODE_STATUS <- ifelse(grepl('1', clin96$LYMPH_NODE_STATUS), 'Positive', 'Negative')
clin96$DFS.status <- as.integer(clin96$DFS.status=='1')
clin96$DRFS.status <- as.integer(clin96$DRFS.status=='1')
clin96$TUMOR_SIZE <- 10 * clin96$TUMOR_SIZE

ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ ER_STATUS, clin96), pval=T)
ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ LYMPH_NODE_STATUS, clin96), pval=T)

ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ ER_STATUS, clin96), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ LYMPH_NODE_STATUS, clin96), pval=T)

clinList[['GSE6532.GPL96']] <- clin96

rownames(clin570) <- clin570$geo_accession

clin570 <- clin570[,-1]

clin570 <- na.omit(clin570)

colnames(clin570) <- c('ER_STATUS', 'AGE_AT_DIAGNOSIS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'DRFS.status', 'DRFS.time', 'DFS.status', 'DFS.time')

clin570 <- clin570 %>% mutate_all(as.numeric)

clin570$ER_STATUS <- ifelse(grepl(1, clin570$ER_STATUS), 'Positive', 'Negative')
clin570$LYMPH_NODE_STATUS <- ifelse(grepl('1', clin570$LYMPH_NODE_STATUS), 'Positive', 'Negative')
clin570$DFS.status <- as.integer(clin570$DFS.status=='1')
clin570$DRFS.status <- as.integer(clin570$DRFS.status=='1')
clin570$TUMOR_SIZE <- 10 * clin570$TUMOR_SIZE

ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ ER_STATUS, clin570), pval=T)
ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ LYMPH_NODE_STATUS, clin570), pval=T)

ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ ER_STATUS, clin570), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ LYMPH_NODE_STATUS, clin570), pval=T)

clinList[['GSE6532.GPL570']] <- clin570


GSE6532.hvar <- cbind(t(ex.ave), clin[colnames(ex.ave),])

dim(GSE6532.hvar)

ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ quantcut(UBE2C, 2), GSE6532.hvar), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ quantcut(UBE2C, 2), GSE6532.hvar), pval=T)

write.csv(GSE6532.hvar, 'external_validation/GSE6532/GSE6532_hvar.csv')


#### GSE7390 Transbig cohort ####
##
##  Loi S, Haibe-Kains B, Desmedt C, Wirapati P et al. Predicting
##  prognosis using molecular profiling in estrogen receptor-positive
##  breast cancer treated with tamoxifen. BMC Genomics 2008 May
##  22;9:239. PMID: 18498629
##
##  Loi S, Haibe-Kains B, Majjaj S, Lallemand F et al. PIK3CA
##  mutations associated with gene signature of low mTORC1 signaling
##  and better outcomes in estrogen receptor-positive breast
##  cancer. Proc Natl Acad Sci U S A 2010 Jun
##  1;107(22):10208-13. PMID: 20479250

gset <- getGEO("GSE7390", destdir = 'external_validation/GSE7390', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE7390", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
## sapply(strsplit(gset@featureData@data[['Gene Symbol']], ' /// ')
ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

rownames(ex.ave) <- sapply(strsplit(rownames(ex.ave), ' /// '), function(x) unique(x))

hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex.ave <- ex.ave[hvar.genes,]


title <- paste("GSE7390", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE7390")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

## gene.subset <- rownames(ex.ave)[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) any(x %in% common.genes))]

## ex.ave <- ex.ave[gene.subset,]

clin <- gset@phenoData@data %>% dplyr::select(`er:ch1`, `age:ch1`, `node:ch1`, `size:ch1`, `e.os:ch1`, `t.os:ch1`, `e.rfs:ch1`, `t.rfs:ch1`, `e.dmfs:ch1`, `t.dmfs:ch1`) %>% mutate_at(c('er:ch1', 'age:ch1', 'node:ch1', 'size:ch1', 'e.os:ch1', 't.os:ch1', 'e.rfs:ch1', 't.rfs:ch1', 'e.dmfs:ch1', 't.dmfs:ch1'), as.numeric)

clin <- na.omit(clin)

colnames(clin) <- c('ER_STATUS', 'AGE_AT_DIAGNOSIS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'OS.status', 'OS.time', 'DFS.status', 'DFS.time', 'DRFS.status', 'DRFS.time')

clin$ER_STATUS <- ifelse(grepl(1, clin$ER_STATUS), 'Positive', 'Negative')
clin$LYMPH_NODE_STATUS <- ifelse(grepl('1', clin$LYMPH_NODE_STATUS), 'Positive', 'Negative')
clin$DFS.status <- as.integer(clin$DFS.status=='1')
clin$OS.status <- as.integer(clin$OS.status=='1')
clin$DRFS.status <- as.integer(clin$DRFS.status=='1')
clin$TUMOR_SIZE <- 10*clin$TUMOR_SIZE

ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

ggsurvplot(survfit(Surv(OS.time, OS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(OS.time, OS.status) ~ quantcut(AGE_AT_DIAGNOSIS,2), clin), pval=T)
ggsurvplot(survfit(Surv(OS.time, OS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

clinList[['GSE7390']] <- clin


GSE7390.hvar <- cbind(t(ex.ave), clin[colnames(ex.ave),])

dim(GSE7390.hvar)

ggsurvplot(survfit(Surv(OS.time, OS.status) ~ quantcut(UBE2C, 2), GSE7390.hvar), pval=T)
ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ quantcut(UBE2C, 2), GSE7390.hvar), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ quantcut(UBE2C, 2), GSE7390.hvar), pval=T)

write.csv(GSE7390.hvar, 'external_validation/GSE7390/GSE7390_hvar.csv')


#### GSE45255 IRB/JNR/NUH cohort ####
##
##  Nagalla, S. et al. Interactions between immunity, proliferation
##  and molecular subtype in breast cancer prognosis. Genome Biol. 14,
##  R34 (2013).

gset <- getGEO("GSE45255", destdir = 'external_validation/GSE45255', GSEMatrix =TRUE)

gset <- gset[[1]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE45255", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
## sapply(strsplit(gset@featureData@data[['Gene Symbol']], ' /// ')
ex.ave <- avereps(ex, ID=gset@featureData@data[['Gene Symbol']])

ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

rownames(ex.ave) <- sapply(strsplit(rownames(ex.ave), ' /// '), function(x) unique(x))

hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex.ave <- ex.ave[hvar.genes,]


title <- paste("GSE45255", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE45255")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

## gene.subset <- rownames(ex.ave)[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) any(x %in% common.genes))]

## ex.ave <- ex.ave[gene.subset,]

clin <- gset@phenoData@data %>% dplyr::select(`er status:ch1`, `patient age:ch1`, `ln status:ch1`, `size (mm):ch1`, `dss event (defined as death from breast cancer):ch1`, `DSS time:ch1`, `dfs event (defined as any type of recurrence or death from breast cancer):ch1`, `dfs time:ch1`, `dmfs event (defined as distant metastasis or death from breast cancer):ch1`, `dmfs time:ch1`) %>% mutate_at(c('patient age:ch1', 'size (mm):ch1', 'dss event (defined as death from breast cancer):ch1', 'DSS time:ch1', 'dfs event (defined as any type of recurrence or death from breast cancer):ch1', 'dfs time:ch1', 'dmfs event (defined as distant metastasis or death from breast cancer):ch1', 'dmfs time:ch1'), as.numeric)

clin[clin=='NA'] <- NA

clin <- clin[rownames(na.omit(clin[,c('er status:ch1', 'patient age:ch1', 'ln status:ch1', 'size (mm):ch1')])),]
clin <- clin[rowSums(is.na(clin[,c('DSS time:ch1','dfs time:ch1','dmfs time:ch1')]))!=3 &
             rowSums(is.na(clin[,c('dss event (defined as death from breast cancer):ch1','dfs event (defined as any type of recurrence or death from breast cancer):ch1','dmfs event (defined as distant metastasis or death from breast cancer):ch1')]))!=3,]

dim(clin)

colnames(clin) <- c('ER_STATUS', 'AGE_AT_DIAGNOSIS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'DSS.status', 'DSS.time', 'DFS.status', 'DFS.time', 'DRFS.status', 'DRFS.time')

clin$ER_STATUS <- ifelse(grepl('[+]', clin$ER_STATUS), 'Positive', 'Negative')
clin$LYMPH_NODE_STATUS <- ifelse(grepl('[+]', clin$LYMPH_NODE_STATUS), 'Positive', 'Negative')
clin$DFS.status <- as.integer(clin$DFS.status=='1')
clin$DSS.status <- as.integer(clin$DSS.status=='1')
clin$DRFS.status <- as.integer(clin$DRFS.status=='1')

clin$DFS.time <- 365.25*clin$DFS.time
clin$DSS.time <- 365.25*clin$DSS.time
clin$DRFS.time <- 365.25*clin$DRFS.time


ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)


clinList[['GSE45255']] <- clin

GSE45255.hvar <- cbind(t(ex.ave), clin[colnames(ex.ave),])

dim(GSE45255.hvar)

ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ quantcut(UBE2C, 2), GSE45255.hvar), pval=T)
ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ quantcut(UBE2C, 2), GSE45255.hvar), pval=T)
ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ quantcut(UBE2C, 2), GSE45255.hvar), pval=T)

write.csv(GSE45255.hvar, 'external_validation/GSE45255/GSE45255_hvar.csv')


#### GSE3494 Upp Cohort 


gset <- getGEO("GSE3494", destdir = 'external_validation/GSE3494', GSEMatrix=TRUE)

gset96 <- gset[[1]]
annotation(gset96)

ex96 <- exprs(gset96)
# log2 transform
qx96 <- as.numeric(quantile(ex96, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx96[5] > 100) || (qx96[6]-qx96[1] > 50 && qx96[2] > 0)
if (LogC) {
    ex96[which(ex96 <= 0)] <- NaN
    ex96 <- log2(ex96)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE3494", "/", annotation(gset96), " value distribution", sep ="")
plotDensities(ex96, main=title, legend=F)

ex96 <- na.omit(ex96) # eliminate rows with NAs
plotSA(lmFit(ex96), main="Mean variance trend, GSE3494")


gset97 <- gset[[2]]
annotation(gset97)

ex97 <- exprs(gset97)
# log2 transform
qx97 <- as.numeric(quantile(ex97, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx97[5] > 100) || (qx97[6]-qx97[1] > 50 && qx97[2] > 0)
LogC
if (LogC) {
    ex97[which(ex97 <= 0)] <- NaN
    ex97 <- log2(ex97)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE3494", "/", annotation(gset97), " value distribution", sep ="")
plotDensities(ex97, main=title, legend=F)

ex97 <- na.omit(ex97) # eliminate rows with NAs
plotSA(lmFit(ex97), main="Mean variance trend, GSE3494")

sampMap <- read.csv('external_validation/GSE3494/samplePatientMap.csv', header=T, sep='\t')
rownames(sampMap) <- sampMap$GEO.Sample.Accession..
head(sampMap)

colnames(ex96) <- sampMap[colnames(ex96), 'Patient.ID']
colnames(ex97) <- sampMap[colnames(ex97), 'Patient.ID']

ex.joint <- rbind(ex96, ex97[,colnames(ex96)])

dim(ex.joint)

dim(ex.joint)

ex.ave <- avereps(ex.joint, ID=c(gset96@featureData@data[['Gene Symbol']],
                                 gset97@featureData@data[['Gene Symbol']]))

ex.ave <- ex.ave[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) length(unique(x))==1),]

rownames(ex.ave) <- sapply(strsplit(rownames(ex.ave), ' /// '), function(x) unique(x))

hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex.ave <- ex.ave[hvar.genes,]


title <- paste("GSE3494", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex.ave, main=title, legend=F)
plotSA(lmFit(ex.ave), main="Mean variance trend, GSE3494")
dim(ex.ave)

head(gset@assayData$exprs)
head(ex.ave)

## gene.subset <- rownames(ex.ave)[sapply(strsplit(rownames(ex.ave), ' /// '), function(x) any(x %in% common.genes))]

## ex.ave <- ex.ave[gene.subset,]

sampMap <- read.csv('external_validation/GSE3494/samplePatientMap.csv', header=T, sep='\t')
sampMap <- sampMap %>% filter(Affy.platform=='HG-U133A')
rownames(sampMap) <- sampMap$Patient.ID

head(sampMap)

clin <- read.csv('external_validation/GSE3494/clinData.csv', header=T, sep='\t', check.names=F)
rownames(clin) <- sampMap[clin$`INDEX (ID)`,'GEO.Sample.Accession..']
head(clin)


clin <- clin %>% dplyr::select(`ER status`, `age at diagnosis`, `Lymph node status`, `tumor size (mm)`, `DSS EVENT (Disease-Specific Survival EVENT; 1=death from breast cancer, 0=alive or censored )`, `DSS TIME (Disease-Specific Survival Time in years)`) %>% mutate_at(c('age at diagnosis', 'tumor size (mm)', 'DSS EVENT (Disease-Specific Survival EVENT; 1=death from breast cancer, 0=alive or censored )', 'DSS TIME (Disease-Specific Survival Time in years)'), as.numeric)

clin <- na.omit(clin)

colnames(clin) <- c('ER_STATUS', 'AGE_AT_DIAGNOSIS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'DSS.status', 'DSS.time')

clin$ER_STATUS <- ifelse(grepl('[-]', clin$ER_STATUS), 'Negative', 'Positive')
clin$LYMPH_NODE_STATUS <- ifelse(grepl('[-]', clin$LYMPH_NODE_STATUS), 'Negative', 'Positive')
clin$DSS.status <- as.integer(clin$DSS.status=='1')
clin$DSS.time <- 365.25*clin$DSS.time

ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

clinList[['GSE3494']] <- clin

GSE3494.hvar <- na.omit(cbind(t(ex.ave), clin[colnames(ex.ave),]))

dim(GSE3494.hvar)

GSE3494.hvar %>% group_by(ER_STATUS) %>% summarize_at('ESR1', mean)


ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ quantcut(UBE2C, 2), GSE3494.hvar), pval=T)

write.csv(GSE3494.hvar, 'external_validation/GSE3494/GSE3494_hvar.csv')



#### GSE96058 SCAN-B Cohort


gset <- getGEO("GSE96058", destdir = 'external_validation/scan-b', GSEMatrix=TRUE)

## gset <- gset[[1]]
## rownames(gset[[1]]@phenoData@data)
## rownames(gset[[2]]@phenoData@data) %in% rownames(gset[[1]]@phenoData@data)

## samps <- c(gset[[1]]@phenoData@data$geo_accession,
##            gset[[2]]@phenoData@data$geo_accession)

clin <- rbind(gset[[1]]@phenoData@data, gset[[2]]@phenoData@data)

rownames(clin) <- clin$title

## clin.GSE96058 <- read.csv('external_validation/scan-b/file3ea7511f82fb.tsv', sep='\t', skip=19)

## head(clin.GSE96058$SAMPLE)
## rownames(clin.GSE96058) <- clin.GSE96058$SAMPLE

clin <- clin %>% dplyr::select(c('er status:ch1', 'age at diagnosis:ch1',
                          'lymph node status:ch1', 'tumor size:ch1',
                          'overall survival event:ch1', 'overall survival days:ch1')) %>%
    mutate_at(c('er status:ch1', 'age at diagnosis:ch1','tumor size:ch1',
                'overall survival event:ch1', 'overall survival days:ch1'), as.numeric)

clin[clin=="NA"] <- NA

colnames(clin) <- c('ER_STATUS', 'AGE_AT_DIAGNOSIS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'OS.status', 'OS.time')

ex <- read.table('external_validation/scan-b/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv', sep=',', row.names=1, header=T, colClasses=c('character', rep('numeric',3409)))[,rownames(clin)]

for (id in unique(gsub('repl','',colnames(ex)))) {
    if (paste0(id, 'repl') %in% colnames(ex)) {
        ex[,id] <- rowMeans(ex[,c(id, paste0(id, 'repl'))])
    }
}

ex <- ex[,!grepl('repl',colnames(ex))]
clin <- clin[!grepl('repl',rownames(clin)),]

## clin <- clin[!is.na(clin$TUMOR_SIZE),]
## clin <- clin[!is.na(clin$ER_STATUS),]

dim(ex)

## ex <- fread('external_validation/scan-b/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv', sep=',')

## ex <- vroom('external_validation/scan-b/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv', delim=',', col_select=rownames(clin.GSE96058))

## ex <- ex[ERsplit.genes,]
head(ex)

## ex <- ex[,colnames(ex) %in% rownames(clin.GSE96058)]
## ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE96058", " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

ex <- na.omit(ex) # eliminate rows with NAs
hvar.genes[!hvar.genes %in% rownames(ex.ave)]

ex <- ex[hvar.genes,]

clin[rownames(clin)[is.na(clin$ER_STATUS)],'ER_STATUS'] <- as.integer(ex['ESR1',rownames(clin)[is.na(clin$ER_STATUS)]]>4)

clin <- na.omit(clin)

clin$ER_STATUS <- ifelse(grepl('0', clin$ER_STATUS), 'Negative', 'Positive')
clin$LYMPH_NODE_STATUS <- gsub('Node','',clin$LYMPH_NODE_STATUS)
clin$OS.status <- as.integer(clin$OS.status=='1')

ggsurvplot(survfit(Surv(OS.time, OS.status) ~ ER_STATUS, clin), pval=T)
ggsurvplot(survfit(Surv(OS.time, OS.status) ~ LYMPH_NODE_STATUS, clin), pval=T)

clinList[['GSE96058']] <- clin


GSE96058.hvar <- na.omit(cbind(t(ex), clin[colnames(ex),]))

dim(GSE96058.hvar)

GSE96058.hvar %>% group_by(ER_STATUS) %>% summarize_at('ESR1', mean)

ggsurvplot(survfit(Surv(OS.time, OS.status) ~ quantcut(UBE2C, 2), GSE96058.hvar), pval=T)

write.csv(GSE96058.hvar, 'external_validation/scan-b/GSE96058_hvar.csv')




##### ComBat Integration
library(sva)


joint.expr <- metabric[common.genes,rownames(clinList[['METABRIC']])]
dim(joint.expr)

clin.covs <- c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE', 'ER_STATUS', 'LYMPH_NODE_STATUS')

joint.clin <- cbind.data.frame(Batch='METABRIC', clinList[['METABRIC']])

dat.sc <- scale(t(joint.expr))

pca <- svd(dat.sc, nu=30, nv=30)

pc.plot <- cbind.data.frame(PC=dat.sc %*% pca$v, joint.clin)


pdf('individual_dataset_pca_plots.pdf', width=6.5, height=5)

ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=ER_STATUS, shape=LYMPH_NODE_STATUS)) +
    geom_point() +
    geom_hline(aes(yintercept=mean(PC.2)-6*sd(PC.2))) +
    geom_hline(aes(yintercept=mean(PC.2)+6*sd(PC.2))) +
    geom_vline(aes(xintercept=mean(PC.1)-6*sd(PC.1))) +
    geom_vline(aes(xintercept=mean(PC.1)+6*sd(PC.1))) + 
    theme_bw() +
    ggtitle('METABRIC')


joint.expr <- cbind.data.frame(joint.expr, scanb[common.genes,rownames(clinList[['GSE96058']])])
joint.clin <- rbind(joint.clin, cbind.data.frame(Batch='GSE96058', clinList[['GSE96058']][,clin.covs]))

dat.sc <- scale(t(joint.expr[,joint.clin$Batch=='GSE96058']))

pca <- svd(dat.sc, nu=30, nv=30)

pc.plot <- cbind.data.frame(PC=dat.sc %*% pca$v, joint.clin[joint.clin$Batch=='GSE96058',])

ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=ER_STATUS, shape=LYMPH_NODE_STATUS)) +
    geom_point() +
    geom_hline(aes(yintercept=mean(PC.2)-6*sd(PC.2))) +
    geom_hline(aes(yintercept=mean(PC.2)+6*sd(PC.2))) +
    geom_vline(aes(xintercept=mean(PC.1)-6*sd(PC.1))) +
    geom_vline(aes(xintercept=mean(PC.1)+6*sd(PC.1))) + 
    theme_bw() +
    ggtitle('GSE96058')

for (gseid in ma.list) {
    if (length(geneData.list[[gseid]])==1) {
        ex <- exprs(geneData.list[[gseid]][[1]][common.genes,rownames(clinList[[gseid]])])

        joint.expr <- cbind.data.frame(joint.expr, ex)
        joint.clin <- rbind(joint.clin, cbind.data.frame(Batch=gseid, clinList[[gseid]][,clin.covs]))

        dat.sc <- scale(t(ex))

        pca <- svd(dat.sc, nu=30, nv=30)

        pc.plot <- cbind.data.frame(PC=dat.sc %*% pca$v, joint.clin[joint.clin$Batch==gseid,])

        print(
            ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=ER_STATUS, shape=LYMPH_NODE_STATUS)) +
            geom_point() +
            geom_hline(aes(yintercept=mean(PC.2)-6*sd(PC.2))) +
            geom_hline(aes(yintercept=mean(PC.2)+6*sd(PC.2))) +
            geom_vline(aes(xintercept=mean(PC.1)-6*sd(PC.1))) +
            geom_vline(aes(xintercept=mean(PC.1)+6*sd(PC.1))) + 
            theme_bw() +
            ggtitle(gseid)
        )
        
    } else {
        for (geneData in geneData.list[[gseid]]) {
            label <- paste(gseid, annotation(geneData), sep='.')
            ex <- exprs(geneData[common.genes,rownames(clinList[[label]])])

            joint.expr <- cbind.data.frame(joint.expr, ex)
            joint.clin <- rbind(joint.clin,
                                cbind.data.frame(Batch=label, clinList[[label]][,clin.covs]))

            dat.sc <- scale(t(ex))

            pca <- svd(dat.sc, nu=30, nv=30)

            pc.plot <- cbind.data.frame(PC=dat.sc %*% pca$v, joint.clin[joint.clin$Batch==label,])

            print(
                ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=ER_STATUS, shape=LYMPH_NODE_STATUS)) +
                geom_point() +
                geom_hline(aes(yintercept=mean(PC.2)-6*sd(PC.2))) +
                geom_hline(aes(yintercept=mean(PC.2)+6*sd(PC.2))) +
                geom_vline(aes(xintercept=mean(PC.1)-6*sd(PC.1))) +
                geom_vline(aes(xintercept=mean(PC.1)+6*sd(PC.1))) + 
                theme_bw() +
                ggtitle(gseid)
            )
            
        }
    }
}

dev.off()

dat.sc <- scale(t(joint.expr))

pca <- svd(dat.sc, nu=30, nv=30)

pc.plot <- cbind.data.frame(PC=dat.sc %*% pca$v, joint.clin)

pdf('all_datasets_pca_uncorrected.pdf', width=5.5, height=5)
ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=Batch)) +
    geom_point() +
    theme_bw()
dev.off()

f <- as.formula(paste('~', paste(clin.covs, collapse=' + ')))
f

combat.expr <- ComBat(joint.expr, joint.clin$Batch,
                      mod=model.matrix(f, joint.clin)[,-1], ref.batch='METABRIC')


dat.sc <- scale(t(combat.expr))

pca <- svd(dat.sc, nu=30, nv=30)

pc.plot <- cbind.data.frame(PC=dat.sc %*% pca$v, joint.clin)


pdf('all_datasets_pca_combat.pdf', width=5.5, height=5)

ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=Batch)) +
    geom_point() +
    scale_color_brewer(palette='Set3') +
    theme_bw()

ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=ER_STATUS, shape=LYMPH_NODE_STATUS)) +
    geom_point() +
    theme_bw()

ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=log10(TUMOR_SIZE+1))) +
    geom_point() +
    theme_bw()

ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=AGE_AT_DIAGNOSIS)) +
    geom_point() +
    theme_bw()

dev.off()


ump <- umap(pc.plot[,1:30], n_neighbors = 50, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)

head(ump$layout)

ump.plot <- cbind.data.frame(UMAP=ump$layout, joint.clin)

pdf('all_datasets_umap_combat.pdf', width=5.5, height=5)

ggplot(ump.plot, aes(x=UMAP.1, y=UMAP.2, color=Batch)) +
    geom_point() +
    scale_color_brewer(palette='Set3') +
    theme_bw()

ggplot(ump.plot, aes(x=UMAP.1, y=UMAP.2, color=ER_STATUS, shape=LYMPH_NODE_STATUS)) +
    geom_point() +
    theme_bw()

ggplot(ump.plot, aes(x=UMAP.1, y=UMAP.2, color=log10(TUMOR_SIZE+1))) +
    geom_point() +
    theme_bw()

ggplot(ump.plot, aes(x=UMAP.1, y=UMAP.2, color=AGE_AT_DIAGNOSIS)) +
    geom_point() +
    theme_bw()

dev.off()



plotDensities(joint.expr)
plotDensities(combat.expr)

pdf('individual_dataset_pca_plots_combat.pdf', width=6.5, height=5)

for (gseid in c('GSE96058', ma.list)) {
    if (gseid=='GSE96058' || length(geneData.list[[gseid]])==1) {
        ex <- combat.expr[,rownames(clinList[[gseid]])]

        full <- cbind.data.frame(t(ex), clinList[[gseid]])

        fname <- paste0('external_validation/',
                        ifelse(gseid=='GSE96058','scan-b',gseid),
                        '/',gseid, '_combat.csv')

        print(fname)
        print(dim(full))

        write.csv(full, fname)

        dat.sc <- scale(t(ex))

        pca <- svd(dat.sc, nu=30, nv=30)

        pc.plot <- cbind.data.frame(PC=dat.sc %*% pca$v, joint.clin[joint.clin$Batch==gseid,])

        print(
            ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=ER_STATUS, shape=LYMPH_NODE_STATUS)) +
            geom_point() +
            geom_hline(aes(yintercept=mean(PC.2)-6*sd(PC.2))) +
            geom_hline(aes(yintercept=mean(PC.2)+6*sd(PC.2))) +
            geom_vline(aes(xintercept=mean(PC.1)-6*sd(PC.1))) +
            geom_vline(aes(xintercept=mean(PC.1)+6*sd(PC.1))) + 
            theme_bw() +
            ggtitle(gseid)
        )
        
    } else {
        full <- data.frame()
        for (geneData in geneData.list[[gseid]]) {
            label <- paste(gseid, annotation(geneData), sep='.')
            
            ex <- combat.expr[,rownames(clinList[[label]])]

            full <- rbind.data.frame(full, cbind.data.frame(t(ex), clinList[[label]]))
        }

        fname <- paste0('external_validation/',
                        ifelse(gseid=='GSE96058','scan-b',gseid),
                        '/',gseid, '_combat.csv')

        print(fname)
        print(dim(full))

        write.csv(full, fname)

        dat.sc <- scale(full[,common.genes])

        pca <- svd(dat.sc, nu=30, nv=30)

        pc.plot <- cbind.data.frame(PC=dat.sc %*% pca$v, full[,setdiff(colnames(full),common.genes)])

        print(
            ggplot(pc.plot, aes(x=PC.1, y=PC.2, color=ER_STATUS, shape=LYMPH_NODE_STATUS)) +
            geom_point() +
            geom_hline(aes(yintercept=mean(PC.2)-6*sd(PC.2))) +
            geom_hline(aes(yintercept=mean(PC.2)+6*sd(PC.2))) +
            geom_vline(aes(xintercept=mean(PC.1)-6*sd(PC.1))) +
            geom_vline(aes(xintercept=mean(PC.1)+6*sd(PC.1))) + 
            theme_bw() +
            ggtitle(gseid)
        )    
    }
}

dev.off()











GSE96058 <- cbind(t(ex), clin.GSE96058)

row.has.repl <- unlist(lapply(rownames(GSE96058), function(x) str_detect(x, pattern='repl')))
GSE96058 <- GSE96058[!row.has.repl,]

table(GSE96058[,c('ER.Status', 'Pam50.Subtype')])
adj.rand.index(GSE96058$ER.Status, GSE96058$Pam50.Subtype)

fit <- Mclust(GSE96058$ESR1, G=2, model='V')
summary(fit)
GSE96058$ER.Status <- as.character(fit$classification-1)

table(GSE96058[,c('ER.Status', 'Pam50.Subtype')])
adj.rand.index(GSE96058$ER.Status, GSE96058$Pam50.Subtype)

tail(GSE96058)
row.has.na <- apply(GSE96058, 1, function(x){any(x=="NA")})
sum(row.has.na)
GSE96058 <- GSE96058[!row.has.na,]
tail(GSE96058)

GSE96058$Age.at.Diagnosis <- as.numeric(GSE96058$Age.at.Diagnosis)
GSE96058$Tumor.Size <- as.numeric(GSE96058$Tumor.Size)

## GSE96058 <- GSE96058 %>% mutate_if(function(x) { is.numeric(x) },
##                                    function(x) { (x-mean(x))/sd(x) })

GSE96058 <- GSE96058 %>% mutate_if(function(x) { is.numeric(x) },
                                   function(x) { huge.npn((matrix(x)-mean(x))/sd(x), npn.func='truncation') })

GSE96058$survival <- as.numeric(GSE96058$survival) / 365.25
GSE96058$event <- as.numeric(GSE96058$event)

GSE96058$ER.Status[GSE96058$ER.Status=='1'] = 'Positive'
GSE96058$ER.Status[GSE96058$ER.Status=='0'] = 'Negative'

## GSE96058$Lymph.Node.Status[GSE96058$Lymph.Node.Status=='LN+'] = 'NodePositive'
## GSE96058$Lymph.Node.Status[GSE96058$Lymph.Node.Status=='LN-'] = 'NodeNegative'

## t.test(GSE96058$ESR1[fit$classification==1], GSE96058$ESR1[fit$classification==2])

t.test(GSE96058$ESR1[GSE96058$ER.Status=='Negative'], GSE96058$ESR1[GSE96058$ER.Status=='Positive'])

## fit <- Mclust(GSE96058$ESR1, G=2, model='V')
## summary(fit)
## GSE96058$ER.Status <- fit$classification-1

ggsurvplot(survfit(Surv(survival, event) ~ (Tumor.Size > 0), GSE96058), pval=T)

ggsurvplot(survfit(Surv(survival, event) ~ (Age.at.Diagnosis > 0), GSE96058), pval=T)

ggsurvplot(survfit(Surv(survival, event) ~ ER.Status, GSE96058), pval=T)

ggsurvplot(survfit(Surv(survival, event) ~ Pam50.Subtype, GSE96058), pval=T)

write.table(GSE96058, 'external_validation/scan-b/data/ERsplit_GSE96058_survival.csv', sep=',', row.names=F)

