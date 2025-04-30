library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(gtools)

hubSummarize <- function(corMat, single.cutoff, average.cutoff) {
    groups <- list()
    hc.single <- hclust(as.dist(1-abs(corMat)), method='single')

    clusts <- cutree(hc.single, h=single.cutoff)
    nclusts <- length(unique(clusts))

    keep.genes <- c()

    for (idx in seq_len(nclusts)) {
        gene.group <- names(clusts[clusts==idx])
        if (length(gene.group)==1) {
            keep.genes <- c(keep.genes, gene.group)
        } else {
            hc.ave <- hclust(as.dist(1-abs(corMat[gene.group, gene.group])), method='average')

            clusts.group <- cutree(hc.ave, h=average.cutoff)
            nclusts.group <- length(unique(clusts.group))

            print(idx)
            
            for (jdx in seq_len(nclusts.group)) {
                gene.group.group <- names(clusts.group[clusts.group==jdx])
                if (length(gene.group.group) > 1) {
                    group.corMat <- corMat[gene.group.group,gene.group.group]
                    ave.cors <- (colSums(abs(group.corMat))-1) / (length(gene.group.group)-1)
                    if (max(group.corMat[upper.tri(group.corMat)]) < (1-single.cutoff)) {
                        keep.genes <- c(keep.genes, gene.group.group)
                    } else {
                        max.cor.gene <- names(which.max(ave.cors))
                        print(ave.cors)
                        print(max.cor.gene)
                        keep.genes <- c(keep.genes, max.cor.gene)
                        groups[[max.cor.gene]] <- gene.group.group
                    }
                } else {
                    keep.genes <- c(keep.genes, gene.group.group)
                }
            }
        }
    }

    return(list(keep=keep.genes, groups=groups))
}


expr <- read.table('data_expression_median.txt', sep='\t', row.names=1, header=T)

common.genes <- readRDS('meta_cohort_common_genes.rds')

expr <- as.matrix(expr[common.genes,-1])

dim(expr)

expr <- na.omit(expr)

dim(expr)

expr.mean <- apply(expr, 1, mean)
expr.sd <- apply(expr, 1, sd)

gene.plot <- data.frame(Mean = expr.mean, SD = expr.sd,
                        HighVar = (expr.sd >= sort(expr.sd, decreasing=T)[500]))

ggplot(gene.plot, aes(Mean, SD, color=HighVar)) +
    geom_point() +
    scale_color_manual(values=c('black','red')) + 
    theme_bw()


hvar.500 <- names(sort(expr.sd, decreasing=T)[1:500])

hvar.500

corMat <- cor(t(expr[hvar.500,]))

pheatmap::pheatmap(corMat, breaks=seq(-1,1,length.out=101))


ego <- enrichGO(hvar.500,
                universe      = rownames(expr),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

dim(ego)
head(ego, 50)

dotplot(ego, showCategory=50)

hvar.genes <- hvar.500

cor.mat <- cor(t(expr[hvar.genes,]), method='spearman')

pheatmap::pheatmap(cor.mat, clustering_method='ward.D2')

max(cor.mat[upper.tri(cor.mat)])
min(cor.mat[upper.tri(cor.mat)])
quantile((cor.mat[upper.tri(cor.mat)]), probs=seq(0,1,0.01))

## Remove sets of highly correlated genes, replace each set with the single gene that is most higly correlated to all other genes in the set
rna.hubs <- hubSummarize(cor.mat, 0.2, 0.25)

length(rna.hubs$keep)
length(rna.hubs$groups)

rna.hubs$groups

ego.list <- lapply(rna.hubs$groups,
                   enrichGO,
                   universe      = hvar.500,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

lapply(ego.list, head)

rna.hubs$groups

filt.rna.cor.mat <- cor.mat[rna.hubs$keep,rna.hubs$keep]

pheatmap::pheatmap(filt.rna.cor.mat)
pheatmap::pheatmap(filt.rna.cor.mat[names(rna.hubs$groups),names(rna.hubs$groups)])
max(filt.rna.cor.mat[upper.tri(filt.rna.cor.mat)])
min(filt.rna.cor.mat[upper.tri(filt.rna.cor.mat)])

hvar.expr <- expr[rna.hubs$keep,]

dim(hvar.expr)

for (gene in names(rna.hubs$groups)) {
    rownames(hvar.expr)[which(rownames(hvar.expr)==gene)] <- paste0(gene, '_grp')
}

dim(hvar.expr)

kappa(hvar.expr)

saveRDS(rna.hubs$groups, 'metabric_rna_highCor_gene_groups.RDS')
saveRDS(ego.list, 'metabric_rna_highCor_gene_groups_ora.RDS')

write.csv(hvar.expr, 'data/hvar_genes_500_hub_summarized.csv')

#### Survival data

library(survival)

df <- read.delim('data_clinical_patient_outcomes.txt', sep='\t', row.names=1, header=T)

rownames(df) <- make.names(rownames(df))

dim(df)

df <- df[!(!is.na(df$Histological.Type) & df$Histological.Type=='BENIGN'),]

df[df$TLR==0,'TLR'] <- 1

df[is.na(df$Death),c('DeathBreast','Death')] <- 0

dim(df)

min(df$T)
min(df$TLR)
min(df$TDR)

df$DSS <- Surv(df$T, df$DeathBreast)
df$OD <- Surv(df$T, as.integer(df$Death==1 & df$DeathBreast==0))
## df$AllCause <- Surv(df$T, df$Death)
df$LocalRelapse <- Surv(df$TLR, df$LR)
df$DistantRelapse <- Surv(df$TDR, df$DR)
EOS <- pmax(df$DSS[,2], df$OD[,2])
TOS <- df$DSS[,1]
TOS[df$OD[,2]==1] <- ifelse(df$DSS[df$OD[,2]==1,2]==1,
                            pmin(TOS[df$OD[,2]==1], df$OD[df$OD[,2]==1,1]),
                            df$OD[df$OD[,2]==1,1])
TOS[EOS==0] <- pmax(df$DSS[EOS==0,1], df$OD[EOS==0,1])
df$OS <- Surv(TOS, EOS)
EDRFS <- pmax(df$DSS[,2], df$OD[,2], df$DistantRelapse[,2])
TDRFS <- TOS
TDRFS[df$DistantRelapse[,2]==1] <-
    ifelse(EOS[df$DistantRelapse[,2]==1]==1,
           pmin(TDRFS[df$DistantRelapse[,2]==1],
                df$DistantRelapse[df$DistantRelapse[,2]==1,1]),
           df$DistantRelapse[df$DistantRelapse[,2]==1,1])
TDRFS[EDRFS==0] <- pmax(df$DSS[EDRFS==0,1], df$OD[EDRFS==0,1], df$DistantRelapse[EDRFS==0,1])
df$DRFS <- Surv(TDRFS, EDRFS)
EDFS <- pmax(df$DSS[,2], df$OD[,2], df$DistantRelapse[,2], df$LocalRelapse[,2])
TDFS <- TDRFS
TDFS[df$LocalRelapse[,2]==1] <-
    ifelse(EDRFS[df$LocalRelapse[,2]==1]==1,
           pmin(TDFS[df$LocalRelapse[,2]==1],
                df$LocalRelapse[df$LocalRelapse[,2]==1,1]),
           df$LocalRelapse[df$LocalRelapse[,2]==1,1])
TDFS[EDFS==0] <- pmax(df$DSS[EDFS==0,1], df$OD[EDFS==0,1], df$DistantRelapse[EDFS==0,1], df$LocalRelapse[EDFS==0,1])
df$DFS <- Surv(TDFS, EDFS)

colSums(is.na(df))

dim(df[colnames(hvar.expr),])

## Sample clinical data
clin.samp <- read.delim('data_clinical_sample.txt', sep='\t', row.names=1, header=T, skip=4, na.strings=c('NA', ''))

dim(clin.samp)

rownames(clin.samp) <- make.names(rownames(clin.samp))

head(clin.samp)

table(clin.samp$SAMPLE_TYPE)

table(clin.samp[colnames(hvar.expr),'ER_STATUS'])
table(clin.samp[colnames(hvar.expr),'PR_STATUS'])
table(clin.samp[colnames(hvar.expr),'HER2_STATUS'])
table(clin.samp[colnames(hvar.expr),'GRADE'])
table(clin.samp[colnames(hvar.expr),'CANCER_TYPE_DETAILED'])
table(clin.samp[colnames(hvar.expr),'ONCOTREE_CODE'])

colSums(is.na(clin.samp[colnames(hvar.expr),]))

clin.samp[is.na(clin.samp$ONCOTREE_CODE),'ONCOTREE_CODE'] <- 'BREAST'

samp.feats <- c('ER_STATUS', 'PR_STATUS', 'HER2_STATUS', 'GRADE', 'ONCOTREE_CODE', 'TUMOR_SIZE')
dim(na.omit(clin.samp[colnames(hvar.expr),samp.feats]))


## Patient clinical data
clin.pat <- read.delim('data_clinical_patient.txt', sep='\t', row.names=1, header=T, skip=4, na.strings=c('NA', ''))

dim(clin.pat)

rownames(clin.pat) <- make.names(rownames(clin.pat))

head(clin.pat)

table(clin.pat[colnames(hvar.expr),'ER_IHC'])
table(clin.pat[colnames(hvar.expr),'LYMPH_NODES_EXAMINED_POSITIVE'])
table(clin.pat[colnames(hvar.expr),'HER2_SNP6'])
table(clin.pat[colnames(hvar.expr),'CELLULARITY'])
table(clin.pat[colnames(hvar.expr),'INFERRED_MENOPAUSAL_STATE'])
table(clin.pat[colnames(hvar.expr),'HISTOLOGICAL_SUBTYPE'])

table(clin.pat[colnames(hvar.expr),'HISTOLOGICAL_SUBTYPE'], clin.samp[colnames(hvar.expr),'ONCOTREE_CODE'])

table(data.frame(patient=clin.pat[colnames(hvar.expr),'ER_IHC'], sample=clin.samp[colnames(hvar.expr),'ER_STATUS']))

table(data.frame(patient=clin.pat[colnames(hvar.expr),'HER2_SNP6'], sample=clin.samp[colnames(hvar.expr),'HER2_STATUS']))

table(clin.pat[colnames(hvar.expr),'CELLULARITY'], df[colnames(hvar.expr),'Cellularity'])

colSums(is.na(cbind(clin.pat[colnames(hvar.expr),'CELLULARITY'], df[colnames(hvar.expr),'Cellularity'])))

colSums(is.na(clin.pat[colnames(hvar.expr),]))

pat.feats <- c('LYMPH_NODES_EXAMINED_POSITIVE', 'INFERRED_MENOPAUSAL_STATE', 'AGE_AT_DIAGNOSIS')

dim(na.omit(clin.pat[colnames(hvar.expr),pat.feats]))

dim(na.omit(cbind(clin.pat[colnames(hvar.expr),pat.feats], clin.samp[colnames(hvar.expr),samp.feats])))

clin <- na.omit(cbind(clin.pat[colnames(hvar.expr),pat.feats], clin.samp[colnames(hvar.expr),samp.feats]))

dim(clin)

table(clin$ER_STATUS)

table(clin$HER2_STATUS)

clin$LYMPH_NODE_STATUS <- factor(ifelse(clin$LYMPH_NODES_EXAMINED_POSITIVE==0, 'NodeNegative', ifelse(clin$LYMPH_NODES_EXAMINED_POSITIVE<3, '1to3', '4toX')), levels=c('NodeNegative', '1to3', '4toX'))

clin$GRADE <- paste0('G', clin$GRADE)

## clin$HER2_SNP6[clin$HER2_SNP6=='UNDEF'] <- 'NEUTRAL'

clin$ONCOTREE_CODE[clin$ONCOTREE_CODE=='IMMC'] <- 'BREAST'

table(clin$LYMPH_NODE_STATUS)

table(clin$GRADE)

table(clin$HER2_SNP6)

table(clin$ONCOTREE_CODE)

head(clin)

final.clin <- na.omit(cbind(df[rownames(clin),
                               c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse',
                                 'OS', 'DRFS', 'DFS')],
                            clin[,-1]))

head(final.clin)

dim(final.clin)

samp.ids <- rownames(final.clin)

hvar.expr.df <- data.frame(t(hvar.expr))

head(hvar.expr.df)

metabric <- na.omit(cbind(final.clin[samp.ids,],
                          hvar.expr.df[samp.ids,]))

set.seed(20230215)

foldid <- sample(rep(1:10, length.out=nrow(metabric)))

write.csv(metabric %>% mutate_if(is.Surv, as.matrix), 'data/metabric.rna.full.csv')

for (k in 1:10) {
    write.csv(metabric[foldid!=k,] %>% mutate_if(is.Surv, as.matrix),
              paste0('data/metabric.rna.train.cv', k, '.csv'))

    write.csv(metabric[foldid==k,] %>% mutate_if(is.Surv, as.matrix),
              paste0('data/metabric.rna.test.cv', k, '.csv'))

    print(dim(metabric[foldid!=k,]))
    print(table(metabric[foldid!=k,'ER_STATUS']))
}

metabric.rna.erp <- metabric[metabric$ER_STATUS=='Positive',] %>% dplyr::select(-ER_STATUS)

dim(metabric.rna.erp)

write.csv(metabric.rna.erp %>% mutate_if(is.Surv, as.matrix), 'data/metabric.rna.erp.full.csv')

for (k in 1:10) {
    write.csv(metabric.rna.erp[rownames(metabric)[foldid!=k & metabric$ER_STATUS=='Positive'],] %>% mutate_if(is.Surv, as.matrix),
              paste0('data/metabric.rna.erp.train.cv', k, '.csv'))

    write.csv(metabric.rna.erp[rownames(metabric)[foldid==k & metabric$ER_STATUS=='Positive'],] %>% mutate_if(is.Surv, as.matrix),
              paste0('data/metabric.rna.erp.test.cv', k, '.csv'))
}


metabric.rna.ern <- metabric[metabric$ER_STATUS=='Negative',] %>% dplyr::select(-ER_STATUS)

dim(metabric.rna.ern)

write.csv(metabric.rna.ern %>% mutate_if(is.Surv, as.matrix), 'data/metabric.rna.ern.full.csv')

for (k in 1:10) {
    write.csv(metabric.rna.ern[rownames(metabric)[foldid!=k & metabric$ER_STATUS=='Negative'],] %>% mutate_if(is.Surv, as.matrix),
              paste0('data/metabric.rna.ern.train.cv', k, '.csv'))
   
    write.csv(metabric.rna.ern[rownames(metabric)[foldid==k & metabric$ER_STATUS=='Negative'],] %>% mutate_if(is.Surv, as.matrix),
              paste0('data/metabric.rna.ern.test.cv', k, '.csv'))
    
    print(dim(metabric.rna.ern[rownames(metabric)[foldid!=k & metabric$ER_STATUS=='Negative'],]))
}

