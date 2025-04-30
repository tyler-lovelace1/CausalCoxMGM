library(rCausalMGM)
library(tidyverse)
library(survival)
library(survminer)
library(glmnet)
library(clinfun)

loadFullData <- function(erstatus='all') {

    cvString <- 'full'
  
    if (erstatus=='all') {
        train <- read.csv(paste0('data/metabric.rna.', cvString, '.csv'), row.names=1)
    } else {
        train <- read.csv(paste0('data/metabric.rna.', erstatus, '.', cvString, '.csv'),
                          row.names=1)
    }

    if (any(c('DSS.time', 'DSS.status') %in% colnames(train))) {
        train$DSS <- Surv(train$DSS.time, train$DSS.status)
        if (cvString!='full') {
            test$DSS <- Surv(test$DSS.time, test$DSS.status)
        }
    }
    if (any(c('OD.time', 'OD.status') %in% colnames(train))) {
        train$OD <- Surv(train$OD.time, train$OD.status)
    }
    if (any(c('DistantRelapse.time', 'DistantRelapse.status') %in% colnames(train))) {
        train$DistantRelapse <- Surv(train$DistantRelapse.time, train$DistantRelapse.status)
    }
    if (any(c('LocalRelapse.time', 'LocalRelapse.status') %in% colnames(train))) {
        train$LocalRelapse <- Surv(train$LocalRelapse.time, train$LocalRelapse.status)
    }
    if (any(c('OS.time', 'OS.status') %in% colnames(train))) {
        train$OS <- Surv(train$OS.time, train$OS.status)
    }
    if (any(c('DRFS.time', 'DRFS.status') %in% colnames(train))) {
        train$DRFS <- Surv(train$DRFS.time, train$DRFS.status)
    }
    if (any(c('DFS.time', 'DFS.status') %in% colnames(train))) {
        train$DFS <- Surv(train$DFS.time, train$DFS.status)
    }
    if (any(c('BONE.time', 'BONE.status') %in% colnames(train))) {
        train$BONE <- Surv(train$BONE.time, train$BONE.status)
    }
    if (any(c('BRAIN.time', 'BRAIN.status') %in% colnames(train))) {
        train$BRAIN <- Surv(train$BRAIN.time, train$BRAIN.status)
    }
    if (any(c('LIVER.time', 'LIVER.status') %in% colnames(train))) {
        train$LIVER <- Surv(train$LIVER.time, train$LIVER.status)
    }
    if (any(c('LNS.time', 'LNS.status') %in% colnames(train))) {
        train$LNS <- Surv(train$LNS.time, train$LNS.status)
    }
    if (any(c('LOCOREGIONAL.time', 'LOCOREGIONAL.status') %in% colnames(train))) {
        train$LOCOREGIONAL <- Surv(train$LOCOREGIONAL.time, train$LOCOREGIONAL.status)
    }
    if (any(c('LUNG.time', 'LUNG.status') %in% colnames(train))) {
        train$LUNG <- Surv(train$LUNG.time, train$LUNG.status)
    }
    if (any(c('PLEURA.time', 'PLEURA.status') %in% colnames(train))) {
        train$PLEURA <- Surv(train$PLEURA.time, train$PLEURA.status)
    }

    if (erstatus=='all') {
        train$DSS <- stratifySurv(train$DSS, train$ER_STATUS)
        ## train$OD <- stratifySurv(train$OD, train$ER_STATUS)
        train$DistantRelapse <- stratifySurv(train$DistantRelapse, train$ER_STATUS)
        train$LocalRelapse <- stratifySurv(train$LocalRelapse, train$ER_STATUS)
        train$OS <- stratifySurv(train$OS, train$ER_STATUS)
        train$DRFS <- stratifySurv(train$DRFS, train$ER_STATUS)
        train$DFS <- stratifySurv(train$DFS, train$ER_STATUS)
    }

    train <- train[,!colnames(train) %in% c('DSS.time', 'DSS.status',
                                            'OD.time', 'OD.status',
                                            'DistantRelapse.time', 'DistantRelapse.status',
                                            'LocalRelapse.time', 'LocalRelapse.status',
                                            'OS.time', 'OS.status',
                                            'DRFS.time', 'DRFS.status',
                                            'DFS.time', 'DFS.status',
                                            "BONE.time","BONE.status",
                                            "BRAIN.time","BRAIN.status",
                                            "LIVER.time","LIVER.status",
                                            "LNS.time","LNS.status",
                                            "LOCOREGIONAL.time","LOCOREGIONAL.status",
                                            "LUNG.time","LUNG.status",
                                            "PLEURA.time","PLEURA.status")]
    if (erstatus=='ern') {
        train[train$GRADE!='G3','GRADE'] <- 'G1.G2'
        train[!train$ONCOTREE_CODE!='IDC','ONCOTREE_CODE'] <- 'BREAST'
    }

    train$LYMPH_NODE_STATUS <- ifelse(train$LYMPH_NODE_STATUS=='NodeNegative',
                                      'Negative', 'Positive')

    train <- train %>% mutate_if(is.character, factor)

    numeric.mask <- sapply(train, function(x) !is.Surv(x) & is.numeric(x))
    ## numeric.mask[c('AGE_AT_DIAGNOSIS')] <- FALSE

    train.npn <- train
    train.npn[,numeric.mask] <- huge::huge.npn(as.matrix(train[,numeric.mask]))
    ## train.npn[,'TUMOR_SIZE'] <- log10(train[,'TUMOR_SIZE'])

    return(list(train=train, train.npn=train.npn))
}


loadFullGraph <- function(erstatus='all') {
    if (erstatus=='all') {
        g <- loadGraph('out/fcimaxmgmStARS.metabric.rna.split.all.rank.FDR05.full.txt')
    } else if (erstatus=='erp') {
        g <- loadGraph('out/fcimaxmgmStARS.metabric.rna.split.erp.rank.FDR1.full.txt')
    } else if (erstatus=='ern') {
        g <- loadGraph('out/fcimaxmgmStARS.metabric.rna.split.ern.rank.FDR1.full.txt')
    }
    return(g)
}


trainModels <- function(train, graph, erstatus='all') {
    if (erstatus!='relapse') {
        surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')
    } else {
        surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'OS', 'DRFS', 'DFS', 'BONE',
                        'BRAIN', 'LIVER', 'LNS', 'LOCOREGIONAL', 'LUNG', 'PLEURA')
    }

    models <- list()
    for (sf in surv.feats) {
        if (sf %in% c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse')) {
            mb <- sort(setdiff(graph$markov.blankets[[sf]], c('ER_STATUS', surv.feats)))
            if (erstatus=='all' && sf != 'OD') {
                f.mb <- as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                         paste(c(1, mb), collapse=' + ')))
            } else {
                f.mb <- as.formula(paste(sf, '~ ',
                                         paste(c(1, mb), collapse=' + ')))
            }
        } else {
            sf.splits <- surv.feats[1:(which(surv.feats==sf)-3)]
            mb <- paste('risk',sf.splits,sep='.')
            if (erstatus=='all' && sf != 'OD') {
                f.mb <- as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                         paste(c(1, mb), collapse=' + ')))
            } else {
                f.mb <- as.formula(paste(sf, '~ ',
                                         paste(c(1, mb), collapse=' + ')))
            }
        }

        models[[sf]] <- coxph(f.mb, train)
        train[,paste('risk',sf,sep='.')] <- predict(models[[sf]], newdata=train)
    }
    return(models)
}

surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')

fullData.list <- lapply(list(all='all', erp='erp', ern='ern'), loadFullData)
fullGraph.list <- lapply(list(all='all', erp='erp', ern='ern'), loadFullGraph)



fullModels.list <- lapply(list(all='all', erp='erp', ern='ern'),
                          function(x) trainModels(fullData.list[[x]]$train.npn,
                                                  fullGraph.list[[x]], x))

erpModels.list <-  lapply(list(erp='erp', ern='ern'),
                          function(x) trainModels(fullData.list[['erp']]$train.npn,
                                                  fullGraph.list[[x]], x))

ernModels.list <-  lapply(list(erp='erp', ern='ern'),
                          function(x) trainModels(fullData.list[['ern']]$train.npn,
                                                  fullGraph.list[[x]], x))

erpModels.list

coef.confints <- data.frame()
for (erstatus in c('erp', 'ern')) {
    for (sf in c('DSS', 'DistantRelapse', 'LocalRelapse')) {
        erpMod.sum <- summary(erpModels.list[[erstatus]][[sf]])
        erpMod.ci <- erpMod.sum$conf.int

        ernMod.sum <- summary(ernModels.list[[erstatus]][[sf]])
        ernMod.ci <- ernMod.sum$conf.int

        coef.confints <- rbind.data.frame(
            coef.confints,
            data.frame(Dataset=c(rep('ER+', nrow(erpMod.ci)), rep('ER-', nrow(ernMod.ci))),
                       Model=rep(ifelse(erstatus=='erp', 'ER+', 'ER-'),
                                 nrow(erpMod.ci)+nrow(ernMod.ci)),
                       Outcome=rep(sf, nrow(erpMod.ci)+nrow(ernMod.ci)),
                       Feature=c(rownames(erpMod.ci), rownames(ernMod.ci)),
                       HR=c(erpMod.ci[,1], ernMod.ci[,1]),
                       HR.low=c(erpMod.ci[,3], ernMod.ci[,3]),
                       HR.upper=c(erpMod.ci[,4], ernMod.ci[,4])))
                   
    }
}

## for (idx in 1:nrow(coef.confints)) {
##     x <- as.vector(coef.confints[idx,])

##     mask <- x[2] != coef.confints[,2] & x[3]==coef.confints[,3] & x[4]==coef.confints[,4]

##     coef.confints[mask,'Model'] <- 'Both'
## }

coef.confints$Feature <- factor(coef.confints$Feature,
                                levels=rev(c('LYMPH_NODE_STATUSPositive', 'TUMOR_SIZE',
                                             'UBE2C_grp', 'FGFR4', 'POSTN', 'SERPINA1',
                                             'SERPINA6', 'MAPT', 'ELOVL5', 'S100P', 'SYT17',
                                             'VCAM1', 'C1R')),
                                labels=rev(c('LNS:Positive', 'Tumor Size', 'UBE2C', 'FGFR4',
                                             'POSTN', 'SERPINA1', 'SERPINA6', 'MAPT',
                                             'ELOVL5', 'S100P', 'SYT17', 'VCAM1', 'C1R')))

coef.confints$Outcome <- factor(coef.confints$Outcome,
                                levels=c('DSS', 'DistantRelapse', 'LocalRelapse'),
                                labels=c('DSS', 'Distant', 'Local'))

coef.confints$Model <- factor(coef.confints$Model, levels=c('ER+', 'ER-'))
coef.confints$Dataset <- factor(coef.confints$Dataset, levels=c('ER+', 'ER-'))
coef.confints$Tumor <- factor(coef.confints$Dataset, levels=c('ER+', 'ER-'))

library(ggh4x)

ggplot(coef.confints,
       aes(y=Feature, x=HR, xmin=HR.low, xmax=HR.upper, color=Model)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=2) +
    scale_x_log10() +
    scale_color_brewer(palette='Set1') +
    ## coord_cartesian(xlim =c(0.5, 3.5)) +
    facet_grid(cols=vars(Tumor), rows=vars(Outcome),
               scales='free_y', labeller=label_both) +
    geom_vline(xintercept=1, color='darkgray', lty=5) +
    xlab('Standardized Hazard Ratio') +
    force_panelsizes(rows=c(7, 8, 4)) +
    theme_bw()

ggsave('causal_dss_dr_lr_forestplot_model_compare.png', width=5.5, height=6, dpi=600)


ggplot(coef.confints %>% filter(Outcome != "Local"),
       aes(y=Feature, x=HR, xmin=HR.low, xmax=HR.upper, color=Model)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=2) +
    scale_x_log10() +
    scale_color_brewer(palette='Set1') +
    ## coord_cartesian(xlim =c(0.5, 3.5)) +
    facet_grid(cols=vars(Tumor), rows=vars(Outcome),
               scales='free_y', labeller=label_both) +
    geom_vline(xintercept=1, color='darkgray', lty=5) +
    xlab('Standardized Hazard Ratio') +
    force_panelsizes(rows=c(7, 8)) +
    theme_bw()

ggsave('causal_dss_dr_forestplot_model_compare.png', width=5.5, height=6, dpi=600)


ggplot(coef.confints %>% filter(Tumor=='ER+'),
       aes(y=Feature, x=HR, xmin=HR.low, xmax=HR.upper, color=Model)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=2) +
    scale_x_log10() +
    scale_color_brewer(palette='Set1') +
    ## coord_cartesian(xlim =c(0.5, 3.5)) +
    facet_grid(cols=vars(Tumor), rows=vars(Outcome),
               scales='free_y', labeller=label_both) +
    geom_vline(xintercept=1, color='darkgray', lty=5) +
    xlab('Standardized Hazard Ratio') +
    force_panelsizes(rows=c(7, 8, 4)) +
    theme_bw()

ggsave('causal_dss_dr_lr_forestplot_model_compare_erp.png', width=4, height=6, dpi=600)


ggplot(coef.confints %>% filter(Tumor=='ER-'),
       aes(y=Feature, x=HR, xmin=HR.low, xmax=HR.upper, color=Model)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=2) +
    scale_x_log10() +
    scale_color_brewer(palette='Set1') +
    ## coord_cartesian(xlim =c(0.5, 3.5)) +
    facet_grid(cols=vars(Tumor), rows=vars(Outcome),
               scales='free_y', labeller=label_both) +
    geom_vline(xintercept=1, color='darkgray', lty=5) +
    xlab('Standardized Hazard Ratio') +
    force_panelsizes(rows=c(7, 8, 4)) +
    theme_bw()

ggsave('causal_dss_dr_lr_forestplot_model_compare_ern.png', width=4, height=6, dpi=600)




cph.dss.erp <- coxph(DSS ~ LYMPH_NODE_STATUS + TUMOR_SIZE + UBE2C_grp +
                         FGFR4 + POSTN + SERPINA1 + SERPINA6, fullData.list$erp$train.npn)

cph.dss.ern <- coxph(DSS ~ LYMPH_NODE_STATUS + TUMOR_SIZE + UBE2C_grp +
                         FGFR4 + POSTN + SERPINA1 + SERPINA6, fullData.list$ern$train.npn)

cph.dr.erp <- coxph(DSS ~ LYMPH_NODE_STATUS + TUMOR_SIZE + UBE2C_grp + MAPT + ELOVL5 +
                        S100P + SYT17 + VCAM1, fullData.list$erp$train.npn)

cph.dr.ern <- coxph(DSS ~ LYMPH_NODE_STATUS + TUMOR_SIZE + UBE2C_grp + MAPT + ELOVL5 +
                        S100P + SYT17 + VCAM1, fullData.list$ern$train.npn)

forestplot.dss <- rbind(data.frame(Dataset='ER+', Outcome='DSS',
                                   Model=c('Both', 'Both', 'ER+', rep('ER-', 4)),
                                   Feature=c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp',
                                             'FGFR4', 'POSTN', 'SERPINA1', 'SERPINA6'),
                                   Coef=coef(cph.dss.erp),
                                   SE=sqrt(diag(vcov(cph.dss.erp)))),
                        data.frame(Dataset='ER-', Outcome='DSS',
                                   Model=c('Both', 'Both', 'ER+', rep('ER-', 4)),
                                   Feature=c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp',
                                             'FGFR4', 'POSTN', 'SERPINA1', 'SERPINA6'),
                                   Coef=coef(cph.dss.ern),
                                   SE=sqrt(diag(vcov(cph.dss.ern)))))

forestplot.dr <- rbind(data.frame(Dataset='ER+', Outcome='Distant Relapse',
                                  Model=c('Both', rep('ER+', 4), rep('ER-', 3)),
                                  Feature=c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp',
                                            'MAPT', 'ELOVL5', 'S100P', 'SYT17', 'VCAM1'),
                                  Coef=coef(cph.dr.erp),
                                  SE=sqrt(diag(vcov(cph.dr.erp)))),
                       data.frame(Dataset='ER-', Outcome='Distant Relapse',
                                  Model=c('Both', rep('ER+', 4), rep('ER-', 3)),
                                  Feature=c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp',
                                            'MAPT', 'ELOVL5', 'S100P', 'SYT17', 'VCAM1'),
                                  Coef=coef(cph.dr.ern),
                                  SE=sqrt(diag(vcov(cph.dr.ern)))))

forestplot.dss$Feature <- factor(forestplot.dss$Feature,
                                 levels=rev(c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp',
                                              'FGFR4', 'POSTN', 'SERPINA1', 'SERPINA6')))

forestplot.dss$Dataset <- factor(forestplot.dss$Dataset, levels=c('ER+','ER-'))

forestplot.dss$Model <- factor(forestplot.dss$Model, levels=c('Both', 'ER+','ER-'))

fp.dss <- ggplot(forestplot.dss,
                 aes(x=exp(Coef), y=Feature,
                     xmin=exp(Coef-1.96*SE), xmax=exp(Coef+1.96*SE),
                     color=Model)) +
    geom_point() +
    geom_errorbar(width=0.5) +
    scale_x_log10() +
    scale_color_brewer(palette='Set1') +
    facet_grid(rows=vars(Outcome), cols=vars(Dataset), scales='free', labeller=label_both) +
    geom_vline(xintercept=1, color='darkgray', lty=5) +
    xlab('Standardized Hazard Ratio') +
    theme_bw()

fp.dss

forestplot.dr$Feature <- factor(forestplot.dr$Feature,
                                 levels=rev(c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp',
                                            'MAPT', 'ELOVL5', 'S100P', 'SYT17', 'VCAM1')))

forestplot.dr$Dataset <- factor(forestplot.dr$Dataset, levels=c('ER+','ER-'))

forestplot.dr$Model <- factor(forestplot.dr$Model, levels=c('Both', 'ER+','ER-'))

fp.dr <- ggplot(forestplot.dr,
                 aes(x=exp(Coef), y=Feature,
                     xmin=exp(Coef-1.96*SE), xmax=exp(Coef+1.96*SE),
                     color=Model)) +
    geom_point() +
    geom_errorbar(width=0.5) +
    scale_x_log10() +
    scale_color_brewer(palette='Set1') +
    facet_grid(rows=vars(Outcome), cols=vars(Dataset), scales='free', labeller=label_both) +
    geom_vline(xintercept=1, color='darkgray', lty=5) +
    xlab('Standardized Hazard Ratio') +
    theme_bw()

fp.dr

library(cowplot)
fp <- plot_grid(fp.dss, fp.dr, ncol=1)
fp
ggsave('causal_dss_dr_forestplot_ersplit.png', fp, width=8, height=6, dpi=400)


causal.dss.erp <- coxph(DSS ~ LYMPH_NODE_STATUS + TUMOR_SIZE + UBE2C_grp, fullData.list$erp$train.npn)

causal.dr.erp <- coxph(DSS ~ LYMPH_NODE_STATUS + TUMOR_SIZE + UBE2C_grp + MAPT + ELOVL5, fullData.list$erp$train.npn)


forestplot.dss.erp <- data.frame(Dataset='ER+', Outcome='DSS',
                                 Feature=c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp'),
                                 Coef=coef(causal.dss.erp),
                                 SE=sqrt(diag(vcov(causal.dss.erp))))

forestplot.dr.erp <- data.frame(Dataset='ER+', Outcome='DistantRelapse',
                                Feature=c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp',
                                          'MAPT', 'ELOVL5'),
                                Coef=coef(causal.dr.erp),
                                SE=sqrt(diag(vcov(causal.dr.erp))))



repititions <- 100

bayesianGrowShrink <- function(target, data, maxSize, gamma=0.5) {

    featVec <- setdiff(colnames(data), c(target, surv.feats))
    nEvents <- sum(data[,target][,2])
    p <- length(featVec)
    S <- c()

    ## kappa <- log(p) / log(nEvents)

    ## gamma <- max(min(ceiling(100*(1 - 1/(2*kappa))+1)/100, 1), 0)
    ## ## gamma <- 1 - 1/(2*kappa)

    ## print(paste('Kappa:', kappa))
    ## print(paste('Gamma:', gamma))

    ## gamma <- 0.5

    ## gamma <- 0

    ## gamma <- 1

    ## if (kappa < 1) {
    ##     gamma <- 0.5
    ## } else {
    ##     gamma <- 1
    ## }
    
    mod.null <- coxph(as.formula(paste0(target, ' ~ ',
                                       paste(c(1,S), collapse=' + '))),
                      data)

    bic.null <- BIC(mod.null) + 2 * gamma * log(choose(p, length(S)))
    
    while (length(S) < (nEvents-1) && length(S) < maxSize) {
        bics <- rep(0,p)
        names(bics) <- featVec
        mod.old <- coxph(as.formula(paste0(target, ' ~ ',
                                           paste(c(1,S), collapse=' + '))),
                         data)
        bic.old <- BIC(mod.old) + 2 * gamma * log(choose(p, length(S)))
        
        for (col in featVec) {
            if (col %in% S) next
            
            temp.S <- c(S, col)
            mod <- coxph(as.formula(paste0(target, ' ~ ',
                                           paste(c(1, temp.S), collapse=' + '))),
                         data)
            bics[col] <- BIC(mod) + 2 * gamma * log(choose(p, length(temp.S)))    
        }

        bics <- bics[!names(bics) %in% S]

        posterior <- c(stop=1, exp(-(bics - bic.old)/2))
        posterior <- posterior / sum(posterior)
        
        select <- sample(names(posterior), 1, prob=posterior)
        
        ## print(head(sort(posterior, decreasing=T)))
        
        if (select == 'stop') break
        
        S <- c(S, select)

        ## print(S)
        ## print(bics[select])
    }


    while (TRUE) {
        bics <- rep(0,length(S))
        names(bics) <- S
        mod.old <- coxph(as.formula(paste0(target, ' ~ ',
                                           paste(c(1,S), collapse=' + '))),
                         data)
        bic.old <- BIC(mod.old) + 2 * gamma * log(choose(p, length(S)))
        for (col in S) {
            temp.S <- setdiff(S, col)
            mod <- coxph(as.formula(paste0(target, ' ~ ',
                                           paste(c(1, temp.S), collapse=' + '))),
                         data)
            bics[col] <- BIC(mod) + 2 * gamma * log(choose(p, length(temp.S)))
        }

        ## bics <- bics[!names(bics) %in% S]

        posterior <- c(stop=1, exp(-(bics - bic.old)/2))
        posterior <- posterior / sum(posterior)
        
        select <- sample(names(posterior), 1, prob=posterior)

        ## print(head(sort(posterior, decreasing=T)))
        
        if (select == 'stop') break
        
        S <- setdiff(S, select)

        ## print(S)
        ## print(bics[select])
    }
    ## print(paste(S, collapse=', '))
    ## featCounts[S] <- featCounts[S] + 1
    if (is.null(S) | length(S)==0) S <- c()
    mod.final <- coxph(as.formula(paste0(target, ' ~ ',
                                         paste(c(1,S), collapse=' + '))),
                       data)
    if (is.null(S) | length(S)==0) S <- c('empty')
    attr(S, 'BIC') <- BIC(mod.final) + 2 * gamma * log(choose(p, length(S)))
    ## print(length(S))
    return(S)
}

mb <- bayesianGrowShrink('DFS', fullData.list$ern$train.npn, 30)
mb

summary(coxph(DistantRelapse ~ ., fullData.list$erp$train.npn[,c('DistantRelapse', mb)]))

## fullData.list$ern$train.npn$DSS.null <- fullData.list$ern$train.npn$DSS[sample(nrow(fullData.list$ern$train.npn)),]
## mb.selections.exp.null <- mclapply(rep('DSS.null', repititions), bayesianGrowShrink, data=fullData.list$ern$train.npn, maxSize=30, mc.cores=10)

library(parallel)
mb.selections.ebic <- mclapply(rep('DSS', repititions), bayesianGrowShrink, data=fullData.list$ern$train.npn, maxSize=30, mc.cores=10)

## mb.selections.exp <- mb.selections

## mb.selections <- mb.selections.ebic

featPosteriorSimple <- table(unlist(mb.selections)) / repititions

barplot(sort(featPosteriorSimple[featPosteriorSimple>0.05], decreasing=T))

mb.map <- mb.selections[[which.min(sapply(mb.selections, attr, which='BIC'))]]

mb.map

## colnames(data)

hist(sapply(mb.selections.exp, length))
hist(sapply(mb.selections.bic, length))
hist(sapply(mb.selections.ebic, length))

hist(sapply(mb.selections.exp.null, length))
hist(sapply(mb.selections.bic.null, length))

summary(sapply(mb.selections.exp, length))
summary(sapply(mb.selections.bic, length))
summary(sapply(mb.selections.ebic, length))

summary(sapply(mb.selections.exp.null, length))
summary(sapply(mb.selections.bic.null, length))

min.bic <- attr(mb.map, 'BIC')

mb.weights <- exp(-(sapply(mb.selections, attr, which='BIC') - min.bic)/2)

mb.weights <- mb.weights / sum(mb.weights)

featVec <- setdiff(colnames(fullData.list$ern$train.npn), c(surv.feats))
p <- length(featVec)    

featPosterior <- rep(0,p)
names(featPosterior) <- featVec
featPosterior <- c(empty=0, featPosterior)
for (ridx in 1:repititions) {
    featPosterior[mb.selections[[ridx]]] <- featPosterior[mb.selections[[ridx]]] + mb.weights[ridx]
}

barplot(rev(sort(featPosterior[featPosterior>0.1], decreasing=T)), horiz=T)

mbList <- list()
mbMapList <- list()
posteriorList <- list()
## for (erstatus in c('ern', 'erp')) {
##     mbList[[erstatus]] <- list()
##     mbMapList[[erstatus]] <- list()
##     posteriorList[[erstatus]] <- list()
## }

featVec <- setdiff(colnames(fullData.list$ern$train.npn), c(surv.feats))
p <- length(featVec)    

mbList <- list()
mbMapList <- list()
posteriorList <- list()

repititions <- 200

set.seed(42)

for (erstatus in c('ern', 'erp')) {
    mbList[[erstatus]] <- list()
    mbMapList[[erstatus]] <- list()
    posteriorList[[erstatus]] <- list()
    for (sf in surv.feats[c(1,3,4)]) {

        mbList[[erstatus]][[sf]] <- mclapply(
            rep(sf, repititions),
            bayesianGrowShrink,
            data=fullData.list[[erstatus]]$train.npn,
            maxSize=30,
            mc.cores=10)
        
        mbMapList[[erstatus]][[sf]] <- mbList[[erstatus]][[sf]][[which.min(sapply(mbList[[erstatus]][[sf]], attr, which='BIC'))]]

        min.bic <- attr(mbMapList[[erstatus]][[sf]], 'BIC')

        mb.weights <- exp(-(sapply(mbList[[erstatus]][[sf]], attr, which='BIC') - min.bic)/2)

        mb.weights <- mb.weights / sum(mb.weights)
        
        featPosterior <- rep(0,p)
        names(featPosterior) <- featVec
        featPosterior <- c(empty=0, featPosterior)
        for (ridx in 1:repititions) {
            featPosterior[mbList[[erstatus]][[sf]][[ridx]]] <- featPosterior[mbList[[erstatus]][[sf]][[ridx]]] + mb.weights[ridx]
        }

        posteriorList[[erstatus]][[sf]] <- featPosterior

    }
}


barplot(rev(sort(posteriorList[['ern']][['DSS']][posteriorList[['ern']][['DSS']]>0.01], decreasing=T)), horiz=T, main='ER- DSS')

barplot(rev(sort(posteriorList[['ern']][['DistantRelapse']][posteriorList[['ern']][['DistantRelapse']]>0.01], decreasing=T)), horiz=T, main='ER- DistantRelapse')

barplot(rev(sort(posteriorList[['ern']][['LocalRelapse']][posteriorList[['ern']][['LocalRelapse']]>0.01], decreasing=T)), horiz=T, main='ER- LocalRelapse')

barplot(rev(sort(posteriorList[['erp']][['DSS']][posteriorList[['erp']][['DSS']]>0.05], decreasing=T)), horiz=T, main='ER+ DSS')

barplot(rev(sort(posteriorList[['erp']][['DistantRelapse']][posteriorList[['erp']][['DistantRelapse']]>0.05], decreasing=T)), horiz=T, main='ER+ DistantRelapse')

barplot(rev(sort(posteriorList[['erp']][['LocalRelapse']][posteriorList[['erp']][['LocalRelapse']]>0.05], decreasing=T)), horiz=T, main='ER+ LocalRelapse')


for (erstatus in c('ern', 'erp')) {
    print(paste('ER_STATUS:', ifelse(erstatus=='erp', 'Positive', 'Negative')))
    for (sf in surv.feats[c(1,3,4)]) {
        print(paste('  Outcome:', sf))
        f <- as.formula(paste(sf, '~', paste(mbMapList[[erstatus]][[sf]], collapse=' + ')))

        mod <- coxph(f, fullData.list[[erstatus]]$train.npn)

        print(summary(mod))
    }
}


sampMatList <- list()
for (erstatus in c('ern', 'erp')) {
    sampMatList[[erstatus]] <- list()
    print(paste('ER_STATUS:', ifelse(erstatus=='erp', 'Positive', 'Negative')))
    for (sf in surv.feats[c(1,3,4)]) {
        print(paste('  Outcome:', sf))

        min.bic <- attr(mbMapList[[erstatus]][[sf]], 'BIC')

        mb.weights <- exp(-(sapply(mbList[[erstatus]][[sf]], attr, which='BIC') - min.bic)/2)

        uniqueFeats <- unique(unlist(
            lapply(mbList[[erstatus]][[sf]],
                   function(x) {
                       if (!all(x=='empty')) {
                           f <- as.formula(paste('~', paste(x, collapse=' + ')))
                           return(colnames(model.matrix(f, fullData.list[[erstatus]]$train.npn))[-1])
                       }
                   })))
        nSampTot <- sum(round(10000 * mb.weights))

        sampMat <- matrix(0, nSampTot, length(uniqueFeats))
        colnames(sampMat) <- uniqueFeats
        idx <- 1
        for (r in 1:repititions) {
            nSamp <- round(10000 * mb.weights[r])
            if (all(mbList[[erstatus]][[sf]][[r]]=='empty') || nSamp==0) {
                idx <- idx + nSamp
                next
            }
            
            f <- as.formula(paste(sf, '~', paste(mbList[[erstatus]][[sf]][[r]], collapse=' + ')))
            mod <- coxph(f, fullData.list[[erstatus]]$train.npn)
            ## print(summary(mod))

            beta.samps <- rmvnorm(nSamp, coef(mod), vcov(mod))

            sampMat[idx:(idx+nSamp-1), colnames(beta.samps)] <- beta.samps

            idx <- idx + nSamp
        }
        sampMatList[[erstatus]][[sf]] <- sampMat
    }
}

sort(abs(apply(sampMatList[[erstatus]][[sf]], 2, mean)))

for (erstatus in c('ern', 'erp')) {
    print(paste('ER_STATUS:', ifelse(erstatus=='erp', 'Positive', 'Negative')))
    for (sf in surv.feats[c(1,3,4)]) {
        print(paste('  Outcome:', sf))

        beta.med <- apply(sampMatList[[erstatus]][[sf]], 2, median)

        nonzeroFeats <- names(beta.med)[beta.med!=0]

        if (length(nonzeroFeats)>1) {
            print(apply(sampMatList[[erstatus]][[sf]][,nonzeroFeats], 2, quantile, probs=c(0.025, 0.05, 0.5, 0.95, 0.975)))
        } else {
            temp <- as.matrix(quantile(sampMatList[[erstatus]][[sf]][,nonzeroFeats], probs=c(0.025, 0.05, 0.5, 0.95, 0.975)))
            colnames(temp) <- nonzeroFeats
            print(temp)
        }
    }
}



#### ERN GRASP test

dim(fullData.list$ern$train)
head(fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats])

npnMat <- DataSetNPNTest(fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats])

colnames(npnMat) <- colnames(fullData.list$ern$train)[!colnames(fullData.list$ern$train) %in% surv.feats]

head(npnMat)

plot(npnMat[,'AGE_AT_DIAGNOSIS'], fullData.list$ern$train$AGE_AT_DIAGNOSIS)
plot(npnMat[,'AGE_AT_DIAGNOSIS'], fullData.list$ern$train.npn$AGE_AT_DIAGNOSIS)

hist(fullData.list$ern$train$GJA1)
plot(npnMat[,'GJA1'], fullData.list$ern$train$GJA1)
plot(npnMat[,'GJA1'], fullData.list$ern$train.npn$GJA1)

hist(fullData.list$ern$train$AGR2)
plot(npnMat[,'AGR2'], fullData.list$ern$train$AGR2)
plot(npnMat[,'AGR2'], fullData.list$ern$train.npn$AGR2)


table(npnMat[,'HER2_STATUS'], fullData.list$ern$train$HER2_STATUS)

mb1 <- growShrinkMB(fullData.list$ern$train.npn[,!colnames(fullData.list$ern$train) %in% surv.feats], 'LYMPH_NODE_STATUS', verbose=T)

mb1

mb2 <- GrowShrinkTreeTest(fullData.list$ern$train.npn[,!colnames(fullData.list$ern$train) %in% surv.feats], 'LYMPH_NODE_STATUS')

mb2

setdiff(mb1, mb2)

mb2 <- growShrinkMB(fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats], 'LYMPH_NODE_STATUS', rank=T, verbose=T, penalty=2)

mb2

mb3 <- growShrinkMB(fullData.list$ern$train.npn[,!colnames(fullData.list$ern$train) %in% surv.feats], 'LYMPH_NODE_STATUS', rank=T, verbose=T)

mb3


setdiff(mb1, mb2)

setdiff(mb2, mb1)

setdiff(mb3, mb2)

setdiff(mb2, mb3)



mb1 <- growShrinkMB(fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats], 'PR_STATUS', verbose=T, penalty=2)

mb1

mb2 <- growShrinkMB(fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats], 'PR_STATUS', rank=T, verbose=T, penalty=2)

mb2

mb3 <- growShrinkMB(fullData.list$ern$train.npn[,!colnames(fullData.list$ern$train) %in% surv.feats], 'PR_STATUS', rank=F, verbose=T, penalty=2)

mb3

plot(fullData.list$ern$train$PR_STATUS, fullData.list$ern$train$PGR)
plot(fullData.list$ern$train$PR_STATUS, fullData.list$ern$train.npn$PGR)
plot(fullData.list$ern$train$PR_STATUS, npnMat[,'PGR'])

apply(fullData.list$ern$train.npn, 2, sd)

setdiff(mb1, mb2)

setdiff(mb2, mb1)

setdiff(mb3, mb2)

setdiff(mb2, mb3)

res2 <- glm(paste('PR_STATUS ~', paste(mb2, collapse=" + ")),
            family=binomial(), npnMat %>% as.data.frame %>% mutate_at("PR_STATUS", factor))

summary(res2)



res3 <- glm(paste('PR_STATUS ~', paste(mb3, collapse=" + ")),
            family=binomial(), fullData.list$ern$train.npn)

summary(res3)


kappa <- log(ncol(fullData.list$ern$train)-7) / log(nrow(fullData.list$ern$train))

kappa

gamma <- (1-(1/(4*kappa)))

gamma

## gamma <- (1/(4*kappa))

gamma <- 0.25

ebicPenalty <- 1 + 4 * gamma * kappa

ebicPenalty

system.time({
    set.seed(43)
    g.grasp.ernd1r <- grasp(fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats],
                            depth=1, numStarts=3, penalty=2, rank=T, verbose=T)
})

g.grasp.ern2d2
attr(g.grasp.ern2d2, "Score")

degree <- sapply(g.grasp.ern2$neighbors, length)

hist(degree)

system.time({
    g.boss.ern <- annealboss(fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats], numStarts=1, penalty=2, rank=T, verbose=T)
})


system.time({
    g.boss.ern <- boss(fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats], numStarts=5, penalty=2, rank=T, verbose=T)
})

system.time({
    g.boss.erp <- boss(fullData.list$erp$train[,!colnames(fullData.list$erp$train) %in% surv.feats], numStarts=5, penalty=2, rank=T, verbose=T)
})

g.boss.ern
attr(g.boss.ern2v2, "Score")

degree <- sapply(g.boss.ern$neighbors, length)

hist(degree)
sort(degree)

g.boss.erp
attr(g.boss.erp2, "Score")

degree <- sapply(g.boss.erp$neighbors, length)

hist(degree)
sort(degree)

par(mfrow=c(3,1))
plot(g.boss.ern, c("GRADE", g.boss.ern$neighbors$GRADE), list(fontsize=32))
plot(g.boss.ern, c("LYMPH_NODE_STATUS", g.boss.ern$neighbors$LYMPH_NODE_STATUS), list(fontsize=32))
plot(g.boss.ern, c("TUMOR_SIZE", g.boss.ern$neighbors$TUMOR_SIZE), list(fontsize=32))


plot(g.boss.ern, c("HER2_STATUS", g.boss.ern$neighbors$HER2_STATUS), list(fontsize=32))
plot(g.boss.ern, c("FGFR4", g.boss.ern$neighbors$FGFR4), list(fontsize=32))
plot(g.boss.ern, c("POSTN", g.boss.ern$neighbors$POSTN), list(fontsize=32))
plot(g.boss.ern, c("VCAM1", g.boss.ern$neighbors$VCAM1), list(fontsize=24))

plot(g.boss.ern, c("S100A9_grp", g.boss.ern$neighbors$S100A9_grp), list(fontsize=32))

common <- intersect(g.boss.ern$neighbors$UBE2C_grp, g.boss.erp$neighbors$UBE2C_grp)

par(mfrow=c(3,1))
plot(g.boss.erp, c("GRADE", g.boss.erp$neighbors$GRADE), list(fontsize=32))
plot(g.boss.erp, c("LYMPH_NODE_STATUS", g.boss.erp$neighbors$LYMPH_NODE_STATUS), list(fontsize=32))
plot(g.boss.erp, c("TUMOR_SIZE", g.boss.erp$neighbors$TUMOR_SIZE), list(fontsize=32))


plot(g.boss.erp, c("HER2_STATUS", g.boss.erp$neighbors$HER2_STATUS), list(fontsize=32))
plot(g.boss.erp, c("FGFR4", g.boss.erp$neighbors$FGFR4), list(fontsize=32))
plot(g.boss.erp, c("POSTN", g.boss.erp$neighbors$POSTN), list(fontsize=32))
plot(g.boss.erp, c("VCAM1", g.boss.erp$neighbors$VCAM1), list(fontsize=24))

plot(g.boss.erp, c("S100A9_grp", g.boss.erp$neighbors$S100A9_grp), list(fontsize=32))


SHD(g.grasp.ern2, g.boss.ern2)
SHD(g.grasp.ern2d2, g.boss.ern2)
SHD(g.grasp.ern2, g.grasp.ern2d2)

setdiff(g.grasp.ern2d2$edges, g.grasp.ern2$edges)

setdiff(g.grasp.ern2$edges, g.boss.ern2$edges)
setdiff(g.boss.ern2$edges, g.grasp.ern2$edges)

par(mfrow=c(2,1))
plot(g.boss.ern, c("AGE_AT_DIAGNOSIS", g.boss.ern$neighbors$AGE_AT_DIAGNOSIS), list(fontsize=64))
plot(g.boss.erp, c("AGE_AT_DIAGNOSIS", g.boss.erp$neighbors$AGE_AT_DIAGNOSIS), list(fontsize=64))

plot(g.boss.ern, c("ESR1", g.boss.ern$neighbors$ESR1), list(fontsize=32))
plot(g.boss.erp, c("ESR1", g.boss.erp$neighbors$ESR1), list(fontsize=80))

plot(g.boss.ern, c("UBE2C_grp", common), list(fontsize=32))
plot(g.boss.erp, c("UBE2C_grp", common), list(fontsize=32))



plot(g.grasp.ern2d2, c("AGE_AT_DIAGNOSIS", g.grasp.ern2d2$neighbors$AGE_AT_DIAGNOSIS), list(fontsize=32))

g.boss.ern

g.grasp.ern3

g.grasp.ern
g.grasp.ern$graphs[[1]]

sort(sapply(g.boss.ern$neighbors, length))

par(mfrow=c(4,1))
plot(g.grasp.ern$graphs[[1]], "PR_STATUS", list(fontsize=32))
plot(g.grasp.ern2, "PR_STATUS", list(fontsize=32))
plot(g.grasp.ern3, "PR_STATUS", list(fontsize=32))
plot(g.boss.ern, c("LYMPH_NODE_STATUS", g.boss.ern$neighbors$LYMPH_NODE_STATUS), list(fontsize=32))

par(mfrow=c(2,1))
plot(g.grasp.ern2, unique(c("HER2_STATUS", "ERBB2", g.grasp.ern2$neighbors$HER2_STATUS, g.grasp.ern2$neighbors$ERBB2)), list(fontsize=32))
plot(g.grasp.ern3, unique(c("HER2_STATUS", "ERBB2", g.grasp.ern3$neighbors$HER2_STATUS, g.grasp.ern3$neighbors$ERBB2)), list(fontsize=32))
plot(g.boss.ern, unique(c("HER2_STATUS", "ERBB2", g.boss.ern$neighbors$HER2_STATUS, g.grasp.ern3$neighbors$ERBB2)), list(fontsize=32))


sort(sapply(g.grasp.ern$graphs[[1]]$neighbors, length))
sort(sapply(g.grasp.ern2$neighbors, length))
sort(sapply(g.grasp.ern3$neighbors, length))

sort(sapply(fullGraph.list$ern$neighbors, length))

fullGraph.list$ern$neighbors$DSS
fullGraph.list$ern$neighbors$DistantRelapse


kappa <- log(ncol(fullData.list$erp$train)-7) / log(nrow(fullData.list$erp$train))

kappa

gamma <- (1-(1/(4*kappa)))

gamma

## gamma <- (1/(4*kappa))

gamma <- 0.5

ebicPenalty <- 1 + 4 * gamma * kappa

ebicPenalty

system.time(
    g.grasp.erp <- grasp(fullData.list$erp$train[,!colnames(fullData.list$erp$train) %in% surv.feats],
                         depth=1, numStarts=2, penalty=ebicPenalty, rank=T, verbose=T)
)

g.grasp.erp

sort(sapply(g.grasp.erp$graphs[[1]]$neighbors, length))

sort(sapply(fullGraph.list$erp$neighbors, length))

fullGraph.list$erp$neighbors$DSS
fullGraph.list$erp$neighbors$DistantRelapse

par(mfrow=c(2,1))
feat <- "HER2_STATUS"
plot(g.grasp.erp$graphs[[1]], feat, list(fontsize=32))
plot(fullGraph.list$erp, feat, list(fontsize=32))


par(mfrow=c(2,1))
feat <- "UBE2C_grp"
plot(g.grasp.erp$graphs[[1]], c(feat, g.grasp.erp$graphs[[1]]$neighbors[[feat]]), list(fontsize=32))
plot(fullGraph.list$erp, c(feat, fullGraph.list$erp$neighbors[[feat]]), list(fontsize=32))




#### Network Topology

deg.erp <- sapply(fullGraph.list$erp$neighbors, length)

deg.ern <- sapply(fullGraph.list$ern$neighbors, length)

plot(deg.erp, deg.ern)

cor.test(deg.erp, deg.ern)

toIGraph <- function(graph, directed=FALSE) {
    require(igraph)
    g <- make_empty_graph(directed=directed)
    for (n in graph$nodes) {
        g <- g + vertex(n)
    }
    for (e in graph$edges) {
        e.split <- strsplit(e, " ")[[1]]
        if (directed && e.split[2]=="-->") {
            g <- g + edge(e.split[1], e.split[3])
        } else if (directed) {
            g <- g + edge(e.split[1], e.split[3])
            g <- g + edge(e.split[3], e.split[1])
        } else {
            g <- g + edge(e.split[1], e.split[3])
        }
    }
    g
}

igr.erp <- toIGraph(fullGraph.list$erp, directed=F)

top.erp <- data.frame(Degree=degree(igr.erp),
                      Closeness=closeness(igr.erp),
                      Betweenness=betweenness(igr.erp),
                      HarmonicCentrality=harmonic_centrality(igr.erp),
                      EigenCentrality=eigen_centrality(igr.erp)$vector,
                      HubScore=hub_score(igr.erp)$vector,
                      AuthorityScore=authority_score(igr.erp)$vector,
                      PageRank=page_rank(igr.erp, directed=F)$vector,
                      Knn=knn(igr.erp)$knn)

top.erp <- scoreRegulators(top.erp)

head(top.erp %>% arrange(desc(Degree)), 10)

head(top.erp %>% arrange(desc(Closeness)), 10)

head(top.erp %>% arrange(desc(Betweenness)), 10)

head(top.erp %>% arrange(desc(HarmonicCentrality)), 10)

head(top.erp %>% arrange(desc(EigenCentrality)), 10)

head(top.erp %>% arrange(desc(HubScore)), 10)

head(top.erp %>% arrange(desc(AuthorityScore)), 10)

head(top.erp %>% arrange(desc(PageRank)), 10)

head(top.erp %>% arrange(Knn), 10)

head(top.erp %>% arrange(desc(RegulatorScore)), 10)

ggsurvplot(survfit(DRFS ~ quantcut(UBE2C_grp,2), fullData.list$erp$train), pval=T, conf.int=T)
ggsurvplot(survfit(DRFS ~ quantcut(CYBRD1,2), fullData.list$erp$train), pval=T, conf.int=T)
ggsurvplot(survfit(DRFS ~ quantcut(NFIX,2), fullData.list$erp$train), pval=T, conf.int=T)
ggsurvplot(survfit(DRFS ~ quantcut(STC2,2), fullData.list$erp$train), pval=T, conf.int=T)


plot(top.erp$PageRank, top.ern$PageRank)
cor.test(top.erp$PageRank, top.ern$PageRank)

plot(sqrt(top.erp$RegulatorScore), sqrt(top.ern$RegulatorScore))
cor.test(sqrt(top.erp$RegulatorScore), sqrt(top.ern$RegulatorScore))

top.erp[top.erp$RegulatorScore >= 12 & top.ern$RegulatorScore < 12,] %>% arrange(desc(RegulatorScore))

top.ern[top.ern$RegulatorScore >= 12 & top.erp$RegulatorScore < 12,] %>% arrange(desc(RegulatorScore))

top.erp[top.erp$RegulatorScore >= 12 & top.ern$RegulatorScore >= 12,] %>% arrange(desc(RegulatorScore))

top.ern[top.ern$RegulatorScore >= 12 & top.erp$RegulatorScore >= 12,] %>% arrange(desc(RegulatorScore))

top.erp[top.erp$PageRank >= sqrt(1e-5) & top.ern$PageRank < sqrt(1e-5),] %>% arrange(desc(PageRank))

top.ern[top.ern$PageRank >= sqrt(1e-5) & top.erp$PageRank < sqrt(1e-5),] %>% arrange(desc(PageRank))

top.erp[top.erp$PageRank >= sqrt(1e-5) & top.ern$PageRank >= sqrt(1e-5),] %>% arrange(desc(PageRank))

top.ern[top.ern$PageRank >= sqrt(1e-5) & top.erp$PageRank >= sqrt(1e-5),] %>% arrange(desc(PageRank))



sdBin <- function(x) {
    z <- scale(x)
    bins <- ifelse(z <= -2, "A",
            ifelse(z <= -1, "B",
            ifelse(z <= 0, "C",
            ifelse(z <= 1, "D",
            ifelse(z <= 2, "E", "F")))))
    return(factor(bins, levels=c("A","B","C","D","E","F")))
}

classifyRegulators <- function(top.df) {
    KnnBins <- sdBin(top.df$Knn)
    DegreeBins <- sdBin(top.df$Degree)
    PageRankBins <- sdBin(top.df$PageRank)

    preds <- rep(NA, nrow(top.df))
    names(preds) <- rownames(top.df)
    preds[KnnBins%in%c("A","B")] <- "Regulator"
    preds[KnnBins%in%c("D","E","F")] <- "Target"

    preds[KnnBins=="C" & PageRankBins%in%c("D","E","F")] <- "Regulator"
    preds[KnnBins=="C" & PageRankBins%in%c("A", "B")] <- "Target"

    preds[KnnBins=="C" & PageRankBins == "C" & DegreeBins%in%c("D","E","F")] <- "Regulator"
    preds[KnnBins=="C" & PageRankBins == "C" & DegreeBins%in%c("A","B","C")] <- "Target"

    return(preds)
}

scoreRegulators <- function(top.df) {
    n <- nrow(top.df)
    ## KnnPvals <- rank(top.df$Knn) / (n+1)
    DegreePvals <- rank(-top.df$Degree) / (n+1)
    PageRankPvals <- rank(-top.df$PageRank) / (n+1)

    score <- -2 * (log(DegreePvals) + log(PageRankPvals)) ## log(KnnPvals) + 
    names(score) <- rownames(top.df)

    top.df$RegulatorScore <- score
    top.df$RegulatorPval <- pchisq(score, 4, lower.tail=F)
    top.df$RegulatorFDR <- p.adjust(top.df$RegulatorPval, method='fdr')
    return(top.df)
}



reg.res <- scoreRegulators(top.erp)

reg.res$padj <- p.adjust(reg.res$pval, method='fdr')

reg.res$padj[reg.res$padj<0.05]

regs <- classifyRegulators(top.erp)
table(regs)
regs[regs=="Regulator"]




hist(top.erp$Degree)
hist(top.erp$Knn)
hist(top.erp$PageRank)

hist(floor(scale(top.erp$Degree))+2)

plot(igr.erp, vertex.size=top.erp$EigenCentrality * 10)
plot(igr.erp, vertex.size=top.erp$Betweenness / max(top.erp$Betweenness) * 10)
plot(igr.erp, vertex.size=top.erp$PageRank / max(top.erp$PageRank) * 10)

ceb.erp <- cluster_edge_betweenness(igr.erp)

dendPlot(ceb.erp, mode="hclust")

modules.erp <- list()
for (mod in 1:length(ceb.erp)) {
    modules.erp[[mod]] <- sort(ceb.erp$names[ceb.erp$membership==mod])
}

modules.erp


igr.ern <- toIGraph(fullGraph.list$ern, directed=F)

top.ern <- data.frame(Degree=degree(igr.ern),
                      Closeness=closeness(igr.ern),
                      Betweenness=betweenness(igr.ern),
                      HarmonicCentrality=harmonic_centrality(igr.ern),
                      EigenCentrality=eigen_centrality(igr.ern)$vector,
                      HubScore=hub_score(igr.ern)$vector,
                      AuthorityScore=authority_score(igr.ern)$vector,
                      PageRank=page_rank(igr.ern)$vector,
                      Knn=knn(igr.ern)$knn)

top.ern <- scoreRegulators(top.ern)

head(top.ern %>% arrange(desc(Degree)), 10)

head(top.ern %>% arrange(desc(Closeness)), 10)

head(top.ern %>% arrange(desc(Betweenness)), 10)

head(top.ern %>% arrange(desc(HarmonicCentrality)), 10)

head(top.ern %>% arrange(desc(EigenCentrality)), 10)

head(top.ern %>% arrange(desc(HubScore)), 10)

head(top.ern %>% arrange(desc(AuthorityScore)), 10)

head(top.ern %>% arrange(desc(PageRank)), 10)

head(top.ern %>% arrange(Knn), 10)

head(top.ern %>% arrange(desc(RegulatorScore)), 10)


ggsurvplot(survfit(DRFS ~ quantcut(CD3D_grp,2), fullData.list$ern$train), pval=T, conf.int=T)
ggsurvplot(survfit(DRFS ~ quantcut(PIP,2), fullData.list$ern$train), pval=T, conf.int=T)
ggsurvplot(survfit(DRFS ~ quantcut(AGR2,2), fullData.list$ern$train), pval=T, conf.int=T)
ggsurvplot(survfit(DRFS ~ quantcut(VTCN1,2), fullData.list$ern$train), pval=T, conf.int=T)

plot(igr.ern, vertex.size=top.ern$EigenCentrality * 10)
plot(igr.ern, vertex.size=top.ern$Betweenness / max(top.ern$Betweenness) * 10)
plot(igr.ern, vertex.size=top.ern$PageRank / max(top.ern$PageRank) * 10)

ceb.ern <- cluster_edge_betweenness(igr.ern)

dendPlot(ceb.ern, mode="hclust")

modules.ern <- list()
for (mod in 1:length(ceb.ern)) {
    modules.ern[[mod]] <- sort(ceb.ern$names[ceb.ern$membership==mod])
}

modules.ern


getPa <- function(graph, target) {
    pa <- c()
    for (edge in graph$edges) {
        edgesplit <- strsplit(edge, " ")[[1]]
        if (edgesplit[3] == target) {
            pa <- c(pa, edgesplit[1])
        } else if (edgesplit[2] == target && edgesplit[1] == target) {
            pa <- c(pa, edgesplit[3])
        }
    }
    pa
}


ern.effects <- rep(NA, 437)
names(ern.effects) <- colnames(fullData.list$ern$train.npn)[9:445]

erp.effects <- rep(NA, 437)
names(erp.effects) <- colnames(fullData.list$erp$train.npn)[9:445]

ern.zs <- rep(NA, 437)
names(ern.zs) <- colnames(fullData.list$ern$train.npn)[9:445]

erp.zs <- rep(NA, 437)
names(erp.zs) <- colnames(fullData.list$erp$train.npn)[9:445]


for (gene in colnames(fullData.list$ern$train.npn)[9:445]) {
    f.ern <- as.formula(paste("DFS ~", paste(c(gene, getPa(g.boss.ern, gene)), collapse=" + ")))

    ern.table <- summary(coxph(f.ern, fullData.list$ern$train.npn))$coefficients


    f.erp <- as.formula(paste("DFS ~ +", paste(c(gene, getPa(g.boss.erp, gene)), collapse=" + ")))

    erp.table <- summary(coxph(f.erp, fullData.list$erp$train.npn))$coefficients

    ern.effects[gene] <- ern.table[1,1]

    erp.effects[gene] <- erp.table[1,1]

    ern.zs[gene] <- ern.table[1,4]

    erp.zs[gene] <- erp.table[1,4]
}

tail(sort(abs(ern.effects)), 20)

tail(sort(abs(ern.zs)), 20)

tail(sort(abs(erp.effects)), 20)

tail(sort(abs(erp.zs)), 20)
