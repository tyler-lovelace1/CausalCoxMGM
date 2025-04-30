library(rCausalMGM)
library(tidyverse)
library(survival)
library(survminer)
library(glmnet)
library(clinfun)
library(mstate)
library(parallel)

loadData <- function(erstatus='all', k=-1) {

    if (k < 0) {
        cvString <- 'full'
        test <- NULL
    } else {
        cvString <- paste0('cv', k)
    }

    if (erstatus=='all') {
        if (cvString=='full') {
            train <- read.csv(paste0('data/metabric.rna.', cvString, '.csv'), row.names=1)
        } else {
            train <- read.csv(paste0('data/metabric.rna.train.', cvString, '.csv'), row.names=1)
            test <- read.csv(paste0('data/metabric.rna.test.', cvString, '.csv'), row.names=1)
        }
    } else {
        if (cvString=='full') {
            train <- read.csv(paste0('data/metabric.rna.', erstatus, '.', cvString, '.csv'), row.names=1)
        } else {
            train <- read.csv(paste0('data/metabric.rna.', erstatus, '.train.', cvString, '.csv'), row.names=1)
            test <- read.csv(paste0('data/metabric.rna.', erstatus, '.test.', cvString, '.csv'), row.names=1)
        }
    }

    if (any(c('DSS.time', 'DSS.status') %in% colnames(train))) {
        train$DSS <- Surv(train$DSS.time, train$DSS.status)
        if (cvString!='full') {
            test$DSS <- Surv(test$DSS.time, test$DSS.status)
        }
    }
    if (any(c('OD.time', 'OD.status') %in% colnames(train))) {
        train$OD <- Surv(train$OD.time, train$OD.status)
        if (cvString!='full') {
            test$OD <- Surv(test$OD.time, test$OD.status)
        }
    }
    if (any(c('DistantRelapse.time', 'DistantRelapse.status') %in% colnames(train))) {
        train$DistantRelapse <- Surv(train$DistantRelapse.time, train$DistantRelapse.status)
        if (cvString!='full') {
            test$DistantRelapse <- Surv(test$DistantRelapse.time, test$DistantRelapse.status)
        }
    }
    if (any(c('LocalRelapse.time', 'LocalRelapse.status') %in% colnames(train))) {
        train$LocalRelapse <- Surv(train$LocalRelapse.time, train$LocalRelapse.status)
        if (cvString!='full') {
            test$LocalRelapse <- Surv(test$LocalRelapse.time, test$LocalRelapse.status)
        }
    }
    if (any(c('OS.time', 'OS.status') %in% colnames(train))) {
        train$OS <- Surv(train$OS.time, train$OS.status)
        if (cvString!='full') {
            test$OS <- Surv(test$OS.time, test$OS.status)
        }
    }
    if (any(c('DRFS.time', 'DRFS.status') %in% colnames(train))) {
        train$DRFS <- Surv(train$DRFS.time, train$DRFS.status)
        if (cvString!='full') {
            test$DRFS <- Surv(test$DRFS.time, test$DRFS.status)
        }
    }
    if (any(c('DFS.time', 'DFS.status') %in% colnames(train))) {
        train$DFS <- Surv(train$DFS.time, train$DFS.status)
        if (cvString!='full') {
            test$DFS <- Surv(test$DFS.time, test$DFS.status)
        }
    }
    if (any(c('BONE.time', 'BONE.status') %in% colnames(train))) {
        train$BONE <- Surv(train$BONE.time, train$BONE.status)
        if (cvString!='full') {
            test$BONE <- Surv(test$BONE.time, test$BONE.status)
        }
    }
    if (any(c('BRAIN.time', 'BRAIN.status') %in% colnames(train))) {
        train$BRAIN <- Surv(train$BRAIN.time, train$BRAIN.status)
        if (cvString!='full') {
            test$BRAIN <- Surv(test$BRAIN.time, test$BRAIN.status)
        }
    }
    if (any(c('LIVER.time', 'LIVER.status') %in% colnames(train))) {
        train$LIVER <- Surv(train$LIVER.time, train$LIVER.status)
        if (cvString!='full') {
            test$LIVER <- Surv(test$LIVER.time, test$LIVER.status)
        }
    }
    if (any(c('LNS.time', 'LNS.status') %in% colnames(train))) {
        train$LNS <- Surv(train$LNS.time, train$LNS.status)
        if (cvString!='full') {
            test$LNS <- Surv(test$LNS.time, test$LNS.status)
        }
    }
    if (any(c('LOCOREGIONAL.time', 'LOCOREGIONAL.status') %in% colnames(train))) {
        train$LOCOREGIONAL <- Surv(train$LOCOREGIONAL.time, train$LOCOREGIONAL.status)
        if (cvString!='full') {
            test$LOCOREGIONAL <- Surv(test$LOCOREGIONAL.time, test$LOCOREGIONAL.status)
        }
    }
    if (any(c('LUNG.time', 'LUNG.status') %in% colnames(train))) {
        train$LUNG <- Surv(train$LUNG.time, train$LUNG.status)
        if (cvString!='full') {
            test$LUNG <- Surv(test$LUNG.time, test$LUNG.status)
        }
    }
    if (any(c('PLEURA.time', 'PLEURA.status') %in% colnames(train))) {
        train$PLEURA <- Surv(train$PLEURA.time, train$PLEURA.status)
        if (cvString!='full') {
            test$PLEURA <- Surv(test$PLEURA.time, test$PLEURA.status)
        }
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

    train <- train %>% mutate_if(is.character, factor)

    numeric.mask <- sapply(train, function(x) !is.Surv(x) & is.numeric(x))

    train.npn <- train
    train.npn[,numeric.mask] <- huge::huge.npn(as.matrix(train[,numeric.mask]))

    npnFunctions <- list()
    for (feat in colnames(train)[numeric.mask]) {
        npnFunctions[[feat]] <- approxfun(x=train[,feat], y=train.npn[,feat], method='constant',
                                          yleft=min(train.npn[,feat]), yright=max(train.npn[,feat]))
    }    

    train.means <- apply(train[,numeric.mask], 2, mean)
    train.sds <- apply(train[,numeric.mask], 2, sd)

    train.zscore <- train
    train.zscore[,numeric.mask] <- scale(train.zscore[,numeric.mask])

    if (cvString != 'full') {
        test <- test[,!colnames(test) %in% c('DSS.time', 'DSS.status',
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
            test[test$GRADE!='G3','GRADE'] <- 'G1.G2'
            test[!test$ONCOTREE_CODE!='IDC','ONCOTREE_CODE'] <- 'BREAST'
        }

        factor.mask <- sapply(train, is.factor)

        for (col in colnames(train)[factor.mask]) {
            test[,col] <- factor(test[,col], levels=levels(train[,col]))
        }

        ## test <- test %>% mutate_if(is.character, factor)

        test.npn <- test
        ## test.npn[,numeric.mask] <- huge::huge.npn(as.matrix(test[,numeric.mask]))

        for (feat in colnames(train)[numeric.mask]) {
            test.npn[,feat] <- npnFunctions[[feat]](test[,feat])
        }    


        test.zscore <- test
        test.zscore[,numeric.mask] <- apply(test[,numeric.mask], 1,
                                            function(x) (x-train.means)/train.sds)

    } else {
        test <- NULL
        test.npn <- NULL
        test.zscore <- NULL
    }

    return(list(train=train, train.npn=train.npn, train.zscore=train.zscore,
                test=test, test.npn=test.npn, test.zscore=test.zscore))

}

loadGraphs <- function(erstatus='all', k=-1) {
    splitStatus <- c('split', 'composite')
    mgmSelect <- c('BIC', 'EBIC', 'StEPS', 'StARS')
    orientRules <- c('maxp')
    alphas <- c(0.05, 0.1, 0.2)

    if (k < 0) {
        cvString <- 'full'
    } else {
        cvString <- paste0('cv', k)
    }

    graph.list <- list()
    for (ss in splitStatus) {
        graph.list[[ss]] <- list()
        for (mgms in mgmSelect) {
            graph.list[[ss]][[mgms]] <- list()
            for (rule in orientRules) {
                graph.list[[ss]][[mgms]][[rule]] <- list()
                for (alpha in alphas) {
                    g <- loadGraph(
                        paste0('out/', ifelse(rule=='majority', 'mfci', 'fcimax'),
                               'mgm', mgms, '.metabric.rna.', ss, '.',
                               erstatus, '.rank.FDR',
                               gsub('0[.]', '', as.character(alpha)),
                               '.', cvString, '.txt'))
                    graph.list[[ss]][[mgms]][[rule]][[paste0('FDR',gsub('0[.]','',alpha))]] <- g
                }
            }
        }
    }
    return(graph.list)
}

getNeighbors <- function(graph, target) {
    splitedges <- strsplit(graph$edges[grep(target, graph$edges)], ' ')
    neighbors <- setdiff(unique(as.vector(unlist(sapply(splitedges,
                                                        function(x) {
                                                            if (target %in% x) {
                                                                return(x[c(1,3)])
                                                            }
                                                        })))), target)
    if (is.null(neighbors)) neighbors <- as.character(c())
    return(neighbors)
}


evaluateLASSO <- function(train, test, erstatus, times) {
    require(survival)
    require(pracma)
    require(survcomp)
    require(glmnet)
    require(survAUC)

    for (idx in 1:ncol(train)) {
        if (is.factor(train[,idx])) {
            test[,idx] <- factor(test[,idx], levels=levels(train[,idx]))
        }
    }

    if (erstatus!='relapse') {
        surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')
    } else {
        surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'OS', 'DRFS', 'DFS', 'BONE',
                        'BRAIN', 'LIVER', 'LNS', 'LOCOREGIONAL', 'LUNG', 'PLEURA')
    }
    
    eval.df <- data.frame()

    for (sf in surv.feats) {
        print(sf)
        
        nobs <- nrow(train)
        nvars <- ncol(train)-7
        
        cv.out <- cv.glmnet(model.matrix(~., train[,!colnames(train)%in%surv.feats])[,-1],
                            train[,sf],
                            family='cox',
                            lambda.min.ratio=ifelse(nobs < nvars, 0.01, 1e-03))

        plot(cv.out)

        risk.min <- predict(cv.out,
                            newx=model.matrix(~., test[,!colnames(test)%in%surv.feats])[,-1],
                            s=cv.out$lambda.min)

        test.logpl <- logpl(risk.min,
                            test[,sf][,1],
                            test[,sf][,2])

        
        if (length(unique(risk.min))!=1) {
            ## test.roc <- timeROC(test[,sf][,1], test[,sf][,2],
            ##                     risk.min,
            ##                     cause=1, times=test[,sf][,1],
            ##                     weighting='cox')
            
            ## iauc <- trapz(test.roc$times[!is.na(test.roc$AUC)],
            ##               test.roc$AUC[!is.na(test.roc$AUC)]) /
            ##     (max(test.roc$times[!is.na(test.roc$AUC)]) -
            ##      min(test.roc$times[!is.na(test.roc$AUC)]))

            test.auc.uno <- AUC.uno(train[,sf], test[,sf], risk.min, times=sort(test[,sf][,1]))

            iauc <- test.auc.uno$iauc

            uno.conc <- UnoC(train[,sf], test[,sf], risk.min)
                           
        } else {
            iauc <- 0.5
            uno.conc <- 0.5
        }

        if (erstatus=='all') {
            uno <- concordance(test[,sf] ~ strata(test$ER_STATUS) + risk.min,
                               timewt='n/G2')
        
            harrell <- concordance(test[,sf] ~ strata(test$ER_STATUS) + risk.min)
        } else {
            uno <- concordance(test[,sf] ~ risk.min,
                               timewt='n/G2')
        
            harrell <- concordance(test[,sf] ~ risk.min)
        }

        feats <- dimnames(coef(cv.out, s=cv.out$lambda.min))[[1]][as.vector(coef(cv.out, s=cv.out$lambda.min)!=0)]

        numfeats <- length(feats)

        for (fact.feat in c('GRADE', 'ONCOTREE_CODE', 'LYMPH_NODE_STATUS')) {
            feat.names <- paste0(fact.feat, unique(test[,fact.feat]))
            if (sum(feat.names %in% feats) > 1) {
                numfeats <- numfeats - sum(feat.names %in% feats) + 1
            }
        }
        
        
        eval.df <- rbind(eval.df,
                         data.frame(Outcome = sf,
                                    Features = 'Min',
                                    Time = 'Full',
                                    Metric = "log Partial Likelihood",
                                    Value = test.logpl[1] / test.logpl[2], ## / nrow(test),
                                    check.names=F),
                         data.frame(Outcome = sf,
                                    Features = 'Min',
                                    Time = 'Full',
                                    Metric = "Uno's Concordance",
                                    Value = uno.conc,
                                    check.names=F),
                         data.frame(Outcome = sf,
                                    Features = 'Min',
                                    Time = 'Full',
                                    Metric = "Harrell's Concordance",
                                    Value = 1-harrell$concordance,
                                    check.names=F),
                         data.frame(Outcome = sf,
                                    Features = 'Min',
                                    Time = 'Full',
                                    Metric = "Integrated C/D AUC",
                                    Value = iauc,
                                    check.names=F),
                         data.frame(Outcome = sf,
                                    Features = 'Min',
                                    Time = 'Full',
                                    Metric = "Number of Features",
                                    Value = numfeats,
                                    check.names=F))

        
        for (tt in times) {

            if (erstatus=='all') {
                uno <- concordance(test[,sf] ~ strata(test$ER_STATUS) + risk.min,
                                   timewt='n/G2', ymax=365.25*tt)
                
                harrell <- concordance(test[,sf] ~ strata(test$ER_STATUS) + risk.min,
                                       ymax=365.25*tt)
            } else {
                uno <- concordance(test[,sf] ~ risk.min,
                                   timewt='n/G2', ymax=365.25*tt)
                
                harrell <- concordance(test[,sf] ~ risk.min, ymax=365.25*tt)
            }


            if (length(unique(risk.min))!=1) { 

                ## test.roc <- timeROC(test[,sf][,1], test[,sf][,2],
                ##                     risk.min, cause=1,
                ##                     times=365.25*tt, weighting='cox')

                ## auc <- test.roc$AUC[2]

                test.auc.uno <- AUC.uno(train[,sf], test[,sf], risk.min, times=365.25*tt)

                auc <- test.auc.uno$auc

                ## uno.conc <- UnoC(train[,sf], test[,sf], risk.min, 365.25 * tt)
                                
                
            } else {
                auc <- 0.5
                uno.conc <- 0.5
            }
            
            eval.df <- rbind(eval.df,
                             data.frame(Outcome = sf,
                                        Features = 'Min',
                                        Time = paste(tt, 'Years'),
                                        Metric = "Uno's Concordance",
                                        Value = 1-uno$concordance, ##uno.conc,
                                        check.names=F),
                             data.frame(Outcome = sf,
                                        Features = 'Min',
                                        Time = paste(tt, 'Years'),
                                        Metric = "Harrell's Concordance",
                                        Value = 1-harrell$concordance,
                                        check.names=F),
                             data.frame(Outcome = sf,
                                        Features = 'Min',
                                        Time = paste(tt, 'Years'),
                                        Metric = "C/D AUC",
                                        Value = auc,
                                        check.names=F))
            
        }

        risk.1se <- predict(cv.out,
                            newx=model.matrix(~., test[,!colnames(test)%in%surv.feats])[,-1],
                            s=cv.out$lambda.1se)

        test.logpl <- logpl(risk.1se,
                            test[,sf][,1],
                            test[,sf][,2])

        if (length(unique(risk.1se))!=1) {
            ## test.roc <- timeROC(test[,sf][,1], test[,sf][,2],
            ##                     risk.1se,
            ##                     cause=1, times=test[,sf][,1],
            ##                     weighting='cox')
            
            ## iauc <- trapz(test.roc$times[!is.na(test.roc$AUC)],
            ##               test.roc$AUC[!is.na(test.roc$AUC)]) /
            ##     (max(test.roc$times[!is.na(test.roc$AUC)]) -
            ##      min(test.roc$times[!is.na(test.roc$AUC)]))

            test.auc.uno <- AUC.uno(train[,sf], test[,sf], risk.1se, times=sort(test[,sf][,1]))

            iauc <- test.auc.uno$iauc

            ## uno.conc <- UnoC(train[,sf], test[,sf], risk.1se)
                           
        } else {
            iauc <- 0.5
            uno.conc <- 0.5
        }

        if (erstatus=='all') {
            uno <- concordance(test[,sf] ~ strata(test$ER_STATUS) + risk.1se,
                               timewt='n/G2')
        
            harrell <- concordance(test[,sf] ~ strata(test$ER_STATUS) + risk.1se)
        } else {
            uno <- concordance(test[,sf] ~ risk.1se,
                               timewt='n/G2')
        
            harrell <- concordance(test[,sf] ~ risk.1se)
        }

        feats <- dimnames(coef(cv.out, s=cv.out$lambda.1se))[[1]][as.vector(coef(cv.out, s=cv.out$lambda.1se)!=0)]

        numfeats <- length(feats)

        for (fact.feat in c('GRADE', 'ONCOTREE_CODE', 'LYMPH_NODE_STATUS')) {
            feat.names <- paste0(fact.feat, unique(test[,fact.feat]))
            if (sum(feat.names %in% feats) > 1) {
                numfeats <- numfeats - sum(feat.names %in% feats) + 1
            }
        }
        
        
        eval.df <- rbind(eval.df,
                         data.frame(Outcome = sf,
                                    Features = '1SE',
                                    Time = 'Full',
                                    Metric = "log Partial Likelihood",
                                    Value = test.logpl[1] / test.logpl[2], ## / nrow(test),
                                    check.names=F),
                         data.frame(Outcome = sf,
                                    Features = '1SE',
                                    Time = 'Full',
                                    Metric = "Uno's Concordance",
                                    Value = 1-uno$concordance, ## uno.conc,
                                    check.names=F),
                         data.frame(Outcome = sf,
                                    Features = '1SE',
                                    Time = 'Full',
                                    Metric = "Harrell's Concordance",
                                    Value = 1-harrell$concordance,
                                    check.names=F),
                         data.frame(Outcome = sf,
                                    Features = '1SE',
                                    Time = 'Full',
                                    Metric = "Integrated C/D AUC",
                                    Value = iauc,
                                    check.names=F),
                         data.frame(Outcome = sf,
                                    Features = '1SE',
                                    Time = 'Full',
                                    Metric = "Number of Features",
                                    Value = numfeats,
                                    check.names=F))

        
        for (tt in times) {

            if (erstatus=='all') {
                uno <- concordance(test[,sf] ~ strata(test$ER_STATUS) + risk.1se,
                                   timewt='n/G2', ymax=365.25*tt)
                
                harrell <- concordance(test[,sf] ~ strata(test$ER_STATUS) + risk.1se,
                                       ymax=365.25*tt)
            } else {
                uno <- concordance(test[,sf] ~ risk.1se,
                                   timewt='n/G2', ymax=365.25*tt)
                
                harrell <- concordance(test[,sf] ~ risk.1se, ymax=365.25*tt)
            }


            if (length(unique(risk.1se))!=1) { 

                ## test.roc <- timeROC(test[,sf][,1], test[,sf][,2],
                ##                     risk.1se, cause=1,
                ##                     times=365.25*tt, weighting='cox')

                ## auc <- test.roc$AUC[2]

                test.auc.uno <- AUC.uno(train[,sf], test[,sf], risk.1se, times=365.25*tt)

                auc <- test.auc.uno$auc

                ## uno.conc <- UnoC(train[,sf], test[,sf], risk.1se, 365.25 * tt)
                                
                
            } else {
                auc <- 0.5
                uno.conc <- 0.5
            }
            
            eval.df <- rbind(eval.df,
                             data.frame(Outcome = sf,
                                        Features = '1SE',
                                        Time = paste(tt, 'Years'),
                                        Metric = "Uno's Concordance",
                                        Value = 1-uno$concordance, ## uno.conc,
                                        check.names=F),
                             data.frame(Outcome = sf,
                                        Features = '1SE',
                                        Time = paste(tt, 'Years'),
                                        Metric = "Harrell's Concordance",
                                        Value = 1-harrell$concordance,
                                        check.names=F),
                             data.frame(Outcome = sf,
                                        Features = '1SE',
                                        Time = paste(tt, 'Years'),
                                        Metric = "C/D AUC",
                                        Value = auc,
                                        check.names=F))
            
        }
    }
    return(eval.df)
}

evaluateCausal <- function(train, test, cvGraphs, erstatus, times) {
    require(survival)
    require(pracma)
    require(survcomp)
    require(survAUC)
    ## require(CPE)
    
    splitStatus <- c('split')
    mgmSelect <- c('BIC', 'StARS')
    orientRules <- c('maxp')
    alphas <- c(0.05, 0.1, 0.2)
    
    if (erstatus!='relapse') {
        surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')
    } else {
        surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'OS', 'DRFS', 'DFS', 'BONE',
                        'BRAIN', 'LIVER', 'LNS', 'LOCOREGIONAL', 'LUNG', 'PLEURA')
    }

    eventCounts <- sapply(surv.feats, function(x) sum(train[,x][,2]) + sum(test[,x][,2]))
    names(eventCounts) <- surv.feats
    
    eval.df <- data.frame()

    for (ss in splitStatus) {
        for (mgms in mgmSelect) {
            for (rule in orientRules) {
                for (alpha in alphas) {
                    g <- cvGraphs[[ss]][[mgms]][[rule]][[paste0('FDR',gsub('0[.]','',alpha))]]
                    totalDev <- 0
                    totalDev.split <- 0
                    totalDev.terminal <- 0
                    ## totalEvents <- 0
                    ## split.events <- 0
                    ## terminal.events <- 0
                    cph.models <- list()
                    for (sf in surv.feats) {
                        cph.models[[sf]] <- list()
                        if (sf %in% c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse') ||
                            ss == 'composite') {
                            nb <- sort(getNeighbors(g, sf))
                            mb <- sort(setdiff(g$markov.blankets[[sf]],
                                               c('ER_STATUS', surv.feats)))
                            if (erstatus=='all' && sf != 'OD') {
                                f.nb <- as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                                         paste(c(1, nb), collapse=' + ')))
                                f.mb <- as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                                         paste(c(1, mb), collapse=' + ')))
                            } else {
                                f.nb <- as.formula(paste(sf, '~ ',
                                                         paste(c(1, nb), collapse=' + ')))
                                f.mb <- as.formula(paste(sf, '~ ',
                                                         paste(c(1, mb), collapse=' + ')))
                            }
                            cph.models[[sf]][['NB']] <- coxph(f.nb, train)
                            cph.models[[sf]][['MB']] <- coxph(f.mb, train)

                            train[,paste('risk','NB',sf,sep='.')] <-
                                predict(cph.models[[sf]][['NB']], newdata=train)
                            test[,paste('risk','NB',sf,sep='.')] <-
                                predict(cph.models[[sf]][['NB']], newdata=test)
                            train[,paste('risk','MB',sf,sep='.')] <-
                                predict(cph.models[[sf]][['MB']], newdata=train)
                            test[,paste('risk','MB',sf,sep='.')] <-
                                predict(cph.models[[sf]][['MB']], newdata=test)
                        } else {
                            sf.splits <- surv.feats[1:(which(surv.feats==sf)-3)]
                            nb.union <- sort(
                                setdiff(
                                    unique(
                                        unlist(
                                            lapply(sf.splits,
                                                   getNeighbors, graph=g))),
                                    c('ER_STATUS', surv.feats)))
                            
                            mb.union <- sort(
                                setdiff(
                                    unique(unlist(g$markov.blankets[sf.splits])),
                                    c('ER_STATUS', surv.feats)))
                            
                            nb <- paste('risk','NB',sf.splits,sep='.')
                            mb <- paste('risk','MB',sf.splits,sep='.')
                            if (erstatus=='all') {
                                f.nb <- as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                                         paste(c(1, nb), collapse=' + ')))
                                f.mb <- as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                                         paste(c(1, mb), collapse=' + ')))
                            } else {
                                f.nb <- as.formula(paste(sf, '~ ',
                                                         paste(c(1, nb), collapse=' + ')))
                                f.mb <- as.formula(paste(sf, '~ ',
                                                         paste(c(1, mb), collapse=' + ')))
                            }
                            cph.models[[sf]][['NB']] <- coxph(f.nb, train)
                            cph.models[[sf]][['MB']] <- coxph(f.mb, train)

                            train[,paste('risk','NB',sf,sep='.')] <-
                                predict(cph.models[[sf]][['NB']], newdata=train)
                            test[,paste('risk','NB',sf,sep='.')] <-
                                predict(cph.models[[sf]][['NB']], newdata=test)
                            train[,paste('risk','MB',sf,sep='.')] <-
                                predict(cph.models[[sf]][['MB']], newdata=train)
                            test[,paste('risk','MB',sf,sep='.')] <-
                                predict(cph.models[[sf]][['MB']], newdata=test)
                        }

                        risk.nb <- predict(cph.models[[sf]][['NB']], newdata=test)

                        if (erstatus=='all') {
                            test.logpl <- logpl(risk.nb,
                                                test[,sf][,1],
                                                test[,sf][,2],
                                                test$ER_STATUS)

                            test.logpl.null <- logpl(rep(0, nrow(test)),
                                                     test[,sf][,1],
                                                     test[,sf][,2],
                                                     test$ER_STATUS)
                        } else {
                            test.logpl <- logpl(risk.nb,
                                                test[,sf][,1],
                                                test[,sf][,2])

                            test.logpl.null <- logpl(rep(0, nrow(test)),
                                                     test[,sf][,1],
                                                     test[,sf][,2])
                        }

                        test.dev <- -2 * test.logpl[1] / test.logpl[2]
                        test.dev.null <- -2 * test.logpl.null[1] / test.logpl.null[2]

                        ## split.events <- split.events + ifelse(sf %in% c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse'), test.logpl[2], 0)
                        ## terminal.events <- terminal.events + ifelse(sf %in% c('DSS', 'OS', 'DRFS', 'DFS'), test.dev, 0)

                        ## totalEvents <- totalEvents + test.logpl[2]
                        totalDev <- totalDev + test.dev
                        
                        totalDev.split <- totalDev.split + ifelse(sf %in% c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse'), test.dev, 0)
                        totalDev.terminal <- totalDev.terminal + ifelse(sf %in% c('DSS', 'OS', 'DRFS', 'DFS'), test.dev, 0)

                        
                        if (length(unique(risk.nb))!=1) {
                            ## test.roc <- timeROC(test[,sf][,1], test[,sf][,2],
                            ##                     risk.nb,
                            ##                     cause=1, times=test[,sf][,1],
                            ##                     weighting='cox')

                            
                            ## iauc <- trapz(test.roc$times[!is.na(test.roc$AUC)],
                            ##               test.roc$AUC[!is.na(test.roc$AUC)]) /
                            ##     (max(test.roc$times[!is.na(test.roc$AUC)]) -
                            ##      min(test.roc$times[!is.na(test.roc$AUC)]))

                            test.auc.uno <- AUC.uno(train[,sf], test[,sf], risk.nb, times=sort(test[,sf][,1]))

                            iauc <- test.auc.uno$iauc

                            uno.conc <- UnoC(train[,sf], test[,sf], risk.nb)
                            
                            ## na.mask <-!is.na(coef(cph.models[[sf]][['NB']]))
                            ## coef.names <- names(coef(cph.models[[sf]][['NB']]))

                            ## cpe <- phcpe2(coef(cph.models[[sf]][['NB']])[na.mask],
                            ##               vcov(cph.models[[sf]][['NB']])[na.mask,na.mask],
                            ##               model.matrix(
                            ##                   ~., test)[,coef.names[na.mask]],
                            ##               CPE.SE=TRUE)

                            ## if (erstatus=='all') {
                            ##     er.frac <- table(test$ER_STATUS)
                                
                            ##     er.frac <- er.frac / sum(er.frac)
                                
                            ##     cpe.erp <- phcpe2(
                            ##         coef(cph.models[[sf]][['NB']])[na.mask],
                            ##         vcov(cph.models[[sf]][['NB']])[na.mask,na.mask],
                            ##         model.matrix(
                            ##             ~., test)[test$ER_STATUS=='Positive',
                            ##                       coef.names[na.mask]],
                            ##         CPE.SE=TRUE)

                            ##     cpe.ern <-
                            ##         phcpe2(coef(cph.models[[sf]][['NB']])[na.mask],
                            ##                vcov(cph.models[[sf]][['NB']])[na.mask,na.mask],
                            ##                model.matrix(
                            ##                    ~., test)[test$ER_STATUS=='Negative',
                            ##                              coef.names[na.mask]],
                            ##                CPE.SE=TRUE)

                            ##     cpe[[1]] <- sum(er.frac[c('Negative', 'Positive')] *
                            ##                     c(cpe.ern[[1]], cpe.erp[[1]]))
                                
                            ## }
                        } else {
                            iauc <- 0.5
                            ## cpe <- list(0.5)
                            uno.conc <- 0.5
                        }                        
                        

                        ## uno <- concordance(cph.models[[sf]][['NB']], newdata=test,
                        ##                    timewt='n/G2')
                        
                        harrell <- concordance(cph.models[[sf]][['NB']], newdata=test)

                        eval.df <- rbind(eval.df,
                                         data.frame(`Split Status` = ss,
                                                    `MGM Select` = mgms,
                                                    `Orient Rule` = rule,
                                                    FDR = alpha,
                                                    Outcome = sf,
                                                    Features = 'NB',
                                                    Time = 'Full',
                                                    Metric = "Deviance",
                                                    Value = test.dev,
                                                    check.names=F),
                                         data.frame(`Split Status` = ss,
                                                    `MGM Select` = mgms,
                                                    `Orient Rule` = rule,
                                                    FDR = alpha,
                                                    Outcome = sf,
                                                    Features = 'NB',
                                                    Time = 'Full',
                                                    Metric = "% Deviance Explained",
                                                    Value = 100 * (test.dev.null - test.dev) / test.dev.null,
                                                    check.names=F),
                                         ## data.frame(`Split Status` = ss,
                                         ##            `MGM Select` = mgms,
                                         ##            `Orient Rule` = rule,
                                         ##            FDR = alpha,
                                         ##            Outcome = sf,
                                         ##            Features = 'NB',
                                         ##            Time = 'Full',
                                         ##            Metric = "CPE",
                                         ##            Value = cpe[[1]],
                                         ##            check.names=F),
                                         data.frame(`Split Status` = ss,
                                                    `MGM Select` = mgms,
                                                    `Orient Rule` = rule,
                                                    FDR = alpha,
                                                    Outcome = sf,
                                                    Features = 'NB',
                                                    Time = 'Full',
                                                    Metric = "Uno's Concordance",
                                                    Value = uno.conc,
                                                    check.names=F),
                                         data.frame(`Split Status` = ss,
                                                    `MGM Select` = mgms,
                                                    `Orient Rule` = rule,
                                                    FDR = alpha,
                                                    Outcome = sf,
                                                    Features = 'NB',
                                                    Time = 'Full',
                                                    Metric = "Harrell's Concordance",
                                                    Value = harrell$concordance,
                                                    check.names=F),
                                         data.frame(`Split Status` = ss,
                                                    `MGM Select` = mgms,
                                                    `Orient Rule` = rule,
                                                    FDR = alpha,
                                                    Outcome = sf,
                                                    Features = 'NB',
                                                    Time = 'Full',
                                                    Metric = "Integrated C/D AUC",
                                                    Value = iauc,
                                                    check.names=F),
                                         data.frame(`Split Status` = ss,
                                                    `MGM Select` = mgms,
                                                    `Orient Rule` = rule,
                                                    FDR = alpha,
                                                    Outcome = sf,
                                                    Features = 'NB',
                                                    Time = 'Full',
                                                    Metric = "Number of Features",
                                                    Value =
                                                        ifelse(
                                                            sf%in%c('OS','DRFS','DFS') &
                                                            ss=='split',
                                                            length(nb.union),
                                                            length(nb)),
                                                    check.names=F))

                        
                        for (tt in times) {
                            ## uno <- concordance(cph.models[[sf]][['NB']], newdata=test,
                            ##                    ymax=365.25*tt, timewt='n/G2')
                            
                            harrell <- concordance(cph.models[[sf]][['NB']], newdata=test,
                                                   ymax=365.25*tt)

                            if (length(unique(risk.nb))!=1) { 

                                ## test.roc <- timeROC(test[,sf][,1], test[,sf][,2],
                                ##                     risk.nb, cause=1,
                                ##                     times=365.25*tt, weighting='cox')

                                ## auc <- test.roc$AUC[2]

                                test.auc.uno <- AUC.uno(train[,sf], test[,sf], risk.nb, times=365.25*tt)

                                auc <- test.auc.uno$auc

                                uno.conc <- UnoC(train[,sf], test[,sf], risk.nb, 365.25 * tt)
                                
                            } else {
                                auc <- 0.5
                                uno.conc <- 0.5
                            }
                            
                            eval.df <- rbind(eval.df,
                                             data.frame(`Split Status` = ss,
                                                        `MGM Select` = mgms,
                                                        `Orient Rule` = rule,
                                                        FDR = alpha,
                                                        Outcome = sf,
                                                        Features = 'NB',
                                                        Time = paste(tt, 'Years'),
                                                        Metric = "Uno's Concordance",
                                                        Value = uno.conc,
                                                        check.names=F),
                                             data.frame(`Split Status` = ss,
                                                        `MGM Select` = mgms,
                                                        `Orient Rule` = rule,
                                                        FDR = alpha,
                                                        Outcome = sf,
                                                        Features = 'NB',
                                                        Time = paste(tt, 'Years'),
                                                        Metric = "Harrell's Concordance",
                                                        Value = harrell$concordance,
                                                        check.names=F),
                                             data.frame(`Split Status` = ss,
                                                        `MGM Select` = mgms,
                                                        `Orient Rule` = rule,
                                                        FDR = alpha,
                                                        Outcome = sf,
                                                        Features = 'NB',
                                                        Time = paste(tt, 'Years'),
                                                        Metric = "C/D AUC",
                                                        Value = auc,
                                                        check.names=F))
                            
                        }

                        
                        ## risk.mb <- predict(cph.models[[sf]][['MB']], newdata=test)

                        ## if (erstatus=='all') {
                        ##     test.logpl <- logpl(risk.mb,
                        ##                         test[,sf][,1],
                        ##                         test[,sf][,2],
                        ##                         test$ER_STATUS)
                        ## } else {
                        ##     test.logpl <- logpl(risk.mb,
                        ##                         test[,sf][,1],
                        ##                         test[,sf][,2])
                        ## }

                        ## test.dev <- -2 * test.logpl[1]
                        
                        ## if (length(unique(risk.mb))!=1) {
                        ##     test.roc <- timeROC(test[,sf][,1], test[,sf][,2],
                        ##                         risk.mb,
                        ##                         cause=1, times=test[,sf][,1],
                        ##                         weighting='cox')
                            
                        ##     iauc <- trapz(test.roc$times[!is.na(test.roc$AUC)],
                        ##                   test.roc$AUC[!is.na(test.roc$AUC)]) /
                        ##         (max(test.roc$times[!is.na(test.roc$AUC)]) -
                        ##          min(test.roc$times[!is.na(test.roc$AUC)]))

                        ##     na.mask <-!is.na(coef(cph.models[[sf]][['MB']]))
                        ##     coef.names <- names(coef(cph.models[[sf]][['MB']]))

                        ##     cpe <- phcpe2(coef(cph.models[[sf]][['MB']])[na.mask],
                        ##                   vcov(cph.models[[sf]][['MB']])[na.mask,na.mask],
                        ##                   model.matrix(
                        ##                       ~., test)[,coef.names[na.mask]],
                        ##                   CPE.SE=TRUE)

                        ##     if (erstatus=='all') {
                        ##         er.frac <- table(test$ER_STATUS)
                                
                        ##         er.frac <- er.frac / sum(er.frac)
                                
                        ##         cpe.erp <- phcpe2(
                        ##             coef(cph.models[[sf]][['MB']])[na.mask],
                        ##             vcov(cph.models[[sf]][['MB']])[na.mask,na.mask],
                        ##             model.matrix(
                        ##                 ~., test)[test$ER_STATUS=='Positive',
                        ##                           coef.names[na.mask]],
                        ##             CPE.SE=TRUE)

                        ##         cpe.ern <-
                        ##             phcpe2(coef(cph.models[[sf]][['MB']])[na.mask],
                        ##                    vcov(cph.models[[sf]][['MB']])[na.mask,na.mask],
                        ##                    model.matrix(
                        ##                        ~., test)[test$ER_STATUS=='Negative',
                        ##                                  coef.names[na.mask]],
                        ##                    CPE.SE=TRUE)

                        ##         cpe[[1]] <- sum(er.frac[c('Negative', 'Positive')] *
                        ##                         c(cpe.ern[[1]], cpe.erp[[1]]))
                                
                        ##     }

                        ## } else {
                        ##     iauc <- 0.5
                        ##     cpe <- list(0.5)
                        ## }

                        
                        ## uno <- concordance(cph.models[[sf]][['MB']], newdata=test,
                        ##                    timewt='n/G2')
                        
                        ## harrell <- concordance(cph.models[[sf]][['MB']], newdata=test)

                        ## eval.df <- rbind(eval.df,
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "log Partial Likelihood",
                        ##                             Value = test.logpl[1] / test.logpl[2],
                        ##                             check.names=F),
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "Likelihood Ratio",
                        ##                             Value = test.dev.null - test.dev,
                        ##                             check.names=F),
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "Deviance",
                        ##                             Value = test.dev,
                        ##                             check.names=F),
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "% Deviance Explained",
                        ##                             Value = 100 * (test.dev.null - test.dev) / test.dev.null,
                        ##                             check.names=F),
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "CPE",
                        ##                             Value = cpe[[1]],
                        ##                             check.names=F),
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "Uno's Concordance",
                        ##                             Value = uno$concordance,
                        ##                             check.names=F),
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "Harrell's Concordance",
                        ##                             Value = harrell$concordance,
                        ##                             check.names=F),
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "Integrated C/D AUC",
                        ##                             Value = iauc,
                        ##                             check.names=F),
                        ##                  data.frame(`Split Status` = ss,
                        ##                             `MGM Select` = mgms,
                        ##                             `Orient Rule` = rule,
                        ##                             FDR = alpha,
                        ##                             Outcome = sf,
                        ##                             Features = 'MB',
                        ##                             Time = 'Full',
                        ##                             Metric = "Number of Features",
                        ##                             Value =
                        ##                                 ifelse(
                        ##                                     sf%in%c('OS','DRFS','DFS') &
                        ##                                     ss=='split',
                        ##                                     length(mb.union),
                        ##                                     length(mb)),
                        ##                             check.names=F))

                        
                        ## for (tt in times) {
                        ##     uno <- concordance(cph.models[[sf]][['MB']], newdata=test,
                        ##                        ymax=365.25*tt, timewt='n/G2')
                            
                        ##     harrell <- concordance(cph.models[[sf]][['MB']], newdata=test,
                        ##                            ymax=365.25*tt)

                        ##     if (length(unique(risk.mb))!=1) { 

                        ##         test.roc <- timeROC(test[,sf][,1], test[,sf][,2],
                        ##                             risk.mb, cause=1,
                        ##                             times=365.25*tt, weighting='cox')

                        ##         auc <- test.roc$AUC[2]
                                
                        ##     } else {
                        ##         auc <- 0.5
                        ##     }
                            
                        ##     eval.df <- rbind(eval.df,
                        ##                      data.frame(`Split Status` = ss,
                        ##                                 `MGM Select` = mgms,
                        ##                                 `Orient Rule` = rule,
                        ##                                 FDR = alpha,
                        ##                                 Outcome = sf,
                        ##                                 Features = 'MB',
                        ##                                 Time = paste(tt, 'Years'),
                        ##                                 Metric = "Uno's Concordance",
                        ##                                 Value = uno$concordance,
                        ##                                 check.names=F),
                        ##                      data.frame(`Split Status` = ss,
                        ##                                 `MGM Select` = mgms,
                        ##                                 `Orient Rule` = rule,
                        ##                                 FDR = alpha,
                        ##                                 Outcome = sf,
                        ##                                 Features = 'MB',
                        ##                                 Time = paste(tt, 'Years'),
                        ##                                 Metric = "Harrell's Concordance",
                        ##                                 Value = harrell$concordance,
                        ##                                 check.names=F),
                        ##                      data.frame(`Split Status` = ss,
                        ##                                 `MGM Select` = mgms,
                        ##                                 `Orient Rule` = rule,
                        ##                                 FDR = alpha,
                        ##                                 Outcome = sf,
                        ##                                 Features = 'MB',
                        ##                                 Time = paste(tt, 'Years'),
                        ##                                 Metric = "C/D AUC",
                        ##                                 Value = auc,
                        ##                                 check.names=F))
                            
                        ## }    
                        
                    }

                    eval.df <- rbind(eval.df,
                                     data.frame(`Split Status` = ss,
                                                `MGM Select` = mgms,
                                                `Orient Rule` = rule,
                                                FDR = alpha,
                                                Outcome = 'All',
                                                Features = 'NB',
                                                Time = 'Full',
                                                Metric = "Deviance",
                                                Value = totalDev,
                                                check.names=F),
                                     data.frame(`Split Status` = ss,
                                                `MGM Select` = mgms,
                                                `Orient Rule` = rule,
                                                FDR = alpha,
                                                Outcome = 'All Split',
                                                Features = 'NB',
                                                Time = 'Full',
                                                Metric = "Deviance",
                                                Value = totalDev.split,
                                                check.names=F),
                                     data.frame(`Split Status` = ss,
                                                `MGM Select` = mgms,
                                                `Orient Rule` = rule,
                                                FDR = alpha,
                                                Outcome = 'All Terminal',
                                                Features = 'NB',
                                                Time = 'Full',
                                                Metric = "Deviance",
                                                Value = totalDev.terminal,
                                                check.names=F))
                }
            }
        }
    }
    
    return(eval.df)
}

evaluateModel <- function(model, test, times, target, form=NULL) {
    eval.df <- data.frame()
    for (t in times) {
        if (c('coxph') %in% class(model)) {
            conc.uno <- concordance(model,
                                newdata=test,
                                ymax=365.25*t,
                                timewt='n/G2')$concordance
            conc.harrell <- concordance(model,
                                newdata=test,
                                ymax=365.25*t)$concordance
            eval.df <- rbind(eval.df, data.frame(`Uno Concordance` = conc.uno,
                                                 `Harrell Concordance` = conc.harrell,
                                                 Time = paste(t, 'Years'),
                                                 Target = target,
                                                 check.names=F))
        }
        
        if (c('glmnet') %in% class(model)) {
            risk <- -predict(model, newx=model.matrix(form, test))
            conc.uno <- concordance(
                test[,target] ~ strata(test$ER_STATUS) + risk,
                ymax=365.25*t, timewt='n/G2')$concordance
            conc.harrell <- concordance(
                test[,target] ~ strata(test$ER_STATUS) + risk,
                ymax=365.25*t)$concordance
            eval.df <- rbind(eval.df, data.frame(`Uno Concordance` = conc.uno,
                                                 `Harrell Concordance` = conc.harrell,
                                                 Time = paste(t, 'Years'),
                                                 Target = target,
                                                 check.names=F))
        }

        if (c('cv.glmnet') %in% class(model)) {
            risk <- -predict(model, newx=model.matrix(~., test[,!colnames(test) %in% c('ER_STATUS', 'DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')])[,-1], s=model$lambda.min)
            conc.uno <- concordance(
                test[,target] ~ strata(test$ER_STATUS) + risk,
                ymax=365.25*t, timewt='n/G2')$concordance
            conc.harrell <- concordance(
                test[,target] ~ strata(test$ER_STATUS) + risk,
                ymax=365.25*t)$concordance

            feats <- dimnames(coef(model, s=model$lambda.min))[[1]][as.vector(coef(model, s=model$lambda.min)!=0)]

            numfeats <- length(feats)

            for (fact.feat in c('GRADE', 'ONCOTREE_CODE', 'LYMPH_NODE_STATUS')) {
                feat.names <- paste0(fact.feat, unique(test[,fact.feat]))
                if (sum(feat.names %in% feats) > 1) {
                    numfeats <- numfeats - sum(feat.names %in% feats) + 1
                }
            }
 
            eval.df <- rbind(eval.df,
                             data.frame(
                                 `Uno Concordance` = conc.uno,
                                 `Harrell Concordance` = conc.harrell,
                                 Time = paste(t, 'Years'),
                                 Target = target,
                                 Model = 'LASSO Min',
                                 `Number of Features`=numfeats,
                                 check.names=F))

            
            risk <- -predict(model, newx=model.matrix(~., test[,!colnames(test) %in% c('ER_STATUS', 'DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')])[,-1], s=model$lambda.1se)

            feats <- dimnames(coef(model, s=model$lambda.1se))[[1]][as.vector(coef(model, s=model$lambda.1se)!=0)]

            numfeats <- length(feats)

            for (fact.feat in c('GRADE', 'ONCOTREE_CODE', 'LYMPH_NODE_STATUS')) {
                feat.names <- paste0(fact.feat, unique(test[,fact.feat]))
                if (sum(feat.names %in% feats) > 1) {
                    numfeats <- numfeats - sum(feat.names %in% feats) + 1
                }
            }
            
            conc.uno <- concordance(
                test[,target] ~ strata(test$ER_STATUS) + risk,
                ymax=365.25*t, timewt='n/G2')$concordance
            conc.harrell <- concordance(
                test[,target] ~ strata(test$ER_STATUS) + risk,
                ymax=365.25*t)$concordance
            eval.df <- rbind(eval.df,
                             data.frame(
                                 `Uno Concordance` = conc.uno,
                                 `Harrell Concordance` = conc.harrell,
                                 Time = paste(t, 'Years'),
                                 Target = target,
                                 Model = 'LASSO 1SE',
                                 `Number of Features`=numfeats,
                                 check.names=F))
        }
        
    }

    return(eval.df)
}

tmat <- transMat(list(c(2,3,4,5), c(3,4,5), c(4, 5), c(), c()),
                 c('Surgery', 'LocalRelapse', 'DistantRelapse', 'DSS', 'OD'))

prepDataMSM <- function(train, test, graphs, tmat) {
    data <- rbind.data.frame(train, test)

    msm.surv.feats <- colnames(tmat)[-1]
    for (sf.idx in 1:4) {
        sf <- msm.surv.feats[sf.idx]
        data[,paste(sf, 'status', sep='.')] <- data[,sf][,2]
        for (sf.jdx in (sf.idx+1):4) {
            sf2 <- msm.surv.feats[sf.jdx]
            if (is.na(sf2)) break
            equal.vec <- apply(data[,sf][,1:2]==data[,sf2][,1:2], 1, all)
            data[equal.vec & data[,sf][,2],sf2][,1] <- data[equal.vec & data[,sf][,2],sf2][,1] + 1
        }
        data[,paste(sf, 'time', sep='.')] <- data[,sf][,1]
    }


    erp.covs <- unique(unlist(graphs[['erp']]$neighbors[msm.surv.feats]))

    ern.covs <- unique(unlist(graphs[['ern']]$neighbors[msm.surv.feats]))

    full.covs <- unique(c('ER_STATUS', erp.covs, ern.covs))

    msdata <- msprep(time = c(NA, paste0(msm.surv.feats, '.time')),
                     status = c(NA, paste0(msm.surv.feats, '.status')),
                     data = data, trans = tmat, keep = full.covs,
                     id = rownames(data))

    msdata$Destination <- factor(msdata$to, levels=2:5, labels=c('LocalRelapse', 'DistantRelapse', 'DSS', 'OD'))


    f <- as.formula(paste('~0 + ', paste(paste0('ER_STATUS:Destination:(',paste(setdiff(full.covs, 'ER_STATUS'), collapse=' + '),')'))))

    expanded <- model.matrix(f, msdata)

    head(expanded)

    msdata <- cbind(msdata[,!colnames(msdata) %in% c(erp.covs, ern.covs)], data.frame(expanded))

    head(msdata)

    msdata$strata <- factor(paste(msdata$ER_STATUS, msdata$Destination, sep='.'),
                                 levels=c(paste('Positive', msm.surv.feats, sep='.'),
                                          paste('Negative', msm.surv.feats, sep='.')),
                                 labels=1:8)

    mstrain <- msdata[msdata$id %in% rownames(train),]
    mstest <- msdata[msdata$id %in% rownames(test),]

    return(list(train=mstrain, test=mstest))
}


msmTrain <- function(train, mstrain, graphs, tmat) {
    msm.surv.feats <- colnames(tmat)[-1]

    terms <- c()

    for (erstatus in c('erp', 'ern')) {
        if (erstatus=='erp') erstring <- 'ER_STATUSPositive'
        else if (erstatus=='ern') erstring <- 'ER_STATUSNegative'
        for (sf in msm.surv.feats) {
            tostring <- paste0('Destination',sf)
            for (feat in graphs[[erstatus]]$neighbors[[sf]]) {
                featstring <- paste0(feat, levels(train[,feat])[-1])
                terms <- c(terms, paste(erstring, tostring, featstring, sep='.'))
            }
        }
    }

    f <- as.formula(paste('Surv(Tstart, Tstop, status) ~ strata(strata) + ', paste(terms, collapse=' + ')))

    model <- coxph(f, mstrain %>% filter(time>0))

    print(summary(model))

    return(model)
    
}

patient.chf <- function(model, data, tmat) {
    msf1 <- msfit(model, newdata = data, trans = tmat)
    pt <- probtrans(msf1, predt=0)

    nTime <- nrow(na.omit(pt[[1]]))

    chf.df <- data.frame(Time=pt[[1]][1:nTime,'time'],
                         DSS=-log(1-pt[[1]][1:nTime,'pstate4']),
                         OS=-log(1-pt[[1]][1:nTime,'pstate4']-pt[[1]][1:nTime,'pstate5']),
                         DRFS=-log(1-pt[[1]][1:nTime,'pstate3']-pt[[1]][1:nTime,'pstate4']-pt[[1]][1:nTime,'pstate5']),
                         DFS=-log(pt[[1]][1:nTime,'pstate1']))

    chf.fun <- list()
    for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
        chf.fun[[sf]] <- approxfun(chf.df[,'Time'], chf.df[,sf], method='constant')
    }
    return(list(chf.df=chf.df,
                chf.fun=chf.fun))    
}


msmPredict <- function(model, mstest, tmat, times=NULL) {
    require(parallel)
    
    samps <- unique(mstest$id)

    riskdf <- do.call(rbind.data.frame,
                      mclapply(samps, function(samp) {
                          riskdf <- data.frame()
                          chf.out <- patient.chf(model, mstest[mstest$id==samp,], tmat)
                          single <- cbind.data.frame(id=samp,
                                                     ER_STATUS=mstest$ER_STATUS[1],
                                                     tail(chf.out$chf.df,1))
                          single$Time <- 'Full'
                          riskdf <- rbind(riskdf, single)
                          for (tt in times) {
                              single$Time <- paste(tt,'Years')
                              single$DSS  <- chf.out$chf.fun$DSS(365.25*tt)
                              single$OS   <- chf.out$chf.fun$OS(365.25*tt)
                              single$DRFS <- chf.out$chf.fun$DRFS(365.25*tt)
                              single$DFS  <- chf.out$chf.fun$DFS(365.25*tt)
                              riskdf <- rbind(riskdf, single)
                          }
                          return(riskdf)
                      }, mc.cores=detectCores()))

    return(riskdf)
}


library(mstate)

prepDataMSM <- function(train, test, graphs, tmat) {
    data <- rbind.data.frame(train, test)

    msm.surv.feats <- colnames(tmat)[-1]
    for (sf.idx in 1:4) {
        sf <- msm.surv.feats[sf.idx]
        data[,paste(sf, 'status', sep='.')] <- data[,sf][,2]
        for (sf.jdx in (sf.idx+1):4) {
            sf2 <- msm.surv.feats[sf.jdx]
            if (is.na(sf2)) break
            equal.vec <- apply(data[,sf][,1:2]==data[,sf2][,1:2], 1, all)
            data[equal.vec & data[,sf][,2],sf2][,1] <- data[equal.vec & data[,sf][,2],sf2][,1] + 1
        }
        data[,paste(sf, 'time', sep='.')] <- data[,sf][,1]
    }


    erp.covs <- unique(unlist(graphs[['erp']]$neighbors[msm.surv.feats]))

    ern.covs <- unique(unlist(graphs[['ern']]$neighbors[msm.surv.feats]))

    full.covs <- unique(c('ER_STATUS', erp.covs, ern.covs))

    msdata <- msprep(time = c(NA, paste0(msm.surv.feats, '.time')),
                     status = c(NA, paste0(msm.surv.feats, '.status')),
                     data = data, trans = tmat, keep = full.covs,
                     id = rownames(data))

    msdata$Destination <- factor(msdata$to, levels=2:5, labels=c('LocalRelapse', 'DistantRelapse', 'DSS', 'OD'))


    f <- as.formula(paste('~0 + ', paste(paste0('ER_STATUS:Destination:(',paste(setdiff(full.covs, 'ER_STATUS'), collapse=' + '),')'))))

    expanded <- model.matrix(f, msdata)

    head(expanded)

    msdata <- cbind(msdata[,!colnames(msdata) %in% c(erp.covs, ern.covs)], data.frame(expanded))

    head(msdata)

    msdata$strata <- factor(paste(msdata$ER_STATUS, msdata$Destination, sep='.'),
                                 levels=c(paste('Positive', msm.surv.feats, sep='.'),
                                          paste('Negative', msm.surv.feats, sep='.')),
                                 labels=1:8)

    mstrain <- msdata[msdata$id %in% rownames(train),]
    mstest <- msdata[msdata$id %in% rownames(test),]

    return(list(train=mstrain, test=mstest))
}


msmTrain <- function(train, mstrain, graphs, tmat) {
    msm.surv.feats <- colnames(tmat)[-1]

    terms <- c()

    for (erstatus in c('erp', 'ern')) {
        if (erstatus=='erp') erstring <- 'ER_STATUSPositive'
        else if (erstatus=='ern') erstring <- 'ER_STATUSNegative'
        for (sf in msm.surv.feats) {
            tostring <- paste0('Destination',sf)
            for (feat in graphs[[erstatus]]$neighbors[[sf]]) {
                featstring <- paste0(feat, levels(train[,feat])[-1])
                terms <- c(terms, paste(erstring, tostring, featstring, sep='.'))
            }
        }
    }

    f <- as.formula(paste('Surv(Tstart, Tstop, status) ~ strata(strata) + ', paste(terms, collapse=' + ')))

    model <- coxph(f, mstrain %>% filter(time>0))

    print(summary(model))

    return(model)
    
}

patient.chf <- function(model, data, tmat) {
    msf1 <- msfit(model, newdata = data, trans = tmat)
    pt <- probtrans(msf1, predt=0)

    nTime <- nrow(na.omit(pt[[1]]))

    chf.df <- data.frame(Time=pt[[1]][1:nTime,'time'],
                         DSS=-log(1-pt[[1]][1:nTime,'pstate4']),
                         OS=-log(1-pt[[1]][1:nTime,'pstate4']-pt[[1]][1:nTime,'pstate5']),
                         DRFS=-log(1-pt[[1]][1:nTime,'pstate3']-pt[[1]][1:nTime,'pstate4']-pt[[1]][1:nTime,'pstate5']),
                         DFS=-log(pt[[1]][1:nTime,'pstate1']))

    chf.fun <- list()
    for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
        chf.fun[[sf]] <- approxfun(chf.df[,'Time'], chf.df[,sf], method='constant')
    }
    return(list(chf.df=chf.df,
                chf.fun=chf.fun))    
}


msmPredict <- function(model, mstest, tmat, times=NULL) {
    require(parallel)
    
    samps <- unique(mstest$id)

    riskdf <- do.call(rbind.data.frame,
                      mclapply(samps, function(samp) {
                          riskdf <- data.frame()
                          chf.out <- patient.chf(model, mstest[mstest$id==samp,], tmat)
                          single <- cbind.data.frame(id=samp,
                                                     ER_STATUS=mstest$ER_STATUS[1],
                                                     tail(chf.out$chf.df,1))
                          single$Time <- 'Full'
                          riskdf <- rbind(riskdf, single)
                          for (tt in times) {
                              single$Time <- paste(tt,'Years')
                              single$DSS  <- chf.out$chf.fun$DSS(365.25*tt)
                              single$OS   <- chf.out$chf.fun$OS(365.25*tt)
                              single$DRFS <- chf.out$chf.fun$DRFS(365.25*tt)
                              single$DFS  <- chf.out$chf.fun$DFS(365.25*tt)
                              riskdf <- rbind(riskdf, single)
                          }
                          return(riskdf)
                      }, mc.cores=detectCores()))

    return(riskdf)
}


cvData.list <- list()
cvData.list[['all']] <- lapply(1:10,
                               function(k) {
                                   cvData <- loadData('all', k)
                                   return(cvData)
                               })

cvData.list[['ern']] <- lapply(1:10,
                               function(k) {
                                   cvData <- loadData('ern', k)
                                   return(cvData)
                               })

cvData.list[['erp']] <- lapply(1:10,
                               function(k) {
                                   cvData <- loadData('erp', k)
                                   return(cvData)
                               })

fullData.list <- lapply(c('all', 'ern', 'erp'), loadData, k=-1)
names(fullData.list) <- c('all', 'ern', 'erp')


xmat <- model.matrix(~., fullData.list$ern$train[,!colnames(fullData.list$ern$train) %in% surv.feats])[,-1:-11]


localVar <- function(x, k=6) {
    n <- length(x)
    ordX <- sort(x)
    low <- 1
    high <- k
    localVar <- 0
    localMu <- mean(ordX[low:high])
    for (i in 1:n) {
        while (high < n & ordX[i]-ordX[low] > ordX[high]-ordX[i]) {
            localMu <- localMu - ordX[low]/k
            low <- low+1
            high <- high+1
            localMu <- localMu + ordX[high]/k
        }
        localVar <- localVar + sum((ordX[low:high] - localMu)^2)
    }
    localVar / (n-1)
}

lgdScore <- function(x, k=6) {
    var(x) / localVar(x, k)
}

grmSelect <- function(A, scores, rho=1, mu=1) {
    p <- ncol(A)
    z <- matrix(1/p, p)
    s <- matrix(scores / sum(scores))
    alpha <- matrix(0, p)
    centerMat <- matrix(1, p) %*% matrix(1, 1, p) / p

    oldLambda <- Inf
    lambda <- as.numeric(t(z) %*% A %*% z / t(z) %*% s)

    idx <- 1

    print(paste("Iter ", idx, ": lambda =", lambda))

    while ((oldLambda - lambda)/lambda > 1e-8) {
        oldLambda <- lambda

        v <- solve(A + diag(mu/2, p, p), lambda / 2 * s + mu / 2 * (z + alpha/mu))
        w <- v - alpha / mu
        u <- w - centerMat %*% w + 1 / p

        lambdaBar <- 0
        f <- 1

        while (f > 1e-10) {
            f <- mean(pmax(lambdaBar - u, 0)) - lambdaBar
            posFreq <- mean(lambdaBar - u > 0)
            fprime <- posFreq - 1

            lambdaBar <- lambdaBar - f / fprime
        }

        z <- pmax(u - lambdaBar, 0)

        lambda <- as.numeric(t(z) %*% A %*% z / t(z) %*% s)
        idx <- idx + 1

        alpha <- alpha + mu * (z - v)

        mu <- rho * mu

        print(paste("Iter ", idx, ": lambda =", lambda))
    }
    grmScore <- as.vector(z)
    names(grmScore) <- rownames(z)
    grmScore
}

scores <- apply(xmat, 2, lgdScore, k=6)

hist(scores)

hist(sqrt(scores))

sort(scores)

colnames(xmat)

corMat <- cor(xmat)

grmScores <- grmSelect(corMat^2, scores=sqrt(scores), rho=1.01, mu=1)

pheatmap::pheatmap(corMat, breaks=seq(-1,1,length.out=101))

grmvars3 <- rev(names(sort(grmScores))[sort(grmScores)>0])

grmvars

setdiff(grmvars, grmvars2)

pheatmap::pheatmap(corMat[grmvars3[1:150],grmvars3[1:150]], breaks=seq(-1,1,length.out=101))

(sum(corMat[grmvars3,grmvars3]^2)-length(grmvars3))/(length(grmvars3) * (length(grmvars3)-1))

pheatmap::pheatmap(corMat[grmvars2[1:150],grmvars2[1:150]], breaks=seq(-1,1,length.out=101))

(sum(corMat[grmvars2,grmvars2]^2)-length(grmvars2))/(length(grmvars2) * (length(grmvars2)-1))

redvar <- rep(NA, length(grmvars3))
redsd <- rep(NA, length(grmvars3))

for (k in 2:length(grmvars3)) {
    sdfeats <- grmvars3[1:k]
    redsd[k] <- (sum(corMat[sdfeats,sdfeats]^2)-k)/(k * (k-1))
    if (k <= length(grmvars2)) {
        varfeats <- grmvars2[1:k]
        redvar[k] <- (sum(corMat[varfeats,varfeats]^2)-k)/(k * (k-1))
    }
}

plot(redsd)
points(redvar, col='red')

g.boss.grmern <- boss(xmat[,grmvars], rank=T, verbose=T)

g.boss.grmern
plot(g.boss.grmern, "PIP")

g.grasp.grmern <- grasp(xmat[,grmvars], depth=3, rank=T, verbose=T)

plot(g.grasp.grmern, "PIP")


hist(sapply(g.boss.grmern$neighbors, length

sum(grmScores)

hist(grmScores)

sort(grmScores)

sort(colMeans(abs(corMat) - diag(1, nrow=nrow(corMat))))

p <- ncol(xmat)

parCorMat <- matrix(0, p, p)
colnames(parCorMat) <- colnames(xmat)
rownames(parCorMat) <- colnames(xmat)

for (j in 1:p) {
    res <- glmnet(xmat[,-1*j], y=xmat[,j], family='gaussian', alpha=0, lambda=0.125)
    beta <- coef(res)[-1,]
    parCorMat[j,names(beta)] <- beta
}

pheatmap::pheatmap(parCorMat)

parCorMat <- (parCorMat + t(parCorMat)) / 2

pheatmap::pheatmap(parCorMat)

relevance <- colSums(abs(parCorMat)) / (p-1)

relevance <- rep(0, p)
names(relevance) <- colnames(xmat)

for (j in 1:p) {
    res <- glmnet(xmat[,-1*j], y=xmat[,j], family='gaussian', alpha=0, lambda=0.5)
    beta <- coef(res)[-1,]
    relevance[names(beta)] <- relevance[names(beta)] + abs(beta) / (p-1)
}

corMat <- cor(xmat)

library(rags2ridges)

res <- optPenalty.kCV(xmat, 0.01, 0.1, fold=5, step=50)

ridgePrec <- res$optPrec
diag(ridgePrec) <- 1

ridgeCor <- solve(res$optPrec)
diag(ridgeCor) <- 0

hist(ridgePrec)

relevancea <- colSums(abs(ridgePrec)) / (p-1)

relevanceg <- exp(colSums(log(abs(ridgePrec))) / (p-1))

relevance <- relevancea

relevance <- colSums(abs(ridgeCor)) / (p-1)

relevance <- colSums(abs(corMat) - diag(1, nrow=nrow(corMat))) / (p-1)

relevance <- exp(colSums(log(abs(corMat))) / (p-1))


K <- 100
S <- c()
while (length(S) < K) {
    remaining <- setdiff(colnames(xmat), S)

    ## res <- optPenalty.kCV(xmat[,remaining], 0.01, 0.1, fold=5, step=50)

    ## ridgePrec <- res$optPrec
    ## diag(ridgePrec) <- 0
    
    r <- length(remaining)
    ## mrmr <- colSums(abs(ridgePrec)) / (r-1)
    ## mrmr <- rep(0, r)
    ## names(mrmr) <- remaining

    ## for (j in 1:r) {
    ##     res <- glmnet(xmat[,remaining[-1*j]], y=xmat[,remaining[j]], family='gaussian', alpha=0, lambda=0.5)
    ##     beta <- coef(res)[-1,]
    ##     mrmr[names(beta)] <- relevance[names(beta)] + abs(beta) / (r-1)
    ## }

    mrmr <- relevance[remaining]

    if (length(S) > 0) {
        for (col in remaining) {
            mrmr[col] <- mrmr[col] / mean(abs(ridgeCor[col,S])) ## exp(mean(log(abs(corMat[col,S])))) ##
        }
    }
    S <- c(S, names(mrmr)[which.max(mrmr)])
}

pheatmap::pheatmap(corMat[S,S], breaks=seq(-1,1,length.out=101))

max(abs(corMat[S,S])-diag(1,K))

sum(abs(corMat[S,S])-diag(1,K)) / (K * (K-1))

res <- iCoxBoost(fullData.list$erp$train.npn$DSS ~ .,
                 data=as.data.frame(xmat),
                 mandatory=colnames(xmat)[grepl("TUMOR_SIZE", colnames(xmat)) |
                                          grepl("GRADE", colnames(xmat)) |
                                          grepl("LYMPH_NODE_STATUS", colnames(xmat))])

res

summary(res)

plot(res)


for (k in 1:10) {
    print(paste("Check for train/test overlap for fold", k, "ER-:",
                any(rownames(cvData.list$all[[k]]$train) %in% rownames(cvData.list$ern[[k]]$test))))
        print(paste("Check for train/test overlap for fold", k, "ER+:",
                any(rownames(cvData.list$all[[k]]$train) %in% rownames(cvData.list$erp[[k]]$test))))
}


tmat <- transMat(list(c(2,3,4,5), c(3,4,5), c(4, 5), c(), c()),
                 c('Surgery', 'LocalRelapse', 'DistantRelapse', 'DSS', 'OD'))

paths(tmat)


msmconcdf <- data.frame()

for (k in 1:10) {

    ernGraphs <- loadGraphs('ern', k)
    erpGraphs <- loadGraphs('erp', k)

    msdata <- prepDataMSM(cvData.list[['all']][[k]]$train.npn,
                          cvData.list[['all']][[k]]$test.npn,
                          list(ern=ernGraphs$split$StARS$maxp$FDR1,
                               erp=erpGraphs$split$StARS$maxp$FDR1),
                          tmat)

    msmodel <- msmTrain(cvData.list[['all']][[k]]$train.npn, msdata$train,
                        list(ern=ernGraphs$split$StARS$maxp$FDR1,
                             erp=erpGraphs$split$StARS$maxp$FDR1), tmat)

    mspred <- msmPredict(msmodel, msdata$test, tmat, times=c(3, 5, 8, 10, 15))

    for (tt in c(3, 5, 8, 10, 15)) {

        temp <- mspred %>% filter(Time==paste(tt,'Years')) %>% mutate_if(is.numeric, function(x) -1 * x)
        rownames(temp) <- as.character(temp$id)

        temp <- cbind.data.frame(cvData.list[['all']][[k]]$test.npn[as.character(unique(mspred$id)),], risk=temp)

        for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {

            f <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + risk.', sf))
            
            conc.harrell <- concordance(f, temp, ymax=365.25*tt)
            ## conc.uno <- concordance(f, temp, ymax=365.25*tt, timewt='n/G2')

            single <- data.frame(Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='Stratified',
                                 Fold=k,
                                 `Concordance` = conc.harrell$concordance,
                                 check.names=F)

            msmconcdf <- rbind(msmconcdf, single)

            f <- as.formula(paste0(sf, ' ~ risk.', sf))

            conc.harrell <- concordance(f, temp, ymax=365.25*tt)

            single <- data.frame(Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='All',
                                 Fold=k,
                                 `Concordance` = conc.harrell$concordance,
                                 check.names=F)

            msmconcdf <- rbind(msmconcdf, single)

            conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Positive'), ymax=365.25*tt)
            ## conc.uno <- concordance(f, temp %>% filter(ER_STATUS=='Positive'), ymax=365.25*tt, timewt='n/G2')

            single <- data.frame(Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='Positive',
                                 Fold=k,
                                 `Concordance` = conc.harrell$concordance,
                                 check.names=F)

            msmconcdf <- rbind(msmconcdf, single)

            conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Negative'), ymax=365.25*tt)
            ## conc.uno <- concordance(f, temp %>% filter(ER_STATUS=='Negative'), ymax=365.25*tt, timewt='n/G2')

            single <- data.frame(Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='Negative',
                                 Fold=k,
                                 `Concordance` = conc.harrell$concordance,
                                 check.names=F)

            msmconcdf <- rbind(msmconcdf, single)

        }
    }
}

write.csv(msmconcdf, 'cv_msm_conc.csv')

head(msmconcdf)

msmconcdf.sum  <-  msmconcdf %>% group_by (`ER Status`, Outcome, Time) %>% summarize_at('Concordance', c(`Concordance SE`=sd, `Concordance`=mean))

head(msmconcdf.sum)

msmconcdf.sum$Time <- factor(msmconcdf.sum$Time, paste(c(3,5,8,10,15), 'Years'))

write.csv(msmconcdf.sum, 'cv_msm_conc_sum.csv')

ggplot(msmconcdf.sum, ## %>% filter(`ER Status`!='All'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=`ER Status`)) +
    geom_pointrange(position=position_dodge(width=0.8)) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome)) +
    theme_bw()


surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')


library(randomForestSRC)

options(rf.cores=16, mc.cores=16)

pdf('cvrf_model_train_plots.pdf', width=8, height=8)
## rfModels <- list()
## rfModels[[sf]] <- list()
for (erstatus in c('all')) {

    print(paste('ER Status:', erstatus))

    for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {

        ## rfModels[[sf]][[erstatus]] <- list()

        print(paste('  Outcome:', sf))

        rfModels <- list()
        for (k in 1:10) {
            print(paste('    Fold:', k))
            
            temp <- cvData.list[[erstatus]][[k]]$train.npn[,setdiff(colnames(cvData.list[[erstatus]][[k]]$train), surv.feats)]
            temp$time <- cvData.list[[erstatus]][[k]]$train.npn[,sf][,1]
            temp$status <- cvData.list[[erstatus]][[k]]$train.npn[,sf][,2]

            rfModels[[k]] <- rfsrc(Surv(time, status) ~ ., temp,
                                   ntree=2000, block.size=10)
            
            print(rfModels[[k]])
            plot(rfModels[[k]])

            gc()
        }
        
        saveRDS(rfModels, paste0('rsf_cv_',erstatus,'_',sf,'_models.rds'))
    }
}
dev.off()

## saveRDS(rfModels, 'rsf_cv_models.rds')


rfconcdf <- data.frame()

for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
    
    print(paste('  Outcome:', sf))
    
    rfModels <- readRDS(paste0('rsf_cv_all_',sf,'_models.rds'))

    for (k in 1:10) {
        print(paste('    Fold:', k))

        pred <- predict(rfModels[[k]], newdata=cvData.list[['all']][[k]]$test.npn[,setdiff(colnames(cvData.list[['all']][[k]]$test), surv.feats)])
        
        temp <- cvData.list[['all']][[k]]$test.npn
        temp[,paste0('risk.', sf)] <- -pred$predicted

        for (tt in c(3, 5, 8, 10, 15)) {
            
            f <- as.formula(paste0(sf, ' ~ risk.', sf)) 

            conc.harrell <- concordance(f, temp, ymax=365.25*tt)

            single <- data.frame(Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='All',
                                 Fold=k,
                                 `Concordance` = conc.harrell$concordance,
                                 check.names=F)

            rfconcdf <- rbind(rfconcdf, single)


            f <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + risk.', sf))

            conc.harrell <- concordance(f, temp, ymax=365.25*tt)

            single <- data.frame(Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='Stratified',
                                 Fold=k,
                                 `Concordance` = conc.harrell$concordance,
                                 check.names=F)

            rfconcdf <- rbind(rfconcdf, single)

            conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Positive'), ymax=365.25*tt)

            single <- data.frame(Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='Positive',
                                 Fold=k,
                                 `Concordance` = conc.harrell$concordance,
                                 check.names=F)

            rfconcdf <- rbind(rfconcdf, single)
            
            conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Negative'), ymax=365.25*tt)

            single <- data.frame(Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='Negative',
                                 Fold=k,
                                 `Concordance` = conc.harrell$concordance,
                                 check.names=F)

            rfconcdf <- rbind(rfconcdf, single)
        }
    }
}


write.csv(rfconcdf, 'cv_rf_conc.csv')

head(rfconcdf)

rfconcdf.sum  <-  rfconcdf %>% group_by (`ER Status`, Outcome, Time) %>% summarize_at('Concordance', c(`Concordance SE`=sd, `Concordance`=mean))

head(rfconcdf.sum)

rfconcdf.sum$Time <- factor(rfconcdf.sum$Time, paste(c(3,5,8,10,15), 'Years'))

write.csv(rfconcdf.sum, 'cv_rf_conc_sum.csv')

ggplot(rfconcdf.sum, ## %>% filter(`ER Status`!='All'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=`ER Status`)) +
    geom_pointrange(position=position_dodge(width=0.8)) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome)) +
    theme_bw()


library(glmnet)

pdf('cvcoxlasso_model_train_plots.pdf', width=8, height=8)
lassoModels <- list()
for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
    lassoModels[[sf]] <- list()
    lassoModels[[sf]][['nostratify']] <- list()

    for (erstatus in c('all', 'ern', 'erp')) {
        lassoModels[[sf]][[erstatus]] <- list()

        for (k in 1:10) {

            temp <- cvData.list[[erstatus]][[k]]$train.npn[,setdiff(colnames(cvData.list[[erstatus]][[k]]$train), c('ER_STATUS', surv.feats))]
    
            lassoModels[[sf]][[erstatus]][[k]] <- cv.glmnet(
                x=model.matrix(~., temp)[,-1],
                y=cvData.list[[erstatus]][[k]]$train.npn[,sf],
                family='cox', lambda.min.ratio=0.01)
    
            print(lassoModels[[sf]][[erstatus]][[k]])
            plot(lassoModels[[sf]][[erstatus]][[k]], main=paste('LASSO',sf,'ER Status:',erstatus,'Fold:',k))

            if (erstatus=='all') {

                temp <- cvData.list[[erstatus]][[k]]$train.npn[,setdiff(colnames(cvData.list[[erstatus]][[k]]$train), surv.feats)]

                lassoModels[[sf]][['nostratify']][[k]] <- cv.glmnet(
                    x=model.matrix(~., temp)[,-1],
                    y=Surv(cvData.list[[erstatus]][[k]]$train.npn[,sf][,1],
                           cvData.list[[erstatus]][[k]]$train.npn[,sf][,2]),
                    family='cox', lambda.min.ratio=0.01)
    
                print(lassoModels[[sf]][['nostratify']][[k]])
                plot(lassoModels[[sf]][['nostratify']][[k]], main=paste('LASSO',sf,'No Stratify, Fold:',k))


            }
        }
    }
}
dev.off()

saveRDS(lassoModels, 'coxlasso_cv_models.rds')

lassoModels <- readRDS('coxlasso_cv_models.rds')

lassoconcdf <- data.frame()

for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
    
    print(paste('  Outcome:', sf))
    
    for (k in 1:10) {
        print(paste('    Fold:', k))

        for (erstatus in c('all', 'ern', 'erp')) {

            temp <- cvData.list[[erstatus]][[k]]$test.npn[,setdiff(colnames(cvData.list[[erstatus]][[k]]$train), c('ER_STATUS', surv.feats))]
        
            pred <- predict(lassoModels[[sf]][[erstatus]][[k]],
                            newx=as.matrix(model.matrix(~., temp)[,-1]),
                            s=lassoModels[[sf]][[erstatus]][[k]]$lambda.min)
        
            temp <- cvData.list[[erstatus]][[k]]$test.npn
            temp[,paste0('risk.', sf)] <- -pred

            if (erstatus != 'all') {
                for (tt in c(3, 5, 8, 10, 15)) {
            
                    f <- as.formula(paste0(sf, ' ~ risk.', sf)) 

                    conc.harrell <- concordance(f, temp, ymax=365.25*tt)

                    single <- data.frame(Outcome=sf,
                                         Time=paste(tt,'Years'),
                                         `ER Status`=ifelse(erstatus=='ern', 'Negative', 'Positive'),
                                         Fold=k,
                                         `Concordance` = conc.harrell$concordance,
                                         check.names=F)

                    lassoconcdf <- rbind(lassoconcdf, single)
                }
            } else {
                for (tt in c(3, 5, 8, 10, 15)) {
                    f <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + risk.', sf))

                    conc.harrell <- concordance(f, temp, ymax=365.25*tt)

                    single <- data.frame(Outcome=sf,
                                         Time=paste(tt,'Years'),
                                         `ER Status`='Stratified',
                                         Fold=k,
                                         `Concordance` = conc.harrell$concordance,
                                         check.names=F)

                    lassoconcdf <- rbind(lassoconcdf, single)
                }

                temp <- cvData.list[[erstatus]][[k]]$test.npn[,setdiff(colnames(cvData.list[[erstatus]][[k]]$train), surv.feats)]
        
                pred <- predict(lassoModels[[sf]][['nostratify']][[k]],
                                newx=as.matrix(model.matrix(~., temp)[,-1]),
                                s=lassoModels[[sf]][[erstatus]][[k]]$lambda.min)
        
                temp <- cvData.list[[erstatus]][[k]]$test.npn
                temp[,paste0('risk.', sf)] <- -pred

                for (tt in c(3, 5, 8, 10, 15)) {
                    f <- as.formula(paste0(sf, ' ~ risk.', sf))

                    conc.harrell <- concordance(f, temp, ymax=365.25*tt)

                    single <- data.frame(Outcome=sf,
                                         Time=paste(tt,'Years'),
                                         `ER Status`='All',
                                         Fold=k,
                                         `Concordance` = conc.harrell$concordance,
                                         check.names=F)

                    lassoconcdf <- rbind(lassoconcdf, single)
                }
            }
        }
    }
}


write.csv(lassoconcdf, 'cv_lasso_conc.csv')

head(lassoconcdf)

lassoconcdf.sum  <-  lassoconcdf %>% group_by (`ER Status`, Outcome, Time) %>% summarize_at('Concordance', c(`Concordance SE`=sd, `Concordance`=mean))

head(lassoconcdf.sum)

lassoconcdf.sum$Time <- factor(lassoconcdf.sum$Time, paste(c(3,5,8,10,15), 'Years'))

write.csv(lassoconcdf.sum, 'cv_lasso_conc_sum.csv')

ggplot(lassoconcdf.sum, ## %>% filter(`ER Status`!='All'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=`ER Status`)) +
    geom_pointrange(position=position_dodge(width=0.8)) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome)) +
    theme_bw()





cv.conc.all <- list()
for (sf in c('DSS','OD','DistantRelapse','LocalRelapse','OS','DRFS','DFS')) {
    cv.conc.all[[sf]] <- list()
    for (t in c('2', '5', '10', '15')) {
        cv.conc.all[[sf]][[t]] <- matrix(NA, 5, 8)
    }
}

cvData.list <- list()
cvData.list[['all']] <- lapply(1:10,
                               function(k) {
                                   cvData <- loadData('all', k)
                                   return(cvData)
                               })

cvData.list[['ern']] <- lapply(1:10,
                               function(k) {
                                   cvData <- loadData('ern', k)
                                   return(cvData)
                               })

cvData.list[['erp']] <- lapply(1:10,
                               function(k) {
                                   cvData <- loadData('erp', k)
                                   return(cvData)
                               })

fullData.list <- lapply(c('all', 'ern', 'erp', 'relapse'), loadData, k=-1)
names(fullData.list) <- c('all', 'ern', 'erp', 'relapse')

library(rCausalMGM)
library(survival)

temp <- fullData.list$ern$train[,c('DistantRelapse', 'DSS', 'LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'VCAM1', 'S100P', 'SYT17', 'FGFR4', 'POSTN', 'SERPINA1', 'SERPINA6')]

temp$DistantRelapse.time <- temp$DistantRelapse[,1]
temp$DistantRelapse.status <- temp$DistantRelapse[,2]

temp$DSS.time <- temp$DSS[,1]
temp$DSS.status <- temp$DSS[,2]

dr.res.cut <- surv_cutpoint(temp, 'DistantRelapse.time', 'DistantRelapse.status', variables=c('VCAM1', 'S100P', 'SYT17'), minprop=0.1)

plot(dr.res.cut, 'VCAM1')
plot(dr.res.cut, 'S100P')
plot(dr.res.cut, 'SYT17')

dr.res.cat <- surv_categorize(dr.res.cut)

dr.res.cat$LYMPH_NODE_STATUS <- relevel(temp$LYMPH_NODE_STATUS, 'NodeNegative')

lapply(dr.res.cat[,c('VCAM1', 'S100P', 'SYT17')], function(x) prop.table(table(x)))

png('VCAM1_DistantRelapse_metabric_ern.png', height=5, width=5, units='in', res=400)
ggsurvplot(survfit(Surv(DistantRelapse.time, DistantRelapse.status) ~ VCAM1, dr.res.cat), pval=T, conf.int=T)
dev.off()

ggsurvplot(survfit(Surv(DistantRelapse.time, DistantRelapse.status) ~ S100P, dr.res.cat), pval=T, conf.int=T, palette='Set1')
ggsurvplot(survfit(Surv(DistantRelapse.time, DistantRelapse.status) ~ SYT17, dr.res.cat), pval=T, conf.int=T, palette='Set1')

ggsurvplot(survfit(Surv(DistantRelapse.time, DistantRelapse.status) ~ VCAM1 + S100P + SYT17, dr.res.cat), pval=T)

dr.res.cat <- data.frame(dr.res.cat) %>%  mutate_if(is.character, factor, levels=c('low', 'high'))

coxph(Surv(DistantRelapse.time, DistantRelapse.status) ~ LYMPH_NODE_STATUS + VCAM1 + S100P + SYT17, dr.res.cat)


dss.res.cut <- surv_cutpoint(temp, 'DSS.time', 'DSS.status', variables=c('TUMOR_SIZE', 'FGFR4', 'POSTN', 'SERPINA1', 'SERPINA6', 'VCAM1'), minprop=0.1)

plot(dss.res.cut, 'TUMOR_SIZE')
plot(dss.res.cut, 'FGFR4')
plot(dss.res.cut, 'POSTN')
plot(dss.res.cut, 'SERPINA1')
plot(dss.res.cut, 'SERPINA6')
plot(dss.res.cut, 'VCAM1')

dss.res.cat <- surv_categorize(dss.res.cut)

lapply(dss.res.cat[,c('TUMOR_SIZE', 'FGFR4', 'POSTN', 'SERPINA1', 'SERPINA6')], function(x) prop.table(table(x)))

ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ TUMOR_SIZE, dss.res.cat), pval=T, conf.int=T)
ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ FGFR4, dss.res.cat), pval=T, conf.int=T)
ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ POSTN, dss.res.cat), pval=T, conf.int=T)
ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ SERPINA1, dss.res.cat), pval=T, conf.int=T)
ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ SERPINA6, dss.res.cat), pval=T, conf.int=T)
ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ VCAM1, dss.res.cat), pval=T, conf.int=T)

ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ TUMOR_SIZE + FGFR4 + POSTN + SERPINA1 + SERPINA6, dss.res.cat), pval=T)

dss.res.cat$LYMPH_NODE_STATUS <- relevel(temp$LYMPH_NODE_STATUS, 'NodeNegative')
dss.res.cat <- data.frame(dss.res.cat) %>%  mutate_if(is.character, factor, levels=c('low', 'high'))

coxph(Surv(DSS.time, DSS.status) ~ LYMPH_NODE_STATUS + TUMOR_SIZE + FGFR4 + POSTN + SERPINA1 + SERPINA6 + VCAM1, dss.res.cat)


numeric.mask <- sapply(fullData.list$erp$train.npn, function(x) !is.Surv(x) & is.numeric(x))

minCor <- 1
for (col in colnames(fullData.list$erp$train.npn)[numeric.mask]) {
    absCor <- abs(cor(fullData.list$erp$train.npn[,c('TUMOR_SIZE', col)]))[1,2]
    if (absCor < minCor & absCor > 0.05) {
        minCor <- absCor
        print(col)
    }
}

BayesIndTestMultiCoxTest('LYMPH_NODE_STATUS', 'GRADE', c('IGF1R'), fullData.list$erp$train.npn)

res <- lm(TUMOR_SIZE ~ MLPH + IGF1R, fullData.list$erp$train.npn)
res0 <- lm(TUMOR_SIZE ~ IGF1R, fullData.list$erp$train.npn)

summary(res)

deltaBIC <- BIC(res) - BIC(res0)

deltaBIC

1 / (1 + exp(-deltaBIC/2))

rf.out.list <- list()
for (idx in 1:10) {
    rf.out.list[[idx]] <- rfsrc(TUMOR_SIZE ~ .,
                                fullData.list$erp$train.npn %>% dplyr::select(-all_of(surv.feats)),
                                ntree=1000, importance=TRUE)
}


rf.out.null.list <- list()
for (idx in 1:10) {
    rf.out.null.list[[idx]] <- rfsrc(TUMOR_SIZE ~ .,
                                     fullData.list$erp$train.npn %>%
                                     dplyr::select(-all_of(c(surv.feats, 'LYMPH_NODE_STATUS'))),
                                     ntree=1000, importance=TRUE)
}


all.eval.df <- data.frame()
for (k in 1:10) {

    cvGraphs <- loadGraphs('all', k)

    fold.eval.df <- evaluateCausal(cvData.list[['all']][[k]]$train.npn,
                                   cvData.list[['all']][[k]]$test.npn,
                                   cvGraphs, 'all', c(3, 5, 8, 10, 15))

    fold.eval.df$Fold <- k

    all.eval.df <- rbind(all.eval.df, fold.eval.df)
}


all.summary <- all.eval.df %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`,
             FDR, Outcome, Features, Time, Metric) %>%
    summarize_at('Value', c(mean=mean, se=sd)) %>% as.data.frame

head(all.summary)

lapply(all.summary %>% filter(Metric=='Deviance', Time=='Full', grepl('All', Outcome)) %>% arrange(mean) %>% group_by(Outcome) %>% group_split(), function(x) head(as.data.frame(x)))

lapply(all.summary %>% filter((`MGM Select`=='BIC' & FDR %in% c(0.05, 0.2)) | (`MGM Select`=='StARS' & FDR==0.05), Time=='15 Years') %>% arrange(desc(mean)) %>% group_by(Outcome, Metric) %>% group_split(), function(x) as.data.frame(x))


all.summary %>%
    filter(Time=='Full', Features=='MB', Metric=="CPE") %>%
    group_by(Outcome) %>%
    arrange(desc(mean), .by_group=TRUE) %>%
    group_split() %>%
    lapply(as.data.frame)

all.total.dev <- all.eval.df %>%
    filter(Time=='Full', Features=='MB', Metric=='Deviance',
           Outcome%in%surv.feats[1:4]) %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`, FDR, Fold) %>%
    summarize_at('Value', c(sum=sum)) %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`, FDR) %>%
    summarize_at('sum', c(mean=mean, se=sd)) %>%
    arrange((mean)) %>%
    as.data.frame

head(all.total.dev)

#### All: StARS, FDR = 0.05

ggplot(all.total.dev,
       aes(x=1:nrow(all.total.dev), y=sum,
           ymin=sum-log(2^(10)), ymax=sum+log(2^(10)),
           color=`MGM Select`, shape=factor(FDR), lty=factor(FDR))) +
    geom_pointrange() +
    theme_bw()

ern.total.dev <- ern.eval.df %>%
    filter(Time=='Full', Features=='MB', Metric=='Deviance',
           Outcome%in%surv.feats[1:4]) %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`, FDR, Fold) %>%
    summarize_at('Value', c(sum=sum)) %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`, FDR) %>%
    summarize_at('sum', c(mean=mean, se=sd)) %>%
    arrange((mean)) %>%
    as.data.frame

ern.summary %>%
    filter(Time=='Full', Features=='MB', Metric=="Likelihood Ratio") %>%
    group_by(Outcome) %>%
    arrange(desc(mean), .by_group=TRUE) %>%
    group_split() %>%
    lapply(as.data.frame)

head(ern.total.dev)

#### ER-: StARS, FDR = 0.1

ggplot(ern.total.dev,
       aes(x=1:nrow(ern.total.dev), y=mean,
           ymin=mean-se, ymax=mean+se,
           color=`MGM Select`, shape=factor(FDR), lty=factor(FDR))) +
    geom_pointrange() +
    theme_bw()

ggplot(ern.total.dev,
       aes(x=1:nrow(ern.total.dev), y=sum,
           ymin=sum-log(2^(10)), ymax=sum+log(2^(10)),
           color=`MGM Select`, shape=factor(FDR), lty=factor(FDR))) +
    geom_pointrange() +
    theme_bw()

erp.total.dev <- erp.eval.df %>%
    filter(Time=='Full', Features=='NB', Metric=='Deviance',
           Outcome%in%surv.feats[1:4]) %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`, FDR, Fold) %>%
    summarize_at('Value', c(sum=sum)) %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`, FDR) %>%
    summarize_at('sum', c(mean=mean, se=sd)) %>%
    arrange((mean)) %>%
    as.data.frame

## erp.summary %>%
##    filter(Time=='Full', Features=='MB', Metric=="CPE") %>%
##    group_by(Outcome) %>%
##    arrange(desc(mean), .by_group=TRUE) %>%
##    group_split() %>%
##    lapply(as.data.frame)

head(erp.total.dev)

#### ER+: StARS, FDR = 0.1

ggplot(erp.total.dev,
       aes(x=1:nrow(erp.total.dev), y=mean,
           ymin=mean-se, ymax=mean+se,
           color=`MGM Select`, shape=factor(FDR), lty=factor(FDR))) +
    geom_pointrange() +
    theme_bw()

    

all.model.ranks <- list()
all.model.ranks.df <- data.frame()
for (sf in surv.feats[1:4]) {
    all.sum.split <- all.summary %>%
        filter(Time=='Full', Outcome==sf, Features=='MB', !Metric%in%c('Number of Features', 'log Partial Likelihood', "Uno's Concordance")) %>%
        group_by(Outcome, Time, Metric) %>%
        group_split()

    all.sum.split <- lapply(all.sum.split,
                            function(x) {
                                if (x$Metric[1] %in% c('Number of Features',
                                                       'Deviance')) {
                                    x$rank <- rank(x$mean)
                                } else {
                                    x$rank <- rank(-x$mean)
                                }
                                return(as.data.frame(x))
                            })

    compositeRank <- rep(0,nrow(all.sum.split[[1]]))
    for (idx in 1:length(all.sum.split)) {
        compositeRank <- compositeRank + all.sum.split[[idx]]$rank
    }

    final.all.sum.ranks <- data.frame(all.sum.split[[1]][,1:7],
                                      `Composite Rank`=compositeRank/length(all.sum.split),
                                      check.names=F)

    all.model.ranks[[sf]] <- final.all.sum.ranks
    if (sf=='DSS') {
        all.model.ranks.df <- final.all.sum.ranks[,-5]
    } else {
        all.model.ranks.df$`Composite Rank` <- all.model.ranks.df$`Composite Rank` + final.all.sum.ranks$`Composite Rank`
    }
}

lapply(all.model.ranks, function(x) head(x[order(x$`Composite Rank`),]))

all.model.ranks.df$`Composite Rank` <- all.model.ranks.df$`Composite Rank`/length(surv.feats)

all.model.ranks.df <- all.model.ranks.df[order(all.model.ranks.df$`Composite Rank`),]

head(all.model.ranks.df)

#### EBIC, FDR = 0.1

ern.eval.df <- data.frame()
for (k in 1:10) {

    cvGraphs <- loadGraphs('ern', k)

    fold.eval.df <- evaluateCausal(cvData.list[['ern']][[k]]$train.npn,
                                   cvData.list[['ern']][[k]]$test.npn,
                                   cvGraphs, 'ern', c(3, 5, 8, 10, 15))

    fold.eval.df$Fold <- k

    ern.eval.df <- rbind(ern.eval.df, fold.eval.df)
}


ern.summary <- ern.eval.df %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`,
             FDR, Outcome, Features, Time, Metric) %>%
    summarize_at('Value', c(mean=mean, se=sd)) %>% as.data.frame

head(ern.summary)

lapply(ern.summary %>% filter(Metric=='Deviance', Time=='Full', grepl('All', Outcome)) %>% arrange(mean) %>% group_by(Outcome) %>% group_split(), function(x) head(as.data.frame(x), 10))

lapply(ern.summary %>% filter(Metric=='Deviance', Time=='Full') %>% arrange(mean) %>% group_by(Outcome) %>% group_split(), function(x) head(as.data.frame(x), 10))

lapply(ern.summary %>% filter(`Split Status`=='split', (`MGM Select`=='EBIC' & FDR %in% c(0.2)) | (`MGM Select`=='StARS' & FDR==0.1) | (`MGM Select`=='BIC' & FDR==0.05), Time=='15 Years') %>% arrange(desc(mean)) %>% group_by(Outcome, Metric) %>% group_split(), function(x) as.data.frame(x))

lapply(ern.summary %>% filter(Metric=='% Deviance Explained') %>% arrange(desc(mean)) %>% group_by(Outcome) %>% group_split(), function(x) as.data.frame(x))


ern.model.ranks <- list()
ern.model.ranks.df <- data.frame()
for (sf in surv.feats[1:4]) {
    ern.sum.split <- ern.summary %>%
        filter(Time=='Full', Outcome==sf, Features=='MB', !Metric%in%c('Number of Features', 'log Partial Likelihood', "Uno's Concordance")) %>%
        group_by(Outcome, Time, Metric) %>%
        group_split()

    ern.sum.split <- lapply(ern.sum.split,
                            function(x) {
                                if (x$Metric[1] %in% c('Number of Features',
                                                       'Deviance')) {
                                    x$rank <- rank(x$mean)
                                } else {
                                    x$rank <- rank(-x$mean)
                                }
                                return(as.data.frame(x))
                            })

    compositeRank <- rep(0,nrow(ern.sum.split[[1]]))
    for (idx in 1:length(ern.sum.split)) {
        compositeRank <- compositeRank + ern.sum.split[[idx]]$rank
    }

    final.ern.sum.ranks <- data.frame(ern.sum.split[[1]][,1:7],
                                      `Composite Rank`=compositeRank/length(ern.sum.split),
                                      check.names=F)

    ern.model.ranks[[sf]] <- final.ern.sum.ranks
    if (sf=='DSS') {
        ern.model.ranks.df <- final.ern.sum.ranks[,-5]
    } else {
        ern.model.ranks.df$`Composite Rank` <- ern.model.ranks.df$`Composite Rank` + final.ern.sum.ranks$`Composite Rank`
    }
}

lapply(ern.model.ranks, function(x) head(x[order(x$`Composite Rank`),]))

ern.model.ranks.df$`Composite Rank` <- ern.model.ranks.df$`Composite Rank`/length(surv.feats)

ern.model.ranks.df <- ern.model.ranks.df[order(ern.model.ranks.df$`Composite Rank`),]

head(ern.model.ranks.df)

#### BIC, FDR = 0.2


erp.eval.df <- data.frame()
for (k in 1:10) {

    cvGraphs <- loadGraphs('erp', k)

    fold.eval.df <- evaluateCausal(cvData.list[['erp']][[k]]$train.npn,
                                   cvData.list[['erp']][[k]]$test.npn,
                                   cvGraphs, 'erp', c(3, 5, 8, 10, 15))
    fold.eval.df$Fold <- k

    erp.eval.df <- rbind(erp.eval.df, fold.eval.df)
}

erp.summary <- erp.eval.df %>%
    group_by(`Split Status`, `MGM Select`, `Orient Rule`,
             FDR, Outcome, Features, Time, Metric) %>%
    summarize_at('Value', c(mean=mean, se=sd)) %>% as.data.frame

head(erp.summary)

lapply(erp.summary %>% filter(Metric=='Deviance', Time=='Full', grepl('All', Outcome)) %>% arrange(mean) %>% group_by(Outcome) %>% group_split(), function(x) head(as.data.frame(x)))

lapply(erp.summary %>% filter(Metric=='Deviance', Time=='Full') %>% arrange(mean) %>% group_by(Outcome) %>% group_split(), function(x) head(as.data.frame(x)))


lapply(erp.summary %>% filter(Metric=="Uno's Concordance", Time=='15 Years') %>% arrange(mean) %>% group_by(Outcome) %>% group_split(), function(x) head(as.data.frame(x)))


lapply(erp.summary %>% filter((`MGM Select`=='StARS' & FDR==0.05) | (`MGM Select`=='StARS' & FDR==0.1), Time=='15 Years') %>% arrange(desc(mean)) %>% group_by(Outcome, Metric) %>% group_split(), function(x) as.data.frame(x))


erp.model.ranks <- list()
erp.model.ranks.df <- data.frame()
for (sf in surv.feats[1:4]) {
    erp.sum.split <- erp.summary %>%
        filter(Time=='Full', Outcome==sf, Features=='MB', !Metric%in%c('Number of Features', 'log Partial Likelihood', "Uno's Concordance")) %>%
        group_by(Outcome, Time, Metric) %>%
        group_split()

    erp.sum.split <- lapply(erp.sum.split,
                            function(x) {
                                if (x$Metric[1] %in% c('Number of Features',
                                                       'Deviance')) {
                                    x$rank <- rank(x$mean)
                                } else {
                                    x$rank <- rank(-x$mean)
                                }
                                return(as.data.frame(x))
                            })

    compositeRank <- rep(0,nrow(erp.sum.split[[1]]))
    for (idx in 1:length(erp.sum.split)) {
        compositeRank <- compositeRank + erp.sum.split[[idx]]$rank
    }

    final.erp.sum.ranks <- data.frame(erp.sum.split[[1]][,1:7],
                                      `Composite Rank`=compositeRank/length(erp.sum.split),
                                      check.names=F)

    erp.model.ranks[[sf]] <- final.erp.sum.ranks
    if (sf=='DSS') {
        erp.model.ranks.df <- final.erp.sum.ranks[,-5]
    } else {
        erp.model.ranks.df$`Composite Rank` <- erp.model.ranks.df$`Composite Rank` + final.erp.sum.ranks$`Composite Rank`
    }
}

lapply(erp.model.ranks, function(x) head(x[order(x$`Composite Rank`),]))

erp.model.ranks.df$`Composite Rank` <- erp.model.ranks.df$`Composite Rank`/length(surv.feats)

erp.model.ranks.df <- erp.model.ranks.df[order(erp.model.ranks.df$`Composite Rank`),]

head(erp.model.ranks.df)

#### StARS, FDR = 0.05

eval.df <- rbind(## cbind(`ER Status`='All', all.eval.df),
                 cbind(`ER Status`='ER-', ern.eval.df),
                 cbind(`ER Status`='ER+', erp.eval.df))

write.csv(eval.df, '10fold_Causal_results_v3.csv')

eval.df <- read.csv('10fold_Causal_results_v3.csv', row.names=1, check.names=F)

head(eval.df)

allModel.eval.df <- rbind(cbind(Model='CausalCoxMGM',
                                eval.df %>%
                                filter(`MGM Select`=='StARS',
                                       FDR==0.1, `Orient Rule`=='maxp',
                                       `Split Status`=='split',
                                       `ER Status` != "All",
                                       Outcome %in% surv.feats[1:4]) %>%
                          select(-c('MGM Select', 'FDR', 'Orient Rule', 'Split Status'))),
                          cbind(Model='LASSO Cox Min',
                                eval.lasso.df %>% filter(Features=='Min',
                                                         `ER Status` != "All",
                                                         Outcome %in% surv.feats[1:4])),
                          cbind(Model='LASSO Cox 1SE',
                                eval.lasso.df %>% filter(Features=='1SE',
                                                         `ER Status` != "All",
                                                         Outcome %in% surv.feats[1:4])))

allModel.eval.df$Time <- factor(allModel.eval.df$Time,
                                c('Full', paste(c(3,5,8,10,15), 'Years')))

allModel.eval.df$`ER Status` <- factor(allModel.eval.df$`ER Status`,
                                       c('ER+', 'ER-'))

allModel.eval.df$`Outcome` <- factor(allModel.eval.df$`Outcome`,
                                     surv.feats[1:4],
                                     labels=c("DSS", "OD", "Distant", "Local"))

allModel.eval.df$`Model` <- factor(allModel.eval.df$`Model`,
                                   c('CausalCoxMGM', 'LASSO Cox Min', 'LASSO Cox 1SE'),
                                   labels=c('Causal', 'LASSO Min', 'LASSO 1SE'))

harrell.cv.wald.df <- allModel.eval.df %>%
  filter(Metric=="Harrell's Concordance",
         Time %in% c(paste(c(5,10,15), 'Years')),
         !is.na(Outcome)) %>%
    group_by(Fold, `ER Status`, Outcome, Time, Metric) %>%
    mutate(Delta = Value - Value[1]) %>%
    filter(Model != 'Causal') %>%
    group_by(Model, `ER Status`, Outcome, Time, Metric) %>%
    summarize(mean=mean(Delta, na.rm=T),
              se=sqrt(1/10 + 1/9) * ifelse(sd(Delta, na.rm=T)==0, 1, sd(Delta, na.rm=T))) %>%
    mutate(tstat=mean/se,
           pval=2*pt(-abs(tstat), 9)) %>%
    as.data.frame

harrell.cv.wald.df <- harrell.cv.wald.df %>%
    mutate(padj=p.adjust(pval, method='fdr', n=n()),
           signif=ifelse(padj<0.05,
                  ifelse(padj<0.01,
                  ifelse(padj<0.001, '***', '**'), '*'), '')) %>%
    as.data.frame

harrell.cv.wald.df
min(harrell.cv.wald.df$padj)

allModel.summary <- allModel.eval.df %>%
    filter(Metric=="Harrell's Concordance",
           Time %in% c(paste(c(5,10,15), 'Years')),
           !is.na(Outcome)) %>%
    group_by(Model, `ER Status`, Outcome, Time, Metric) %>% 
    summarize(mean=mean(Value, na.rm=T),
              se=sqrt(1/10 + 1/9) * sd(Value, na.rm=T)) %>%
    as.data.frame

## allModel.summary$padj <- NA
## allModel.summary$signif <- ''

all(allModel.summary[25:nrow(allModel.summary),1:5] == harrell.cv.wald.df[,1:5])

allModel.summary[25:nrow(allModel.summary),c('pval')] <- harrell.cv.wald.df[,c('pval')]


time <- '10 Years'

ggplot(allModel.summary %>% filter(Time==time),
       aes(x=Model, y=mean, ymin=mean-1.96*se, ymax=mean+1.96*se,
           color=Model)) +
geom_pointrange(fatten=2) +
facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
scale_color_brewer(palette='Set1') +
ylab("Concordance") +
geom_hline(yintercept=0.5, lty=5, color='darkgrey') + 
theme_bw(base_size=7) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
geom_bracket(aes(xmax=Model, xmin='Causal', y.position=0.8, color='black', label=paste('p =', round(padj,4)), color='black'),
             allModel.summary %>% filter(Time==time) %>% mutate(padj=p.adjust(pval, method='fdr')) %>% filter(padj < 0.05), step.increase=0.075)


ggplot(allModel.summary,
       aes(x=Time, y=mean, ymin=mean-1.96*se, ymax=mean+1.96*se,
           color=Model)) +
geom_pointrange(fatten=1, position=position_dodge(width=0.9)) +
facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
scale_color_brewer(palette='Set1') +
ylab("Concordance") +
geom_hline(yintercept=0.5, lty=5, color='darkgrey') + 
theme_bw(base_size=7) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('cv_harrell_conc_results_signif.png', width=5, height=2.25)

ggplot(allModel.eval.df %>%
       filter(Metric=="Harrell's Concordance",
              Time %in% paste(c(5,10,15), 'Years'),
              !is.na(Outcome)),
              ## Time=='Full',),
       aes(x=Time, y=Value, color=Model)) +
    geom_pointrange(stat='summary', fun.data='mean_sdl',
                    fun.args=list(mult=1.96, na.rm=T),
                    position=position_dodge(width=0.9),
                    fatten=1) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    scale_color_brewer(palette='Set1') +
    ylab("Concordance") +
    geom_hline(yintercept=0.5, lty=5, color='darkgrey') + 
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
geom_bracket(aes(xmin='Causal', xmax=Model, y=1, label=signif), data=harrell.cv.wald.df)

ggsave('cv_harrell_conc_results_signif.png', width=5, height=3)


ggplot(allModel.eval.df %>%
       filter(Metric=="Uno's Concordance",
              Time %in% paste(c(5,10,15), 'Years'),
              !is.na(Outcome)),
       aes(x=Time, y=Value, color=Model)) +
    geom_pointrange(stat='summary', fun.data='mean_sdl', fatten=2,
                    fun.args=list(mult=1.96, na.rm=T), position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    scale_color_brewer(palette='Set1') +
    ylab("Uno's Concordance") +
    geom_hline(yintercept=0.5, lty=5, color='darkgrey') +
    theme_bw(base_size=7) + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

ggsave('cv_uno_conc_results.png', width=10, height=4.5)



numfeats.cv.wald.df <- allModel.eval.df %>%
  filter(Metric=="Number of Features",
         !is.na(Outcome)) %>%
    group_by(Fold, `ER Status`, Outcome, Time, Metric) %>%
    mutate(Delta = Value - Value[1]) %>%
    filter(Model != 'Causal') %>%
    group_by(Model, `ER Status`, Outcome, Time, Metric) %>%
    summarize(mean=mean(Delta, na.rm=T),
              se=sqrt(1/10 + 1/9) * ifelse(sd(Delta, na.rm=T)==0, 1, sd(Delta, na.rm=T))) %>%
    mutate(tstat=mean/se,
           pval=2*pt(-abs(tstat), 9)) %>%
    as.data.frame

numfeats.cv.wald.df <- numfeats.cv.wald.df %>%
    mutate(padj=p.adjust(pval, method='fdr', n=n()),
           signif=ifelse(padj<0.05,
                  ifelse(padj<0.01,
                  ifelse(padj<0.001, '***', '**'), '*'), '')) %>%
    as.data.frame

numfeats.cv.wald.df
min(numfeats.cv.wald.df$padj)

allModel.summary <- allModel.eval.df %>%
    filter(Metric=="Number of Features",
           !is.na(Outcome)) %>%
    group_by(Model, `ER Status`, Outcome, Time, Metric) %>% 
    summarize(mean=mean(Value, na.rm=T),
              se=sqrt(1/10 + 1/9) * sd(Value, na.rm=T)) %>%
    as.data.frame

## allModel.summary$padj <- NA
## allModel.summary$signif <- ''

all(allModel.summary[9:nrow(allModel.summary),1:5] == numfeats.cv.wald.df[,1:5])

allModel.summary[9:nrow(allModel.summary),c('pval', 'padj', 'signif')] <- numfeats.cv.wald.df[,c('pval', 'padj', 'signif')]


ggplot(allModel.summary,
       aes(x=Model, y=mean, ymin=mean-1.96*se, ymax=mean+1.96*se,
           fill=Model)) +
geom_bar(stat='identity') +
geom_errorbar(stat='identity', width=0.5) +
facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
scale_fill_brewer(palette='Set1') +
ylab("Number of Features") +
ylim(c(0, 70)) +
geom_bracket(aes(xmax=Model, xmin='Causal', y.position=ceiling(mean + 2 * se) + 2, color='black', label=paste('p =', round(padj,4)), color='black'),
             allModel.summary %>% filter(padj < 0.05) %>% group_by(Outcome, `ER Status`) %>% mutate_at(c('mean', 'se'), function(x) max(x)), step.increase=0.05, label.size=1.8) +
theme_bw(base_size=7) +
theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

ggsave('cv_numfeats_results_signif.png', width=5, height=1.75)

ggplot(allModel.eval.df %>%
       filter(Metric=="Number of Features",
              Time=='Full',
              !is.na(Outcome)),
       aes(x=Model, y=Value, fill=Model)) +
    geom_bar(stat='summary') + 
    geom_errorbar(stat='summary', fun.data='mean_sdl',
                    fun.args=list(mult=1.96, na.rm=T), width=0.5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    scale_fill_brewer(palette='Set1') +
    ylab("Number of Features") +
    theme_bw(base_size=7) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

ggsave('cv_numfeats_results_signif.png', width=7, height=3.5)

ggplot(ern.eval.df %>% filter(Time=='5 Years', `Split Status`=='split'),
       aes(x=factor(FDR), y=Value, fill=`MGM Select`)) +
    geom_boxplot() +
    facet_grid(cols=vars(Outcome), rows=vars(Metric), scales='free') +
    theme_bw()

ern.eval.df %>% filter(Time=='3 Years', `Split Status`=='split', Metric=="Uno's Concordance", Value<0.3)

ggplot(erp.eval.df %>% filter(Time!='Full', Features=='MB',
                              !Metric%in%c('log Partial Likelihood', 'Number of Features')),
       aes(x=factor(Time, levels=paste(times, 'Years')), y=Value,
           color=`MGM Select`,
           shape=factor(FDR), lty=factor(FDR))) +
    geom_pointrange(stat='summary',
                    position=position_dodge(width=0.8),
                    fun.data='mean_sdl',
                    fun.args=list(mult=1)) +
    facet_grid(cols=vars(Outcome), rows=vars(Metric), scales='free') +
    theme_bw()


#### LASSO Evaluation

all.eval.lasso.df <- data.frame()
for (k in 1:10) {
    cvData <- loadData('all', k)
    ## train.npn <- cvData$train %>%
    ##     mutate_if(function(x) !is.Surv(x) & is.numeric(x),
    ##               function(x) huge::huge.npn(as.matrix(x)))
    ## test.npn <- cvData$test %>%
    ##     mutate_if(function(x) !is.Surv(x) & is.numeric(x),
    ##               function(x) huge::huge.npn(as.matrix(x)))

    fold.eval.df <- evaluateLASSO(train.npn, test.npn, 'all', c(5, 10, 15))

    fold.eval.df$Fold <- k

    all.eval.lasso.df <- rbind(all.eval.lasso.df, fold.eval.df)
}


all.lasso.summary <- all.eval.lasso.df %>%
    group_by(Outcome, Features, Time, Metric) %>%
    summarize_at('Value', c(mean=mean, se=sd)) %>% as.data.frame

head(all.lasso.summary)


ern.eval.lasso.df.v2 <- data.frame()
for (k in 1:10) {
    cvData <- loadData('ern', k)
    ## train.npn <- cvData$train %>%
    ##     mutate_if(function(x) !is.Surv(x) & is.numeric(x),
    ##               function(x) huge::huge.npn(as.matrix(x)))
    ## test.npn <- cvData$test %>%
    ##     mutate_if(function(x) !is.Surv(x) & is.numeric(x),
    ##               function(x) huge::huge.npn(as.matrix(x)))

    fold.eval.df <- evaluateLASSO(cvData$train.npn, cvData$test.npn,
                                  'ern', c(3, 5, 8, 10, 15))

    fold.eval.df$Fold <- k

    ern.eval.lasso.df.v2 <- rbind(ern.eval.lasso.df.v2, fold.eval.df)
}


ern.lasso.summary <- ern.eval.lasso.df %>%
    group_by(Outcome, Features, Time, Metric) %>%
    summarize_at('Value', c(mean=mean, se=sd)) %>% as.data.frame

(ern.lasso.summary %>% filter(Outcome=='DSS', Metric=="Uno's Concordance"))

ern.lasso.summary.v2 <- ern.eval.lasso.df.v2 %>%
    group_by(Outcome, Features, Time, Metric) %>%
    summarize_at('Value', c(mean=mean, se=sd)) %>% as.data.frame

(ern.lasso.summary.v2 %>% filter(Outcome=='DSS', Metric=="Uno's Concordance"))


lapply(ern.lasso.summary %>% filter(Time=='Full', Outcome %in% surv.feats[1:4]) %>% arrange(desc(mean)) %>% group_by(Outcome, Metric) %>% group_split(), function(x) head(as.data.frame(x)))


lapply(ern.lasso.summary %>% filter(Time=='15 Years', Outcome %in% surv.feats[1:4]) %>% arrange(desc(mean)) %>% group_by(Outcome, Metric) %>% group_split(), function(x) head(as.data.frame(x)))

lapply(ern.lasso.summary.v2 %>% filter(Time=='15 Years', Outcome %in% surv.feats[1:4]) %>% arrange(desc(mean)) %>% group_by(Outcome, Metric) %>% group_split(), function(x) head(as.data.frame(x)))




erp.eval.lasso.df <- data.frame()
for (k in 1:10) {
    cvData <- loadData('erp', k)
    ## train.npn <- cvData$train %>%
    ##     mutate_if(function(x) !is.Surv(x) & is.numeric(x),
    ##               function(x) huge::huge.npn(as.matrix(x)))
    ## test.npn <- cvData$test %>%
    ##     mutate_if(function(x) !is.Surv(x) & is.numeric(x),
    ##               function(x) huge::huge.npn(as.matrix(x)))

    fold.eval.df <- evaluateLASSO(cvData$train.npn, cvData$test.npn,
                                  'erp', c(5, 10, 15))

    fold.eval.df$Fold <- k

    erp.eval.lasso.df <- rbind(erp.eval.lasso.df, fold.eval.df)
}


erp.lasso.summary <- erp.eval.lasso.df %>%
    group_by(Outcome, Features, Time, Metric) %>%
    summarize_at('Value', c(mean=mean, se=sd)) %>% as.data.frame

head(erp.lasso.summary)

lapply(erp.lasso.summary %>% filter(Time=='Full', Outcome %in% surv.feats[1:4]) %>% arrange(desc(mean)) %>% group_by(Outcome, Metric) %>% group_split(), function(x) head(as.data.frame(x)))


eval.lasso.df <- rbind(## cbind(`ER Status`='All', all.eval.lasso.df),
                       cbind(`ER Status`='ER-', ern.eval.lasso.df),
                       cbind(`ER Status`='ER+', erp.eval.lasso.df))

write.csv(eval.lasso.df, '10fold_LASSO_results.csv')

eval.lasso.df <- read.csv('10fold_LASSO_results.csv', row.names=1, check.names=FALSE)

ggplot(eval.lasso.df %>% filter(Time!='Full'),
       aes(x=factor(Time, levels=c(paste(times, 'Years'))),
           y=Value, color=`ER Status`, shape=Features, lty=Features)) +
    geom_pointrange(stat='summary', fun.data='mean_sdl',
                    fun.args=list(mult=1), position=position_dodge(width=0.8)) +
    facet_grid(rows=vars(Metric), cols=vars(Outcome)) +
    theme_bw()

all.conc.df <- data.frame()
for (k in 1:10) {
    cvData <- loadData('all', k)
    train.npn <- cvData$train %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))
    test.npn <- cvData$test %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))

    cvGraphs <- loadGraphs('all', k)

    for (sf in c('DSS','OD','DistantRelapse','LocalRelapse','OS','DRFS','DFS')) {
        for (sel.meth in c('StEPS', 'StARS')) {
            if (sf %in% surv.feats[1:4]) {
                nb <- sort(getNeighbors(cvData$graph[[sel.meth]], sf))
                mb <- sort(setdiff(cvData$graph[[sel.meth]]$markov.blankets[[sf]],
                                   surv.feats))
                mb.bic <- sort(growShrinkMB(cvData$train, sf,
                                            graph=cvData$graph[[sel.meth]],
                                            rank=T, verbose=T))
                ## mb.aic <- sort(growShrinkMB(cvData$train, sf,
                ##                             graph=cvData$graph[[sel.meth]],
                ##                             penalty=2/log(nrow(cvData$train)),
                ##                             rank=T, verbose=T))
            } else {
                ## nb <- paste0(surv.feats[1:(which(sf==surv.feats)-3)], '.risk.nb')
                ## mb <- paste0(surv.feats[1:(which(sf==surv.feats)-3)], '.risk.mb')
                
                nb <- sort(unique(unlist(lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                                                function(x) getNeighbors(cvData$graph[[sel.meth]], x)))))
                mb <- setdiff(
                    sort(
                        unique(
                            unlist(
                                cvData$graph[[sel.meth]]$markov.blankets[surv.feats[1:(which(sf==surv.feats)-3)]]
                            )
                        )
                    ),
                    surv.feats)

                mb.bic <- setdiff(
                    sort(
                        unique(
                            unlist(
                                lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                                       function(x) growShrinkMB(cvData$train, x, graph=cvData$graph[[sel.meth]], rank=T, verbose=T))
                            )
                        )
                    ),
                    surv.feats)

                ## mb.aic <- setdiff(
                ##     sort(
                ##         unique(
                ##             unlist(
                ##                 lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                ##                        function(x) growShrinkMB(cvData$train, x, graph=cvData$graph[[sel.meth]], rank=T, penalty=2/log(nrow(cvData$train)), verbose=T))
                ##             )
                ##         )
                ##     ),
                ##     surv.feats)
            }
            
            f.nb <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(c(1, nb), collapse=' + ')))
            f.mb <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(c(1, mb), collapse=' + ')))
            f.mb.bic <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(c(1, mb.bic), collapse=' + ')))

            ## f.mb.aic <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(c(1, mb.aic), collapse=' + ')))

            cph.nb <- coxph(f.nb, train.npn)
            cph.mb <- coxph(f.mb, train.npn)
            cph.mb.bic <- coxph(f.mb.bic, train.npn)
            ## cph.mb.aic <- coxph(f.mb.aic, train.npn)

            f.nb <- as.formula(paste0('~ ', paste(c(1, nb), collapse=' + ')))
            f.mb <- as.formula(paste0('~ ', paste(c(1, mb), collapse=' + ')))
            f.mb.bic <- as.formula(paste0('~ ', paste(c(1, mb.bic), collapse=' + ')))
            ## f.mb.aic <- as.formula(paste0('~ ', paste(c(1, mb.aic), collapse=' + ')))

            if (length(nb)!=0) {
                ridge.nb <- glmnet::glmnet(as.matrix(model.matrix(f.nb, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            if (length(mb)!=0) {
                ridge.mb <- glmnet::glmnet(as.matrix(model.matrix(f.mb, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            if (length(mb.bic)!=0) {
                ridge.mb.bic <- glmnet::glmnet(as.matrix(model.matrix(f.mb.bic, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            ## if (length(mb.aic)!=0) {
            ##     ridge.mb.aic <- glmnet::glmnet(as.matrix(model.matrix(f.mb.aic, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            ## }
            
            cv.lasso <- glmnet::cv.glmnet(model.matrix(~., train.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats, paste0(surv.feats,'.risk.nb'), paste0(surv.feats,'.risk.mb'))])[,-1], train.npn[,sf], family='cox')

            plot(cv.lasso)

            eval.df <- evaluateModel(cph.nb, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, 'NB')
            eval.df$`Number of Features` <- length(nb)
            eval.df$Fold <- k
            all.conc.df <- rbind(all.conc.df, eval.df)

            eval.df <- evaluateModel(cph.mb, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, "MB")
            eval.df$`Number of Features` <- length(mb)
            eval.df$Fold <- k
            all.conc.df <- rbind(all.conc.df, eval.df)

            eval.df <- evaluateModel(cph.mb.bic, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, "BIC MB")
            eval.df$`Number of Features` <- length(mb.bic)
            eval.df$Fold <- k
            all.conc.df <- rbind(all.conc.df, eval.df)

            ## eval.df <- evaluateModel(cph.mb.aic, test.npn, c(2, 5, 10, 15), sf)
            ## eval.df$Model <- paste(sel.meth, "AIC MB")
            ## eval.df$`Number of Features` <- length(mb.aic)
            ## eval.df$Fold <- k
            ## all.conc.df <- rbind(all.conc.df, eval.df)

            if (length(nb)!=0) {
                eval.df <- evaluateModel(ridge.nb, test.npn, c(2, 5, 10, 15), sf, f.nb)
                eval.df$Model <- paste("Ridge", sel.meth, "NB")
                eval.df$`Number of Features` <- length(nb)
                eval.df$Fold <- k
                all.conc.df <- rbind(all.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.nb, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "NB")
                eval.df$`Number of Features` <- length(nb)
                eval.df$Fold <- k
                all.conc.df <- rbind(all.conc.df, eval.df)
            }

            if (length(mb)!=0) {
                eval.df <- evaluateModel(ridge.mb, test.npn, c(2, 5, 10, 15), sf, f.mb)
                eval.df$Model <- paste("Ridge", sel.meth, "MB")
                eval.df$`Number of Features` <- length(mb)
                eval.df$Fold <- k
                all.conc.df <- rbind(all.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.mb, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "MB")
                eval.df$`Number of Features` <- length(mb)
                eval.df$Fold <- k
                all.conc.df <- rbind(all.conc.df, eval.df)
            }

            if (length(mb.bic)!=0) {
                eval.df <- evaluateModel(ridge.mb.bic, test.npn, c(2, 5, 10, 15), sf, f.mb.bic)
                eval.df$Model <- paste("Ridge", sel.meth, "BIC MB")
                eval.df$`Number of Features` <- length(mb.bic)
                eval.df$Fold <- k
                all.conc.df <- rbind(all.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.mb.bic, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "BIC MB")
                eval.df$`Number of Features` <- length(mb.bic)
                eval.df$Fold <- k
                all.conc.df <- rbind(all.conc.df, eval.df)
            }

            ## if (length(mb.aic)!=0) {
            ##     eval.df <- evaluateModel(ridge.mb.aic, test.npn, c(2, 5, 10, 15), sf, f.mb.aic)
            ##     eval.df$Model <- paste("Ridge", sel.meth, "AIC MB")
            ##     eval.df$`Number of Features` <- length(mb.aic)
            ##     eval.df$Fold <- k
            ##     all.conc.df <- rbind(all.conc.df, eval.df)
            ## } else {
            ##     eval.df <- evaluateModel(cph.mb.aic, test.npn, c(2, 5, 10, 15), sf)
            ##     eval.df$Model <- paste("Ridge", sel.meth, "AIC MB")
            ##     eval.df$`Number of Features` <- length(mb.aic)
            ##     eval.df$Fold <- k
            ##     all.conc.df <- rbind(all.conc.df, eval.df)
            ## }
        }

        eval.df <- evaluateModel(cv.lasso, test.npn, c(2, 5, 10, 15), sf)
        eval.df$Fold <- k
        all.conc.df <- rbind(all.conc.df, eval.df)

        
    }
}

## colMeans(cv.conc.dss)
## apply(cv.conc.dss, 2, sd)
all.df <- data.frame()
for (sf in surv.feats) {
    for (t in c('2', '5', '10', '15')) {
        all.df <- rbind(all.df,
                        data.frame(Concordance=colMeans(cv.conc.all[[sf]][[t]]),
                                   SD=apply(cv.conc.all[[sf]][[t]], 2, sd),
                                   Method=factor(c('Causal Neighbors', 'Causal MB',
                                                   'Causal Neighbors Risk',
                                                   'Causal MB Risk',
                                                   'Ridge Causal Neighbors',
                                                   'Ridge Causal MB',
                                                   'LASSO 1SE', 'LASSO Min'),
                                                 levels=c('Causal Neighbors', 'Causal MB',
                                                          'Causal Neighbors Risk',
                                                          'Causal MB Risk',
                                                          'Ridge Causal Neighbors',
                                                          'Ridge Causal MB',
                                                          'LASSO 1SE', 'LASSO Min')),
                                   Time=factor(rep(paste(t, 'Years'), 8),
                                               levels=paste(c(2,5,10,15), 'Years')),
                                   Outcome=factor(rep(sf, 8),
                                                  levels=c('DSS','OD','DistantRelapse','LocalRelapse','OS','DRFS','DFS')),
                                   ER_STATUS=rep('All', 8)))
    }
}

all.conc.df$Target <- factor(all.conc.df$Target, surv.feats)
all.conc.df$Time <- factor(all.conc.df$Time, paste(c(2,5,10,15), 'Years'))

write.csv(all.conc.df, 'cv_all_concordance.csv')

gg <- ggplot(all.conc.df %>% filter(Model %in% c('LASSO 1SE', 'StEPS BIC MB'), Target!='LocalRelapse') %>% group_by(Model, Target, Time), aes(Time, `Uno Concordance`, color=Model)) +
    geom_point(stat='summary', fun.data='mean_sdl', position=position_dodge2(width=0.5)) +
    geom_errorbar(stat='summary', width=0.5, position=position_dodge2(width=0.5), fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    geom_hline(yintercept=0.5, color='darkgray', lty=5) + 
    facet_grid(cols=vars(Target), scales='free_y') +
    ## stat_compare_means(paired=T) + 
    scale_color_brewer(palette='Set1', labels = c("LASSO Cox", "CausalCoxMGM")) +
    ylab("Uno's Concordance") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

ggsave('uno_concordance_cv_all_plot.pdf', gg, height=5, width=10)
ggsave('uno_concordance_cv_all_plot.png', gg, height=5, width=10, dpi=400)

gg <- ggplot(all.conc.df %>% filter(Time=='2 Years', Model %in% c('LASSO 1SE', 'StEPS BIC MB'), Target!='LocalRelapse') %>% group_by(Model, Target, Time), aes(Model, `Number of Features`, fill=Model)) +
    geom_bar(stat='summary', fun.data='mean_sdl') +
    geom_errorbar(stat='summary', width=0.2, fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    facet_grid(cols=vars(Target)) +
    scale_fill_brewer(palette='Set1', labels = c("LASSO Cox", "CausalCoxMGM")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

ggsave('numFeats_cv_all_plot.pdf', gg, height=5, width=10)
ggsave('numFeats_cv_all_plot.png', gg, height=5, width=10, dpi=400)


colMeans(cv.conc.dss.npn)
apply(cv.conc.dss.npn, 2, sd)

dss.df <- data.frame(Mean=colMeans(cv.conc.dss.npn),
                             SD=apply(cv.conc.dss.npn, 2, sd),
                             Method=factor(c('Causal Neigbors', 'Causal MB',
                                             'Ridge Causal Neighbors', 'Ridge Causal MB',
                                             'LASSO 1SE', 'LASSO Min'),
                                           levels=c('Causal Neigbors', 'Causal MB',
                                                    'Ridge Causal Neighbors',
                                                    'Ridge Causal MB',
                                                    'LASSO 1SE', 'LASSO Min')),
                             Outcome=rep('DSS Mortailty', 6),
                             ER_STATUS=rep('All', 6))

ggplot(dss.df, aes(Method, Mean, color=Method)) +
    geom_point() +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD)) +
    ylab('Concordance') +
    facet_grid(vars(ER_STATUS), vars(Outcome)) +
    scale_color_brewer(palette='Paired') +
    theme_bw()



ern.conc.df <- data.frame()
for (k in 1:5) {
    cvData <- loadGraphAndData('ern', k)
    train.npn <- cvData$train %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))
    test.npn <- cvData$test %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))
    
    test.npn$ER_STATUS <- 'Negative'

    for (sf in c('DSS','OD','DistantRelapse','LocalRelapse','OS','DRFS','DFS')) {
        for (sel.meth in c('StEPS', 'StARS')) {
            if (sf %in% surv.feats[1:4]) {
                nb <- sort(getNeighbors(cvData$graph[[sel.meth]], sf))
                mb <- sort(setdiff(cvData$graph[[sel.meth]]$markov.blankets[[sf]],
                                   surv.feats))
                mb.bic <- sort(growShrinkMB(cvData$train, sf,
                                            graph=cvData$graph[[sel.meth]],
                                            rank=T, verbose=T))
                ## mb.aic <- sort(growShrinkMB(cvData$train, sf,
                ##                             graph=cvData$graph[[sel.meth]],
                ##                             penalty=2/log(nrow(cvData$train)),
                ##                             rank=T, verbose=T))
            } else {
                ## nb <- paste0(surv.feats[1:(which(sf==surv.feats)-3)], '.risk.nb')
                ## mb <- paste0(surv.feats[1:(which(sf==surv.feats)-3)], '.risk.mb')
                
                nb <- sort(unique(unlist(lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                                                function(x) getNeighbors(cvData$graph[[sel.meth]], x)))))
                mb <- setdiff(
                    sort(
                        unique(
                            unlist(
                                cvData$graph[[sel.meth]]$markov.blankets[surv.feats[1:(which(sf==surv.feats)-3)]]
                            )
                        )
                    ),
                    surv.feats)

                mb.bic <- setdiff(
                    sort(
                        unique(
                            unlist(
                                lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                                       function(x) growShrinkMB(cvData$train, x, graph=cvData$graph[[sel.meth]], rank=T, verbose=T))
                            )
                        )
                    ),
                    surv.feats)

                ## mb.aic <- setdiff(
                ##     sort(
                ##         unique(
                ##             unlist(
                ##                 lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                ##                        function(x) growShrinkMB(cvData$train, x, graph=cvData$graph[[sel.meth]], rank=T, penalty=2/log(nrow(cvData$train)), verbose=T))
                ##             )
                ##         )
                ##     ),
                ##     surv.feats)
            }
            
            f.nb <- as.formula(paste0(sf, ' ~ ', paste(c(1, nb), collapse=' + ')))
            f.mb <- as.formula(paste0(sf, ' ~ ', paste(c(1, mb), collapse=' + ')))
            f.mb.bic <- as.formula(paste0(sf, ' ~ ', paste(c(1, mb.bic), collapse=' + ')))

            ## f.mb.aic <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(c(1, mb.aic), collapse=' + ')))

            cph.nb <- coxph(f.nb, train.npn)
            cph.mb <- coxph(f.mb, train.npn)
            cph.mb.bic <- coxph(f.mb.bic, train.npn)
            ## cph.mb.aic <- coxph(f.mb.aic, train.npn)

            f.nb <- as.formula(paste0('~ ', paste(c(1, nb), collapse=' + ')))
            f.mb <- as.formula(paste0('~ ', paste(c(1, mb), collapse=' + ')))
            f.mb.bic <- as.formula(paste0('~ ', paste(c(1, mb.bic), collapse=' + ')))
            ## f.mb.aic <- as.formula(paste0('~ ', paste(c(1, mb.aic), collapse=' + ')))

            if (length(nb)!=0) {
                ridge.nb <- glmnet::glmnet(as.matrix(model.matrix(f.nb, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            if (length(mb)!=0) {
                ridge.mb <- glmnet::glmnet(as.matrix(model.matrix(f.mb, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            if (length(mb.bic)!=0) {
                ridge.mb.bic <- glmnet::glmnet(as.matrix(model.matrix(f.mb.bic, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            ## if (length(mb.aic)!=0) {
            ##     ridge.mb.aic <- glmnet::glmnet(as.matrix(model.matrix(f.mb.aic, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            ## }
            
            eval.df <- evaluateModel(cph.nb, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, 'NB')
            eval.df$`Number of Features` <- length(nb)
            eval.df$Fold <- k
            ern.conc.df <- rbind(ern.conc.df, eval.df)

            eval.df <- evaluateModel(cph.mb, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, "MB")
            eval.df$`Number of Features` <- length(mb)
            eval.df$Fold <- k
            ern.conc.df <- rbind(ern.conc.df, eval.df)

            eval.df <- evaluateModel(cph.mb.bic, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, "BIC MB")
            eval.df$`Number of Features` <- length(mb.bic)
            eval.df$Fold <- k
            ern.conc.df <- rbind(ern.conc.df, eval.df)

            ## eval.df <- evaluateModel(cph.mb.aic, test.npn, c(2, 5, 10, 15), sf)
            ## eval.df$Model <- paste(sel.meth, "AIC MB")
            ## eval.df$`Number of Features` <- length(mb.aic)
            ## eval.df$Fold <- k
            ## ern.conc.df <- rbind(ern.conc.df, eval.df)

            if (length(nb)!=0) {
                eval.df <- evaluateModel(ridge.nb, test.npn, c(2, 5, 10, 15), sf, f.nb)
                eval.df$Model <- paste("Ridge", sel.meth, "NB")
                eval.df$`Number of Features` <- length(nb)
                eval.df$Fold <- k
                ern.conc.df <- rbind(ern.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.nb, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "NB")
                eval.df$`Number of Features` <- length(nb)
                eval.df$Fold <- k
                ern.conc.df <- rbind(ern.conc.df, eval.df)
            }

            if (length(mb)!=0) {
                eval.df <- evaluateModel(ridge.mb, test.npn, c(2, 5, 10, 15), sf, f.mb)
                eval.df$Model <- paste("Ridge", sel.meth, "MB")
                eval.df$`Number of Features` <- length(mb)
                eval.df$Fold <- k
                ern.conc.df <- rbind(ern.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.mb, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "MB")
                eval.df$`Number of Features` <- length(mb)
                eval.df$Fold <- k
                ern.conc.df <- rbind(ern.conc.df, eval.df)
            }

            if (length(mb.bic)!=0) {
                eval.df <- evaluateModel(ridge.mb.bic, test.npn, c(2, 5, 10, 15), sf, f.mb.bic)
                eval.df$Model <- paste("Ridge", sel.meth, "BIC MB")
                eval.df$`Number of Features` <- length(mb.bic)
                eval.df$Fold <- k
                ern.conc.df <- rbind(ern.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.mb.bic, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "BIC MB")
                eval.df$`Number of Features` <- length(mb.bic)
                eval.df$Fold <- k
                ern.conc.df <- rbind(ern.conc.df, eval.df)
            }

            ## if (length(mb.aic)!=0) {
            ##     eval.df <- evaluateModel(ridge.mb.aic, test.npn, c(2, 5, 10, 15), sf, f.mb.aic)
            ##     eval.df$Model <- paste("Ridge", sel.meth, "AIC MB")
            ##     eval.df$`Number of Features` <- length(mb.aic)
            ##     eval.df$Fold <- k
            ##     ern.conc.df <- rbind(ern.conc.df, eval.df)
            ## } else {
            ##     eval.df <- evaluateModel(cph.mb.aic, test.npn, c(2, 5, 10, 15), sf)
            ##     eval.df$Model <- paste("Ridge", sel.meth, "AIC MB")
            ##     eval.df$`Number of Features` <- length(mb.aic)
            ##     eval.df$Fold <- k
            ##     ern.conc.df <- rbind(ern.conc.df, eval.df)
            ## }
        }

        cv.lasso <- glmnet::cv.glmnet(model.matrix(~., train.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats, paste0(surv.feats,'.risk.nb'), paste0(surv.feats,'.risk.mb'))])[,-1], train.npn[,sf], family='cox')

        plot(cv.lasso)

        eval.df <- evaluateModel(cv.lasso, test.npn, c(2, 5, 10, 15), sf)
        eval.df$Fold <- k
        ern.conc.df <- rbind(ern.conc.df, eval.df)

        
    }
}

ern.conc.df$Target <- factor(ern.conc.df$Target, surv.feats)
ern.conc.df$Time <- factor(ern.conc.df$Time, paste(c(2,5,10,15), 'Years'))

write.csv(ern.conc.df, 'cv_ern_concordance.csv')

gg <- ggplot(ern.conc.df %>% filter(grepl('StARS', Model) | grepl('LASSO', Model), !grepl('Ridge', Model)) %>% group_by(Model, Target, Time), aes(Model, `Uno Concordance`, color=Model)) +
    geom_point(stat='summary', fun.data='mean_sdl') +
    geom_errorbar(stat='summary', width=0.2, fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    geom_hline(yintercept=0.5, color='darkgray', lty=5) + 
    facet_grid(vars(Target), vars(Time), scales='free_y') +
    ## scale_color_brewer(palette='Paired') +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

gg <- ggplot(ern.conc.df %>% filter(Time=='2 Years',grepl('StARS', Model) | grepl('LASSO', Model), !grepl('Ridge', Model)) %>% group_by(Model, Target, Time), aes(Model, `Number of Features`, fill=Model)) +
    geom_bar(stat='summary', fun.data='mean_sdl') +
    geom_errorbar(stat='summary', width=0.2, fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    facet_grid(cols=vars(Target)) +
    ## scale_fill_brewer(palette='Paired') +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg


gg <- ggplot(ern.conc.df %>% filter(Model %in% c('LASSO 1SE', 'StARS BIC MB')) %>% group_by(Model, Target, Time), aes(Time, `Uno Concordance`, color=Model)) +
    geom_point(stat='summary', fun.data='mean_sdl', position=position_dodge2(width=0.5)) +
    geom_errorbar(stat='summary', width=0.5, position=position_dodge2(width=0.5), fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    geom_hline(yintercept=0.5, color='darkgray', lty=5) + 
    facet_grid(cols=vars(Target), scales='free_y') +
    ## stat_compare_means(paired=T) + 
    scale_color_brewer(palette='Set1', labels = c("LASSO Cox", "CausalCoxMGM")) +
    ylab("Uno's Concordance") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

ggsave('uno_concordance_cv_ern_plot.pdf', gg, height=5, width=10)
ggsave('uno_concordance_cv_ern_plot.png', gg, height=5, width=10, dpi=400)

gg <- ggplot(ern.conc.df %>% filter(Time=='2 Years', Model %in% c('LASSO 1SE', 'StARS BIC MB')) %>% group_by(Model, Target, Time), aes(Model, `Number of Features`, fill=Model)) +
    geom_bar(stat='summary', fun.data='mean_sdl') +
    geom_errorbar(stat='summary', width=0.2, fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    facet_grid(cols=vars(Target)) +
    scale_fill_brewer(palette='Set1', labels = c("LASSO Cox", "CausalCoxMGM")) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

ggsave('numFeats_cv_ern_plot.pdf', gg, height=5, width=10)
ggsave('numFeats_cv_ern_plot.png', gg, height=5, width=10, dpi=400)




erp.conc.df <- data.frame()
for (k in 1:5) {
    cvData <- loadGraphAndData('erp', k)
    train.npn <- cvData$train %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))
    test.npn <- cvData$test %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))

    test.npn$ER_STATUS <- 'Positive'

    for (sf in c('DSS','OD','DistantRelapse','LocalRelapse','OS','DRFS','DFS')) {
        for (sel.meth in c('StEPS', 'StARS')) {
            if (sf %in% surv.feats[1:4]) {
                nb <- sort(getNeighbors(cvData$graph[[sel.meth]], sf))
                mb <- sort(setdiff(cvData$graph[[sel.meth]]$markov.blankets[[sf]],
                                   surv.feats))
                mb.bic <- sort(growShrinkMB(cvData$train, sf,
                                            graph=cvData$graph[[sel.meth]],
                                            rank=T, verbose=T))
                ## mb.aic <- sort(growShrinkMB(cvData$train, sf,
                ##                             graph=cvData$graph[[sel.meth]],
                ##                             penalty=2/log(nrow(cvData$train)),
                ##                             rank=T, verbose=T))
            } else {
                ## nb <- paste0(surv.feats[1:(which(sf==surv.feats)-3)], '.risk.nb')
                ## mb <- paste0(surv.feats[1:(which(sf==surv.feats)-3)], '.risk.mb')
                
                nb <- sort(unique(unlist(lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                                                function(x) getNeighbors(cvData$graph[[sel.meth]], x)))))
                mb <- setdiff(
                    sort(
                        unique(
                            unlist(
                                cvData$graph[[sel.meth]]$markov.blankets[surv.feats[1:(which(sf==surv.feats)-3)]]
                            )
                        )
                    ),
                    surv.feats)

                mb.bic <- setdiff(
                    sort(
                        unique(
                            unlist(
                                lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                                       function(x) growShrinkMB(cvData$train, x, graph=cvData$graph[[sel.meth]], rank=T, verbose=T))
                            )
                        )
                    ),
                    surv.feats)

                ## mb.aic <- setdiff(
                ##     sort(
                ##         unique(
                ##             unlist(
                ##                 lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                ##                        function(x) growShrinkMB(cvData$train, x, graph=cvData$graph[[sel.meth]], rank=T, penalty=2/log(nrow(cvData$train)), verbose=T))
                ##             )
                ##         )
                ##     ),
                ##     surv.feats)
            }
            
            f.nb <- as.formula(paste0(sf, ' ~ ', paste(c(1, nb), collapse=' + ')))
            f.mb <- as.formula(paste0(sf, ' ~ ', paste(c(1, mb), collapse=' + ')))
            f.mb.bic <- as.formula(paste0(sf, ' ~ ', paste(c(1, mb.bic), collapse=' + ')))

            ## f.mb.aic <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(c(1, mb.aic), collapse=' + ')))

            cph.nb <- coxph(f.nb, train.npn)
            cph.mb <- coxph(f.mb, train.npn)
            cph.mb.bic <- coxph(f.mb.bic, train.npn)
            ## cph.mb.aic <- coxph(f.mb.aic, train.npn)

            f.nb <- as.formula(paste0('~ ', paste(c(1, nb), collapse=' + ')))
            f.mb <- as.formula(paste0('~ ', paste(c(1, mb), collapse=' + ')))
            f.mb.bic <- as.formula(paste0('~ ', paste(c(1, mb.bic), collapse=' + ')))
            ## f.mb.aic <- as.formula(paste0('~ ', paste(c(1, mb.aic), collapse=' + ')))

            if (length(nb)!=0) {
                ridge.nb <- glmnet::glmnet(as.matrix(model.matrix(f.nb, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            if (length(mb)!=0) {
                ridge.mb <- glmnet::glmnet(as.matrix(model.matrix(f.mb, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            if (length(mb.bic)!=0) {
                ridge.mb.bic <- glmnet::glmnet(as.matrix(model.matrix(f.mb.bic, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            }
            ## if (length(mb.aic)!=0) {
            ##     ridge.mb.aic <- glmnet::glmnet(as.matrix(model.matrix(f.mb.aic, train.npn)), train.npn[,sf], family='cox', alpha=0, lambda=0.5)
            ## }
            
            eval.df <- evaluateModel(cph.nb, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, 'NB')
            eval.df$`Number of Features` <- length(nb)
            eval.df$Fold <- k
            erp.conc.df <- rbind(erp.conc.df, eval.df)

            eval.df <- evaluateModel(cph.mb, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, "MB")
            eval.df$`Number of Features` <- length(mb)
            eval.df$Fold <- k
            erp.conc.df <- rbind(erp.conc.df, eval.df)

            eval.df <- evaluateModel(cph.mb.bic, test.npn, c(2, 5, 10, 15), sf)
            eval.df$Model <- paste(sel.meth, "BIC MB")
            eval.df$`Number of Features` <- length(mb.bic)
            eval.df$Fold <- k
            erp.conc.df <- rbind(erp.conc.df, eval.df)

            ## eval.df <- evaluateModel(cph.mb.aic, test.npn, c(2, 5, 10, 15), sf)
            ## eval.df$Model <- paste(sel.meth, "AIC MB")
            ## eval.df$`Number of Features` <- length(mb.aic)
            ## eval.df$Fold <- k
            ## erp.conc.df <- rbind(erp.conc.df, eval.df)

            if (length(nb)!=0) {
                eval.df <- evaluateModel(ridge.nb, test.npn, c(2, 5, 10, 15), sf, f.nb)
                eval.df$Model <- paste("Ridge", sel.meth, "NB")
                eval.df$`Number of Features` <- length(nb)
                eval.df$Fold <- k
                erp.conc.df <- rbind(erp.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.nb, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "NB")
                eval.df$`Number of Features` <- length(nb)
                eval.df$Fold <- k
                erp.conc.df <- rbind(erp.conc.df, eval.df)
            }

            if (length(mb)!=0) {
                eval.df <- evaluateModel(ridge.mb, test.npn, c(2, 5, 10, 15), sf, f.mb)
                eval.df$Model <- paste("Ridge", sel.meth, "MB")
                eval.df$`Number of Features` <- length(mb)
                eval.df$Fold <- k
                erp.conc.df <- rbind(erp.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.mb, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "MB")
                eval.df$`Number of Features` <- length(mb)
                eval.df$Fold <- k
                erp.conc.df <- rbind(erp.conc.df, eval.df)
            }

            if (length(mb.bic)!=0) {
                eval.df <- evaluateModel(ridge.mb.bic, test.npn, c(2, 5, 10, 15), sf, f.mb.bic)
                eval.df$Model <- paste("Ridge", sel.meth, "BIC MB")
                eval.df$`Number of Features` <- length(mb.bic)
                eval.df$Fold <- k
                erp.conc.df <- rbind(erp.conc.df, eval.df)
            } else {
                eval.df <- evaluateModel(cph.mb.bic, test.npn, c(2, 5, 10, 15), sf)
                eval.df$Model <- paste("Ridge", sel.meth, "BIC MB")
                eval.df$`Number of Features` <- length(mb.bic)
                eval.df$Fold <- k
                erp.conc.df <- rbind(erp.conc.df, eval.df)
            }

            ## if (length(mb.aic)!=0) {
            ##     eval.df <- evaluateModel(ridge.mb.aic, test.npn, c(2, 5, 10, 15), sf, f.mb.aic)
            ##     eval.df$Model <- paste("Ridge", sel.meth, "AIC MB")
            ##     eval.df$`Number of Features` <- length(mb.aic)
            ##     eval.df$Fold <- k
            ##     erp.conc.df <- rbind(erp.conc.df, eval.df)
            ## } else {
            ##     eval.df <- evaluateModel(cph.mb.aic, test.npn, c(2, 5, 10, 15), sf)
            ##     eval.df$Model <- paste("Ridge", sel.meth, "AIC MB")
            ##     eval.df$`Number of Features` <- length(mb.aic)
            ##     eval.df$Fold <- k
            ##     erp.conc.df <- rbind(erp.conc.df, eval.df)
            ## }
        }

        cv.lasso <- glmnet::cv.glmnet(model.matrix(~., train.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats, paste0(surv.feats,'.risk.nb'), paste0(surv.feats,'.risk.mb'))])[,-1], train.npn[,sf], family='cox')

        plot(cv.lasso)

        eval.df <- evaluateModel(cv.lasso, test.npn, c(2, 5, 10, 15), sf)
        eval.df$Fold <- k
        erp.conc.df <- rbind(erp.conc.df, eval.df)

        
    }
}


erp.conc.df$Target <- factor(erp.conc.df$Target, surv.feats)
erp.conc.df$Time <- factor(erp.conc.df$Time, paste(c(2,5,10,15), 'Years'))

write.csv(erp.conc.df, 'cv_erp_concordance.csv')

gg <- ggplot(erp.conc.df %>% filter(grepl('LASSO', Model) | grepl('StARS', Model), !grepl('Ridge', Model)) %>%
             group_by(Model, Target, Time), aes(Model, `Uno Concordance`, color=Model)) +
    geom_point(stat='summary', fun.data='mean_sdl') +
    geom_errorbar(stat='summary', width=0.2, fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    geom_hline(yintercept=0.5, color='darkgray', lty=5) + 
    facet_grid(vars(Target), vars(Time), scales='free_y') +
    ## scale_color_brewer(palette='Paired') +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

gg <- ggplot(erp.conc.df %>% filter(Time=='2 Years', grepl('LASSO', Model) | grepl('StARS', Model), !grepl('Ridge', Model)) %>% group_by(Model, Target, Time), aes(Model, `Number of Features`, fill=Model)) +
    geom_bar(stat='summary', fun.data='mean_sdl') +
    geom_errorbar(stat='summary', width=0.2, fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    facet_grid(cols=vars(Target)) +
    ## scale_fill_brewer(palette='Paired') +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

gg <- ggplot(erp.conc.df %>% filter(Model %in% c('LASSO 1SE', 'StARS BIC MB'), Target!='LocalRelapse') %>% group_by(Model, Target, Time), aes(Time, `Uno Concordance`, color=Model)) +
    geom_point(stat='summary', fun.data='mean_sdl', position=position_dodge2(width=0.5)) +
    geom_errorbar(stat='summary', width=0.5, position=position_dodge2(width=0.5), fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    geom_hline(yintercept=0.5, color='darkgray', lty=5) + 
    facet_grid(cols=vars(Target), scales='free_y') +
    ## stat_compare_means(paired=T) + 
    scale_color_brewer(palette='Set1', labels = c("LASSO Cox", "CausalCoxMGM")) +
    ylab("Uno's Concordance") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

ggsave('uno_concordance_cv_erp_plot.pdf', gg, height=5, width=10)
ggsave('uno_concordance_cv_erp_plot.png', gg, height=5, width=10, dpi=400)

gg <- ggplot(erp.conc.df %>% filter(Time=='2 Years', Model %in% c('LASSO 1SE','StARS BIC MB'), Target!='LocalRelapse') %>% group_by(Model, Target, Time), aes(Model, `Number of Features`, fill=Model)) +
    geom_bar(stat='summary', fun.data='mean_sdl') +
    geom_errorbar(stat='summary', width=0.2, fun.data='mean_sdl', fun.args=list(mult=1.96)) +
    facet_grid(cols=vars(Target)) + ##
    scale_fill_brewer(palette='Set1', labels = c("LASSO Cox", "CausalCoxMGM")) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))
gg

ggsave('numFeats_cv_erp_plot.pdf', gg, height=5, width=10)
ggsave('numFeats_cv_erp_plot.png', gg, height=5, width=10, dpi=400)




cv.conc.ern.dss.npn <- matrix(NA, 5, 6)
for (k in 1:5) {
    cvData <- loadGraphAndData('ern', k)
    train.npn <- cvData$train %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))
    test.npn <- cvData$test %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))
    
    nb <- sort(getNeighbors(cvData$graph, 'DSS'))
    mb <- sort(setdiff(cvData$graph$markov.blankets[['DSS']], surv.feats))

    f.nb <- as.formula(paste0('DSS ~ ', paste(nb, collapse=' + ')))
    f.mb <- as.formula(paste0('DSS ~ ', paste(mb, collapse=' + ')))

    ## cph.nb <- coxph(f.nb, cvData$train)
    ## cph.mb <- coxph(f.mb, cvData$train)

    ## cv.conc.ern.dss[k,1] <- concordance(cph.nb, newdata=cvData$test)$concordance
    ## cv.conc.ern.dss[k,2] <- concordance(cph.mb, newdata=cvData$test)$concordance

    ## ridge.nb <- glmnet::glmnet(model.matrix(~., cvData$train[,nb])[,-1], cvData$train$DSS, family='cox', lambda=0.5, alpha=0)
    ## ridge.mb <- glmnet::glmnet(model.matrix(~., cvData$train[,mb])[,-1], cvData$train$DSS, family='cox', lambda=0.5, alpha=0)

    ## risk.nb <- -predict(ridge.nb, newx=model.matrix(~., cvData$test[,nb])[,-1])
    ## risk.mb <- -predict(ridge.mb, newx=model.matrix(~., cvData$test[,mb])[,-1])

    ## cv.conc.ern.dss[k,3] <- concordance(cvData$test$DSS ~ strata(cvData$test$ER_STATUS) + risk.nb)$concordance
    ## cv.conc.ern.dss[k,4] <- concordance(cvData$test$DSS ~ strata(cvData$test$ER_STATUS) + risk.mb)$concordance

    cph.nb <- coxph(f.nb, train.npn)
    cph.mb <- coxph(f.mb, train.npn)

    cv.conc.ern.dss.npn[k,1] <- concordance(cph.nb, newdata=test.npn)$concordance
    cv.conc.ern.dss.npn[k,2] <- concordance(cph.mb, newdata=test.npn)$concordance

    ridge.nb <- glmnet::glmnet(model.matrix(~., train.npn[,nb])[,-1], train.npn$DSS, family='cox', lambda=0.5, alpha=0)
    ridge.mb <- glmnet::glmnet(model.matrix(~., train.npn[,mb])[,-1], train.npn$DSS, family='cox', lambda=0.5, alpha=0)

    risk.nb <- -predict(ridge.nb, newx=model.matrix(~., test.npn[,nb])[,-1])
    risk.mb <- -predict(ridge.mb, newx=model.matrix(~., test.npn[,mb])[,-1])

    cv.conc.ern.dss.npn[k,3] <- concordance(test.npn$DSS ~ risk.nb)$concordance
    cv.conc.ern.dss.npn[k,4] <- concordance(test.npn$DSS ~ risk.mb)$concordance

    cv.lasso <- glmnet::cv.glmnet(model.matrix(~., train.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], train.npn$DSS, family='cox')

    risk.min <- -predict(cv.lasso, newx=model.matrix(~., test.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], s=cv.lasso$lambda.min)
    risk.1se <- -predict(cv.lasso, newx=model.matrix(~., test.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], s=cv.lasso$lambda.1se)

    cv.conc.ern.dss.npn[k,5] <- concordance(test.npn$DSS ~ risk.1se)$concordance
    cv.conc.ern.dss.npn[k,6] <- concordance(test.npn$DSS ~ risk.min)$concordance

}

colMeans(cv.conc.ern.dss)
apply(cv.conc.dss, 2, sd)

colMeans(cv.conc.ern.dss.npn)
apply(cv.conc.ern.dss.npn, 2, sd)

ern.dss.df <- data.frame(Mean=colMeans(cv.conc.ern.dss.npn),
                             SD=apply(cv.conc.ern.dss.npn, 2, sd),
                             Method=factor(c('Causal Neigbors', 'Causal MB',
                                             'Ridge Causal Neighbors', 'Ridge Causal MB',
                                             'LASSO 1SE', 'LASSO Min'),
                                           levels=c('Causal Neigbors', 'Causal MB',
                                                    'Ridge Causal Neighbors',
                                                    'Ridge Causal MB',
                                                    'LASSO 1SE', 'LASSO Min')),
                             Outcome=rep('DSS Mortailty', 6),
                             ER_STATUS=rep('ER-', 6))

ggplot(ern.dss.df, aes(Method, Mean, color=Method)) +
    geom_point() +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD)) +
    ylab('Concordance') +
    facet_grid(vars(ER_STATUS), vars(Outcome)) +
    scale_color_brewer(palette='Paired') +
    theme_bw()


cv.conc.erp.dss.npn <- matrix(NA, 5, 6)
for (k in 1:5) {
    cvData <- loadGraphAndData('erp', k)
    train.npn <- cvData$train %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))
    test.npn <- cvData$test %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) huge::huge.npn(as.matrix(x)))
    
    nb <- sort(getNeighbors(cvData$graph, 'DSS'))
    mb <- sort(setdiff(cvData$graph$markov.blankets[['DSS']], surv.feats))

    f.nb <- as.formula(paste0('DSS ~ ', paste(nb, collapse=' + ')))
    f.mb <- as.formula(paste0('DSS ~ ', paste(mb, collapse=' + ')))

    ## cph.nb <- coxph(f.nb, cvData$train)
    ## cph.mb <- coxph(f.mb, cvData$train)

    ## cv.conc.erp.dss[k,1] <- concordance(cph.nb, newdata=cvData$test)$concordance
    ## cv.conc.erp.dss[k,2] <- concordance(cph.mb, newdata=cvData$test)$concordance

    ## ridge.nb <- glmnet::glmnet(model.matrix(~., cvData$train[,nb])[,-1], cvData$train$DSS, family='cox', lambda=0.5, alpha=0)
    ## ridge.mb <- glmnet::glmnet(model.matrix(~., cvData$train[,mb])[,-1], cvData$train$DSS, family='cox', lambda=0.5, alpha=0)

    ## risk.nb <- -predict(ridge.nb, newx=model.matrix(~., cvData$test[,nb])[,-1])
    ## risk.mb <- -predict(ridge.mb, newx=model.matrix(~., cvData$test[,mb])[,-1])

    ## cv.conc.erp.dss[k,3] <- concordance(cvData$test$DSS ~ strata(cvData$test$ER_STATUS) + risk.nb)$concordance
    ## cv.conc.erp.dss[k,4] <- concordance(cvData$test$DSS ~ strata(cvData$test$ER_STATUS) + risk.mb)$concordance

    cph.nb <- coxph(f.nb, train.npn)
    cph.mb <- coxph(f.mb, train.npn)

    cv.conc.erp.dss.npn[k,1] <- concordance(cph.nb, newdata=test.npn)$concordance
    cv.conc.erp.dss.npn[k,2] <- concordance(cph.mb, newdata=test.npn)$concordance

    ridge.nb <- glmnet::glmnet(model.matrix(~., train.npn[,nb])[,-1], train.npn$DSS, family='cox', lambda=0.5, alpha=0)
    ridge.mb <- glmnet::glmnet(model.matrix(~., train.npn[,mb])[,-1], train.npn$DSS, family='cox', lambda=0.5, alpha=0)

    risk.nb <- -predict(ridge.nb, newx=model.matrix(~., test.npn[,nb])[,-1])
    risk.mb <- -predict(ridge.mb, newx=model.matrix(~., test.npn[,mb])[,-1])

    cv.conc.erp.dss.npn[k,3] <- concordance(test.npn$DSS ~ risk.nb)$concordance
    cv.conc.erp.dss.npn[k,4] <- concordance(test.npn$DSS ~ risk.mb)$concordance

    cv.lasso <- glmnet::cv.glmnet(model.matrix(~., train.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], train.npn$DSS, family='cox')

    risk.min <- -predict(cv.lasso, newx=model.matrix(~., test.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], s=cv.lasso$lambda.min)
    risk.1se <- -predict(cv.lasso, newx=model.matrix(~., test.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], s=cv.lasso$lambda.1se)

    cv.conc.erp.dss.npn[k,5] <- concordance(test.npn$DSS ~ risk.1se)$concordance
    cv.conc.erp.dss.npn[k,6] <- concordance(test.npn$DSS ~ risk.min)$concordance

}

## colMeans(cv.conc.erp.dss)
## apply(cv.conc.dss, 2, sd)

colMeans(cv.conc.erp.dss.npn)
apply(cv.conc.erp.dss.npn, 2, sd)

erp.dss.df <- data.frame(Mean=colMeans(cv.conc.erp.dss.npn),
                             SD=apply(cv.conc.erp.dss.npn, 2, sd),
                             Method=factor(c('Causal Neigbors', 'Causal MB',
                                             'Ridge Causal Neighbors', 'Ridge Causal MB',
                                             'LASSO 1SE', 'LASSO Min'),
                                           levels=c('Causal Neigbors', 'Causal MB',
                                                    'Ridge Causal Neighbors',
                                                    'Ridge Causal MB',
                                                    'LASSO 1SE', 'LASSO Min')),
                             Outcome=rep('DSS Mortailty', 6),
                             ER_STATUS=rep('ER+', 6))

ggplot(erp.dss.df, aes(Method, Mean, color=Method)) +
    geom_point() +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD)) +
    ylab('Concordance') +
    facet_grid(vars(ER_STATUS), vars(Outcome)) +
    scale_color_brewer(palette='Paired') +
    theme_bw()


gg <- ggplot(rbind(dss.df, erp.dss.df, ern.dss.df),
       aes(Method, Mean, color=Method)) +
    geom_point() +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.5) +
    ylab('Concordance') +
    facet_grid(cols=vars(ER_STATUS), rows=vars(Outcome)) +
    scale_color_brewer(palette='Paired') +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust=NULL, hjust=1))
gg

ggsave('concordance_cv_plot.pdf', gg, height=5, width=10)
ggsave('concordance_cv_plot.png', gg, height=5, width=10, dpi=400)





metabric <- read.csv('data/metabric.full.csv')

if (any(c('DSS.time', 'DSS.status') %in% colnames(metabric))) {
   metabric$DSS <- Surv(metabric$DSS.time, metabric$DSS.status)
   attr(metabric$DSS, 'strata') <- metabric$ER_STATUS
}
if (any(c('OD.time', 'OD.status') %in% colnames(metabric))) {
   metabric$OD <- Surv(metabric$OD.time, metabric$OD.status)
   attr(metabric$OD, 'strata') <- metabric$ER_STATUS
}
if (any(c('DistantRelapse.time', 'DistantRelapse.status') %in% colnames(metabric))) {
   metabric$DistantRelapse <- Surv(metabric$DistantRelapse.time, metabric$DistantRelapse.status)
   attr(metabric$DistantRelapse, 'strata') <- metabric$ER_STATUS   
}
if (any(c('LocalRelapse.time', 'LocalRelapse.status') %in% colnames(metabric))) {
   metabric$LocalRelapse <- Surv(metabric$LocalRelapse.time, metabric$LocalRelapse.status)
   attr(metabric$LocalRelapse, 'strata') <- metabric$ER_STATUS
}


metabric <- metabric[,!colnames(metabric) %in% c('DSS.time', 'DSS.status',
                                                 'OD.time', 'OD.status',
                                                 'DistantRelapse.time',
                                                 'DistantRelapse.status',
                                                 'LocalRelapse.time',
                                                 'LocalRelapse.status')]

metabric$OS <- Surv(pmin(metabric$DSS[,1], metabric$OD[,1]), pmax(metabric$DSS[,2], metabric$OD[,2]))
attr(metabric$OS, 'strata') <- metabric$ER_STATUS

metabric$DRFS <- Surv(pmin(metabric$DSS[,1], metabric$OD[,1], metabric$DistantRelapse[,1]), pmax(metabric$DSS[,2], metabric$OD[,2], metabric$DistantRelapse[,2]))
attr(metabric$DRFS, 'strata') <- metabric$ER_STATUS

metabric$DFS <- Surv(pmin(metabric$DSS[,1], metabric$OD[,1], metabric$DistantRelapse[,1], metabric$LocalRelapse[,1]), pmax(metabric$DSS[,2], metabric$OD[,2], metabric$DistantRelapse[,2], metabric$LocalRelapse[,2]))
attr(metabric$DFS, 'strata') <- metabric$ER_STATUS

metabric$LYMPH_NODE_STATUS <- factor(metabric$LYMPH_NODE_STATUS, levels=c('NodeNegative','1to3','4toX'))
metabric$CNA_MIS12 <- factor(metabric$CNA_MIS12, levels=c('NEUTRAL','LOSS','GAIN'))

g <- loadGraph('out/mfcimgmStEPS.metabric.all..rank.full.txt')

surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse')

ggsurvplot(survfit(DSS ~ ER_STATUS, metabric), pval=T, risk.table=T)
ggsurvplot(survfit(OD ~ ER_STATUS, metabric), pval=T, risk.table=T)
ggsurvplot(survfit(OS ~ ER_STATUS, metabric), pval=T, risk.table=T)
ggsurvplot(survfit(DistantRelapse ~ ER_STATUS, metabric), pval=T, risk.table=T)
ggsurvplot(survfit(DRFS ~ ER_STATUS, metabric), pval=T, risk.table=T)
ggsurvplot(survfit(LocalRelapse ~ ER_STATUS, metabric), pval=T, risk.table=T)
ggsurvplot(survfit(DFS ~ ER_STATUS, metabric), pval=T, risk.table=T)

metabric.zscore <- metabric %>%
    mutate_if(function(x) !is.Surv(x) & is.numeric(x),
              function(x) scale(x)[,1])

mod.mat.dss <- model.matrix(## ~ RNA_GSTO2 + RNA_SERPINA6 + RNA_S100P + T.cells.gamma.delta + LYMPH_NODE_STATUS + TUMOR_SIZE,
    ~ RNA_IGFBP5 + RNA_NQO1 + RNA_S100P + RNA_UBE2C + LYMPH_NODE_STATUS + TUMOR_SIZE,
                            metabric.npn)[,-1]

p.lymph <- table(metabric.zscore$LYMPH_NODE_STATUS) / nrow(metabric)

mod.mat.dss[,c("LYMPH_NODE_STATUS1to3", "LYMPH_NODE_STATUS4toX")] <- mod.mat.dss[,c("LYMPH_NODE_STATUS1to3", "LYMPH_NODE_STATUS4toX")] / sqrt(sum(p.lymph * (1-p.lymph)))

best.bic <- Inf
for (k in 1:10) {
    mcl.model <- Mclust(mod.mat.dss, k)
    mcl.model$parameters
    mcl.model$BIC
    if (mcl.model$bic < best.bic) {
        best.bic <- mcl.model$bic
        best.model <- mcl.model
        best.class <- mcl.model$classification
        best.k <- k
    }
}

best.k
best.bic
best.model$parameters
best.class

pc.dss <- prcomp(t(mod.mat.dss))

pc.dss.df <- data.frame(pc.dss$rotation,
                        metabric.zscore[,c('LYMPH_NODE_STATUS', 'ER_STATUS', 'HER2_STATUS', 'PR_STATUS', 'TUMOR_SIZE')],
                        risk=risk.dss,
                        class=best.class)

ggplot(pc.dss.df, aes(PC1, PC2, color=ER_STATUS)) +
    geom_point() +
    theme_bw()

ggplot(pc.dss.df, aes(PC1, PC2, color=HER2_STATUS)) +
    geom_point() +
    theme_bw()

ggplot(pc.dss.df, aes(PC1, PC2, color=PR_STATUS)) +
    geom_point() +
    theme_bw()

ggplot(pc.dss.df, aes(PC1, PC2, color=LYMPH_NODE_STATUS)) +
    geom_point(size=5) +
    theme_bw()

ggplot(pc.dss.df, aes(PC1, PC2, color=factor(class))) +
    geom_point(size=5) +
    theme_bw()

ggplot(pc.dss.df, aes(PC1, PC2, color=risk)) +
    geom_point(size=5) +
    theme_bw()

ggplot(pc.dss.df, aes(PC1, PC2, color=TUMOR_SIZE)) +
    geom_point(size=5) +
    theme_bw()

metabric.npn <-  metabric %>%
    mutate_if(function(x) !is.Surv(x) & is.numeric(x),
              function(x) huge::huge.npn(as.matrix(x)))

full.npn <- fullData$train %>% mutate_if(function(x) !is.Surv(x) & is.numeric(x), function(x) as.vector(huge::huge.npn(as.matrix(x))))

f.dss <- DSS ~ strata(ER_STATUS) + RNA_IGFBP5 + RNA_NQO1 + RNA_S100P + RNA_UBE2C + LYMPH_NODE_STATUS + TUMOR_SIZE
f.dss.mb <- as.formula(paste0('DSS ~ strata(ER_STATUS) + ', paste(setdiff(g$markov.blankets[['DSS']], surv.feats), collapse=' + ')))
f.od <- OD ~ strata(ER_STATUS) + AGE_AT_DIAGNOSIS
f.os <- OS ~ strata(ER_STATUS) + RNA_IGFBP5 + RNA_NQO1 + RNA_S100P + RNA_UBE2C + LYMPH_NODE_STATUS + TUMOR_SIZE + AGE_AT_DIAGNOSIS
f.dr <- DistantRelapse ~ strata(ER_STATUS) + CNA_MIS12 + RNA_MMP11 + RNA_S100P + RNA_SERPINA1 + RNA_UBE2C + LYMPH_NODE_STATUS + TUMOR_SIZE
f.dr.mb <- as.formula(paste0('DistantRelapse ~ strata(ER_STATUS) + ', paste(setdiff(g$markov.blankets[['DistantRelapse']], surv.feats), collapse=' + ')))
f.drfs <- DRFS ~ strata(ER_STATUS) + RNA_IGFBP5 + RNA_NQO1 + RNA_S100P + RNA_UBE2C + LYMPH_NODE_STATUS + TUMOR_SIZE + AGE_AT_DIAGNOSIS + CNA_MIS12 + RNA_MMP11 + RNA_SERPINA1
f.dfs <- DFS ~ strata(ER_STATUS) + RNA_IGFBP5 + RNA_NQO1 + RNA_S100P + RNA_UBE2C + LYMPH_NODE_STATUS + TUMOR_SIZE + AGE_AT_DIAGNOSIS + CNA_MIS12 + RNA_MMP11 + RNA_SERPINA1

cphridge.dss <- glmnet(x=model.matrix(~., metabric.npn[,setdiff(g$markov.blankets[['DSS']],
                                                                surv.feats)])[,-1],
                       y=stratifySurv(metabric.npn[,'DSS'], metabric.npn$ER_STATUS),
                       alpha=0, lambda=0.5, family='cox')
cph.dss <- coxph(f.dss.mb, metabric.npn)
coef(cph.dss)
summary(cph.dss)
concordance(cph.dss, ymax=5 * 365.25)
cox.zph(cph.dss)
risk.dss.mb <- -predict(cph.dss)
risk.dss.mb.ridge <- -predict(cphridge.dss,
                             newx=model.matrix(~., metabric.npn[,setdiff(g$markov.blankets[['DSS']],
                                                                         surv.feats)])[,-1])
concordance(metabric.npn$DSS ~ strata(metabric.npn$ER_STATUS) + risk.dss.mb, ymax=5 * 365.25)
concordance(metabric.npn$DSS ~ strata(metabric.npn$ER_STATUS) + risk.dss.mb.ridge, ymax=5 * 365.25)
ggsurvplot(survfit(DSS ~ ER_STATUS + quantcut(risk.dss,2), metabric), pval=T, risk.table=T)
## ggsurvplot(survfit(DSS ~ best.class, metabric), pval=T, risk.table=T)

cph.od <- coxph(f.od, metabric.npn)
summary(cph.od)
concordance(cph.od, ymax=5 * 365.25)
cox.zph(cph.od)
risk.od <- predict(cph.od)
ggsurvplot(survfit(OD ~ ER_STATUS + quantcut(risk.od,2), metabric), pval=T, risk.table=T)

cph.os <- coxph(f.os, metabric.npn)
summary(cph.os)
concordance(cph.os, ymax=5 * 365.25)
cox.zph(cph.os)
risk.os <- predict(cph.os)
ggsurvplot(survfit(OS ~ ER_STATUS + quantcut(risk.os,2), metabric), pval=T, risk.table=T)


cph.dr <- coxph(f.dr.mb, metabric.npn)
summary(cph.dr)
concordance(cph.dr, ymax=5 * 365.25)
cox.zph(cph.dr)

cph.drfs <- coxph(f.drfs, metabric.npn)
summary(cph.drfs)
concordance(cph.drfs, ymax=5 * 365.25)
cox.zph(cph.drfs)

cph.dfs <- coxph(f.dfs, metabric.npn)
summary(cph.dfs)
concordance(cph.dfs, ymax=5 * 365.25)
cox.zph(cph.dfs)

cna.gene.list <- readRDS('../metabric_cna_highCor_gene_groups.RDS')
cna.ego.list <- readRDS('../metabric_cna_highCor_gene_groups_ora.RDS')

library(clusterProfiler)
head(cna.ego.list[['MIS12']])
cna.ego.list[['MIS12']][grep('BREAST', cna.ego.list[['MIS12']][,'ID']),]
cna.gene.list[['MIS12']]


metabric.ern <- read.csv('data/metabric.ern.full.csv')

if (any(c('DSS.time', 'DSS.status') %in% colnames(metabric.ern))) {
   metabric.ern$DSS <- Surv(metabric.ern$DSS.time, metabric.ern$DSS.status)
   attr(metabric.ern$DSS, 'strata') <- metabric.ern$ER_STATUS
}
if (any(c('OD.time', 'OD.status') %in% colnames(metabric.ern))) {
   metabric.ern$OD <- Surv(metabric.ern$OD.time, metabric.ern$OD.status)
   attr(metabric.ern$OD, 'strata') <- metabric.ern$ER_STATUS
}
if (any(c('DistantRelapse.time', 'DistantRelapse.status') %in% colnames(metabric.ern))) {
   metabric.ern$DistantRelapse <- Surv(metabric.ern$DistantRelapse.time, metabric.ern$DistantRelapse.status)
   attr(metabric.ern$DistantRelapse, 'strata') <- metabric.ern$ER_STATUS   
}
if (any(c('LocalRelapse.time', 'LocalRelapse.status') %in% colnames(metabric.ern))) {
   metabric.ern$LocalRelapse <- Surv(metabric.ern$LocalRelapse.time, metabric.ern$LocalRelapse.status)
   attr(metabric.ern$LocalRelapse, 'strata') <- metabric.ern$ER_STATUS
}


metabric.ern <- metabric.ern[,!colnames(metabric.ern) %in% c('DSS.time', 'DSS.status',
                                                             'OD.time', 'OD.status',
                                                             'DistantRelapse.time',
                                                             'DistantRelapse.status',
                                                             'LocalRelapse.time',
                                                             'LocalRelapse.status')]

metabric.ern$OS <- Surv(pmin(metabric.ern$DSS[,1], metabric.ern$OD[,1]), pmax(metabric.ern$DSS[,2], metabric.ern$OD[,2]))
attr(metabric.ern$OS, 'strata') <- metabric.ern$ER_STATUS

metabric.ern$DRFS <- Surv(pmin(metabric.ern$DSS[,1], metabric.ern$OD[,1], metabric.ern$DistantRelapse[,1]), pmax(metabric.ern$DSS[,2], metabric.ern$OD[,2], metabric.ern$DistantRelapse[,2]))
attr(metabric.ern$DRFS, 'strata') <- metabric.ern$ER_STATUS

metabric.ern$DFS <- Surv(pmin(metabric.ern$DSS[,1], metabric.ern$OD[,1], metabric.ern$DistantRelapse[,1], metabric.ern$LocalRelapse[,1]), pmax(metabric.ern$DSS[,2], metabric.ern$OD[,2], metabric.ern$DistantRelapse[,2], metabric.ern$LocalRelapse[,2]))
attr(metabric.ern$DFS, 'strata') <- metabric.ern$ER_STATUS

metabric.ern$LYMPH_NODE_STATUS <- factor(metabric.ern$LYMPH_NODE_STATUS, levels=c('NodeNegative','1to3','4toX'))
metabric.ern$CNA_MIS12 <- factor(metabric.ern$CNA_MIS12, levels=c('NEUTRAL','LOSS','GAIN'))


ggsurvplot(survfit(DSS ~ 1, metabric.ern), pval=T, risk.table=T)
ggsurvplot(survfit(OD ~ 1, metabric.ern), pval=T, risk.table=T)
ggsurvplot(survfit(OS ~ 1, metabric.ern), pval=T, risk.table=T)
ggsurvplot(survfit(DistantRelapse ~ 1, metabric.ern), pval=T, risk.table=T)
ggsurvplot(survfit(DRFS ~ 1, metabric.ern), pval=T, risk.table=T)
ggsurvplot(survfit(LocalRelapse ~ 1, metabric.ern), pval=T, risk.table=T)
ggsurvplot(survfit(DFS ~ 1, metabric.ern), pval=T, risk.table=T)

metabric.ern.npn <-  metabric.ern %>%
    mutate_if(function(x) !is.Surv(x) & is.numeric(x),
              function(x) huge::huge.npn(as.matrix(x)))

f.ern.dss <- DSS ~ RNA_GSTO2 + RNA_SERPINA6 + RNA_S100P + T.cells.gamma.delta + LYMPH_NODE_STATUS + TUMOR_SIZE
f.ern.od <- OD ~ AGE_AT_DIAGNOSIS
f.ern.os <- OS ~ RNA_GSTO2 + RNA_SERPINA6 + RNA_S100P + T.cells.gamma.delta + LYMPH_NODE_STATUS + TUMOR_SIZE + AGE_AT_DIAGNOSIS
f.ern.dr <- DistantRelapse ~ RNA_GSTO2 + RNA_SERPINA6 + RNA_S100P + T.cells.gamma.delta + LYMPH_NODE_STATUS + RNA_C19orf33 + RNA_DCD + RNA_HLA.DRB4
f.ern.drfs <- DRFS ~ RNA_GSTO2 + RNA_SERPINA6 + RNA_S100P + T.cells.gamma.delta + LYMPH_NODE_STATUS + TUMOR_SIZE + AGE_AT_DIAGNOSIS + RNA_C19orf33 + RNA_DCD + RNA_HLA.DRB4
f.ern.lr <- LocalRelapse ~ LYMPH_NODE_STATUS + RNA_DCD + RNA_CHAD
f.ern.dfs <- DFS ~ RNA_GSTO2 + RNA_SERPINA6 + RNA_S100P + T.cells.gamma.delta + LYMPH_NODE_STATUS + TUMOR_SIZE + AGE_AT_DIAGNOSIS + RNA_C19orf33 + RNA_DCD + RNA_HLA.DRB4 + RNA_CHAD

cph.ern.dss <- coxph(f.ern.dss, metabric.ern.npn)
summary(cph.ern.dss)
concordance(cph.ern.dss, ymax=5 * 365.25)
cox.zph(cph.ern.dss)
risk.ern.dss <- predict(cph.ern.dss)

cph.ern.od <- coxph(f.ern.od, metabric.ern.npn)
summary(cph.ern.od)
concordance(cph.ern.od, ymax=5 * 365.25)
cox.zph(cph.ern.od)

cph.ern.os <- coxph(f.ern.os, metabric.ern.npn)
summary(cph.ern.os)
concordance(cph.ern.os, ymax=5 * 365.25)
cox.zph(cph.ern.os)

cph.ern.dr <- coxph(f.ern.dr, metabric.ern.npn)
summary(cph.ern.dr)
concordance(cph.ern.dr, ymax=5 * 365.25)
cox.zph(cph.ern.dr)

cph.ern.drfs <- coxph(f.ern.drfs, metabric.ern.npn)
summary(cph.ern.drfs)
concordance(cph.ern.drfs, ymax=5 * 365.25)
cox.zph(cph.ern.drfs)

cph.ern.lr <- coxph(f.ern.lr, metabric.ern.npn)
summary(cph.ern.lr)
concordance(cph.ern.lr, ymax=5 * 365.25)
cox.zph(cph.ern.lr)

cph.ern.dfs <- coxph(f.ern.dfs, metabric.ern.npn)
summary(cph.ern.dfs)
concordance(cph.ern.dfs, ymax=5 * 365.25)
cox.zph(cph.ern.dfs)
