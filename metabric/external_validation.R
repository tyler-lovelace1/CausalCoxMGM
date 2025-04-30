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

trainLassoModels <- function(train, feats, erstatus) {
    feats <- setdiff(feats, 'ER_STATUS')
    train$LYMPH_NODE_STATUS <- factor(train$LYMPH_NODE_STATUS,
                                     levels=c('Negative','Positive'))
    models <- list()
    for (sf in surv.feats[c(1,5:7)]) {
        models[[sf]] <- cv.glmnet(x=model.matrix(~.,train[,feats])[,-1],
                                  y=train[,sf], family='cox',
                                  lambda.min.ratio=0.01)
        plot(models[[sf]])
    }
    return(models)
}


evaluateCausalModels <- function(model.list, ext.data.list, times) {
    surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')
    eval.df <- data.frame()
    for (ds.name in names(ext.data.list)) {
        test <- ext.data.list[[ds.name]]
        if (!'AGE_AT_DIAGNOSIS' %in% colnames(test)) {
            test$AGE_AT_DIAGNOSIS <- 0
        }

        if (ds.name=='GSE11121') {
            test[,'DRFS'][,1] <- 30.437 * test[,'DRFS'][,1]
        }

        test$LYMPH_NODE_STATUS <- factor(test$LYMPH_NODE_STATUS,
                                         levels=c('Negative','Positive'))
        
        test$ER_STATUS <- factor(test$ER_STATUS,
                                 levels=c('Negative','Positive'))
        
        erPos <- any(test$ER_STATUS=='Positive')
        erNeg <- any(test$ER_STATUS=='Negative')

        for (model in c('Stratified', 'Split')) {
            for (sf in surv.feats) {
                if (model == 'Stratified') {
                    erstatus <- 'all'
                    test[,paste('risk',sf,sep='.')] <- predict(model.list[[erstatus]][[sf]],
                                                               test)
                } else {
                    test[,paste('risk',sf,sep='.')] <- rep(0, nrow(test))
                    erstatus <- 'erp'
                    if (erPos) {
                        test[test$ER_STATUS=='Positive',paste('risk',sf,sep='.')] <-
                            predict(model.list[[erstatus]][[sf]],
                                    test[test$ER_STATUS=='Positive',])
                    }
                    erstatus <- 'ern'
                    if (erPos) {
                        test[test$ER_STATUS=='Negative',paste('risk',sf,sep='.')] <-
                            predict(model.list[[erstatus]][[sf]],
                                    test[test$ER_STATUS=='Negative',])
                    }
                }
                
                if (sf %in% colnames(test)) {
                    if (sf %in% c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse')) {
                        mb.all <- names(model.list[['all']][[sf]]$coefficients)
                        mb.erp <- names(model.list[['erp']][[sf]]$coefficients)
                        mb.ern <- names(model.list[['ern']][[sf]]$coefficients)
                    } else {
                        sf.splits <- surv.feats[1:(which(surv.feats==sf)-3)]
                        mb.all <- unique(
                            unlist(
                                lapply(
                                    sf.splits,
                                    function(x) {
                                        names(model.list[['all']][[x]]$coefficients)
                                    })))
                        mb.erp <- unique(
                            unlist(
                                lapply(
                                    sf.splits,
                                    function(x) {
                                        names(model.list[['erp']][[x]]$coefficients)
                                    })))
                        mb.ern <- unique(
                            unlist(
                                lapply(
                                    sf.splits,
                                    function(x) {
                                        names(model.list[['ern']][[x]]$coefficients)
                                    })))
                    }
                    for (t in times) {
                        uno <- concordance(
                            as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                             paste('risk',sf,sep='.'))),
                            test,
                            ymax=365.25*t,
                            timewt='n/G2')
                        
                        harrell <- concordance(
                            as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                             paste('risk',sf,sep='.'))),
                            test,
                            ymax=365.25*t)
                        
                        eval.df <- rbind(
                            eval.df,
                            data.frame(
                                `Uno Concordance` = 1-uno$concordance,
                                `Harrell Concordance` = 1-harrell$concordance,
                                `Uno SE` = sqrt(uno$var),
                                `Harrell SE` = sqrt(harrell$var),
                                Model = paste(model, 'Causal'),
                                `ER Status` = 'All',
                                Time = paste(t, 'Years'),
                                Outcome = sf,
                                Dataset = ds.name,
                                `Number of Features`=ifelse(
                                    model=='Split',
                                    length(union(mb.erp, mb.ern)),
                                    length(mb.all)),
                                check.names=F))

                        if (erPos) {
                            uno <- concordance(
                                as.formula(paste(sf, '~',
                                                 paste('risk',sf,sep='.'))),
                                test[test$ER_STATUS=='Positive',],
                                ymax=365.25*t,
                                timewt='n/G2')
                            
                            harrell <- concordance(
                                as.formula(paste(sf, '~',
                                                 paste('risk',sf,sep='.'))),
                                test[test$ER_STATUS=='Positive',],
                                ymax=365.25*t)
                            
                            eval.df <- rbind(
                                eval.df,
                                data.frame(
                                    `Uno Concordance` = 1-uno$concordance,
                                    `Harrell Concordance` = 1-harrell$concordance,
                                    `Uno SE` = sqrt(uno$var),
                                    `Harrell SE` = sqrt(harrell$var),
                                    Model = paste(model, 'Causal'),
                                    `ER Status` = 'ER+',
                                    Time = paste(t, 'Years'),
                                    Outcome = sf,
                                    Dataset = ds.name,
                                    `Number of Features`=ifelse(
                                        model=='Split',
                                        length(mb.erp),
                                        length(mb.all)),
                                    check.names=F))

                        }

                        if (erNeg) {
                            uno <- concordance(
                                as.formula(paste(sf, '~',
                                                 paste('risk',sf,sep='.'))),
                                test[test$ER_STATUS=='Negative',],
                                ymax=365.25*t,
                                timewt='n/G2')
                            
                            harrell <- concordance(
                                as.formula(paste(sf, '~',
                                                 paste('risk',sf,sep='.'))),
                                test[test$ER_STATUS=='Negative',],
                                ymax=365.25*t)
                            
                            eval.df <- rbind(
                                eval.df,
                                data.frame(
                                    `Uno Concordance` = 1-uno$concordance,
                                    `Harrell Concordance` = 1-harrell$concordance,
                                    `Uno SE` = sqrt(uno$var),
                                    `Harrell SE` = sqrt(harrell$var),
                                    Model = paste(model, 'Causal'),
                                    `ER Status` = 'ER-',
                                    Time = paste(t, 'Years'),
                                    Outcome = sf,
                                    Dataset = ds.name,
                                    `Number of Features`=ifelse(
                                        model=='Split',
                                        length(mb.ern),
                                        length(mb.all)),
                                    check.names=F))
                        }
                    }
                }
            }
        }
    }
    return(eval.df)
}


evaluateLassoModels <- function(model.list, ext.data.list, feats, times, select='Min') {
    feats <- setdiff(feats, 'ER_STATUS')
    eval.df <- data.frame()
    for (ds.name in names(ext.data.list)) {
        test <- ext.data.list[[ds.name]]
        if (!'AGE_AT_DIAGNOSIS' %in% colnames(test)) {
            test$AGE_AT_DIAGNOSIS <- 0
        }

        if (ds.name=='GSE11121') {
            test[,'DRFS'][,1] <- 30.437 * test[,'DRFS'][,1]
        }

        test$LYMPH_NODE_STATUS <- factor(test$LYMPH_NODE_STATUS,
                                         levels=c('Negative','Positive'))
        
        test$ER_STATUS <- factor(test$ER_STATUS,
                                 levels=c('Negative','Positive'))
        
        erPos <- any(test$ER_STATUS=='Positive')
        erNeg <- any(test$ER_STATUS=='Negative')

        for (model in c('Stratified', 'Split')) {
            for (sf in surv.feats[c(1,5:7)]) {
                if (model == 'Stratified') {
                    erstatus <- 'all'
                    test[,paste('risk',sf,sep='.')] <- predict(
                        model.list[[erstatus]][[sf]],
                        newx=model.matrix(~., test[,feats])[,-1],
                        s=ifelse(select=='Min',
                                 model.list[[erstatus]][[sf]]$lambda.min,
                                 model.list[[erstatus]][[sf]]$lambda.1se))
                    featSet <- dimnames(coef(model.list[[erstatus]][[sf]]))[[1]][as.vector(coef(model.list[[erstatus]][[sf]], s=ifelse(select=='Min', model.list[[erstatus]][[sf]]$lambda.min, model.list[[erstatus]][[sf]]$lambda.1se)))!=0]
                } else {
                    test[,paste('risk',sf,sep='.')] <- rep(0, nrow(test))
                    erstatus <- 'erp'
                    if (erPos) {
                        test[test$ER_STATUS=='Positive',paste('risk',sf,sep='.')] <-
                            predict(
                                model.list[[erstatus]][[sf]],
                                newx=model.matrix(
                                    ~.,test[test$ER_STATUS=='Positive',feats])[,-1],
                                s=ifelse(select=='Min',
                                         model.list[[erstatus]][[sf]]$lambda.min,
                                         model.list[[erstatus]][[sf]]$lambda.1se))
                    }
                    featSet.erp <- dimnames(coef(model.list[[erstatus]][[sf]]))[[1]][as.vector(coef(model.list[[erstatus]][[sf]], s=ifelse(select=='Min', model.list[[erstatus]][[sf]]$lambda.min, model.list[[erstatus]][[sf]]$lambda.1se)))!=0]
                    erstatus <- 'ern'
                    if (erPos) {
                        test[test$ER_STATUS=='Negative',paste('risk',sf,sep='.')] <-
                            predict(model.list[[erstatus]][[sf]],
                                    newx=model.matrix(
                                        ~.,test[test$ER_STATUS=='Negative',feats])[,-1],
                                    s=ifelse(select=='Min',
                                             model.list[[erstatus]][[sf]]$lambda.min,
                                             model.list[[erstatus]][[sf]]$lambda.1se))
                    }
                    featSet.ern <- dimnames(coef(model.list[[erstatus]][[sf]]))[[1]][as.vector(coef(model.list[[erstatus]][[sf]], s=ifelse(select=='Min', model.list[[erstatus]][[sf]]$lambda.min, model.list[[erstatus]][[sf]]$lambda.1se)))!=0]
                }
                
                if (sf %in% colnames(test)) {
                    for (t in times) {
                        uno <- concordance(
                            as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                             paste('risk',sf,sep='.'))),
                            test,
                            ymax=365.25*t,
                            timewt='n/G2')
                        
                        harrell <- concordance(
                            as.formula(paste(sf, '~ strata(ER_STATUS) +',
                                             paste('risk',sf,sep='.'))),
                            test,
                            ymax=365.25*t)
                        
                        eval.df <- rbind(
                            eval.df,
                            data.frame(
                                `Uno Concordance` = 1-uno$concordance,
                                `Harrell Concordance` = 1-harrell$concordance,
                                `Uno SE` = sqrt(uno$var),
                                `Harrell SE` = sqrt(harrell$var),
                                Model = paste(model, 'LASSO', select),
                                `ER Status` = 'All',
                                Time = paste(t, 'Years'),
                                Outcome = sf,
                                Dataset = ds.name,
                                `Number of Features`=ifelse(
                                    model=='Split',
                                    length(union(featSet.erp, featSet.ern)),
                                    length(featSet)),
                                check.names=F))

                        if (erPos) {
                            uno <- concordance(
                                as.formula(paste(sf, '~',
                                                 paste('risk',sf,sep='.'))),
                                test[test$ER_STATUS=='Positive',],
                                ymax=365.25*t,
                                timewt='n/G2')
                            
                            harrell <- concordance(
                                as.formula(paste(sf, '~',
                                                 paste('risk',sf,sep='.'))),
                                test[test$ER_STATUS=='Positive',],
                                ymax=365.25*t)
                            
                            eval.df <- rbind(
                                eval.df,
                                data.frame(
                                    `Uno Concordance` = 1-uno$concordance,
                                    `Harrell Concordance` = 1-harrell$concordance,
                                    `Uno SE` = sqrt(uno$var),
                                    `Harrell SE` = sqrt(harrell$var),
                                    Model = paste(model, 'LASSO', select),
                                    `ER Status` = 'ER+',
                                    Time = paste(t, 'Years'),
                                    Outcome = sf,
                                    Dataset = ds.name,
                                    `Number of Features`=ifelse(
                                        model=='Split',
                                        length(featSet.erp),
                                        length(featSet)),
                                    check.names=F))

                        }

                        if (erNeg) {
                            uno <- concordance(
                                as.formula(paste(sf, '~',
                                                 paste('risk',sf,sep='.'))),
                                test[test$ER_STATUS=='Negative',],
                                ymax=365.25*t,
                                timewt='n/G2')
                            
                            harrell <- concordance(
                                as.formula(paste(sf, '~',
                                                 paste('risk',sf,sep='.'))),
                                test[test$ER_STATUS=='Negative',],
                                ymax=365.25*t)
                            
                            eval.df <- rbind(
                                eval.df,
                                data.frame(
                                    `Uno Concordance` = 1-uno$concordance,
                                    `Harrell Concordance` = 1-harrell$concordance,
                                    `Uno SE` = sqrt(uno$var),
                                    `Harrell SE` = sqrt(harrell$var),
                                    Model = paste(model, 'LASSO', select),
                                    `ER Status` = 'ER-',
                                    Time = paste(t, 'Years'),
                                    Outcome = sf,
                                    Dataset = ds.name,
                                    `Number of Features`=ifelse(
                                        model=='Split',
                                        length(featSet.ern),
                                        length(featSet)),
                                    check.names=F))
                        }
                    }
                }
            }
        }
    }
    return(eval.df)
}


loadGraphAndData <- function(erstatus='all', k=-1) {

    if (k < 0) {
        cvString <- 'full'
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

    if (erstatus=='all') {
        train$DSS <- stratifySurv(train$DSS, train$ER_STATUS)
        ## train$OD <- stratifySurv(train$OD, train$ER_STATUS)
        train$DistantRelapse <- stratifySurv(train$DistantRelapse, train$ER_STATUS)
        train$LocalRelapse <- stratifySurv(train$LocalRelapse, train$ER_STATUS)
        train$OS <- stratifySurv(train$OS, train$ER_STATUS)
        train$DRFS <- stratifySurv(train$DRFS, train$ER_STATUS)
        train$DFS <- stratifySurv(train$DFS, train$ER_STATUS)
        if (cvString != 'full') {
            test$DSS <- stratifySurv(test$DSS, test$ER_STATUS)
            ## test$OD <- stratifySurv(test$OD, test$ER_STATUS)
            test$DistantRelapse <- stratifySurv(test$DistantRelapse, test$ER_STATUS)
            test$LocalRelapse <- stratifySurv(test$LocalRelapse, test$ER_STATUS)
            test$OS <- stratifySurv(test$OS, test$ER_STATUS)
            test$DRFS <- stratifySurv(test$DRFS, test$ER_STATUS)
            test$DFS <- stratifySurv(test$DFS, test$ER_STATUS)
        }
    }

    train <- train[,!colnames(train) %in% c('DSS.time', 'DSS.status',
                                            'OD.time', 'OD.status',
                                            'DistantRelapse.time', 'DistantRelapse.status',
                                            'LocalRelapse.time', 'LocalRelapse.status',
                                            'OS.time', 'OS.status',
                                            'DRFS.time', 'DRFS.status',
                                            'DFS.time', 'DFS.status')]
    if (erstatus=='ern') {
        train[train$GRADE!='G3','GRADE'] <- 'G1.G2'
        train[!train$ONCOTREE_CODE!='IDC','ONCOTREE_CODE'] <- 'BREAST'
    }

    train <- train %>% mutate_if(is.character, factor)

    if (cvString != 'full') {
        test <- test[,!colnames(test) %in% c('DSS.time', 'DSS.status',
                                             'OD.time', 'OD.status',
                                             'DistantRelapse.time', 'DistantRelapse.status',
                                             'LocalRelapse.time', 'LocalRelapse.status',
                                             'OS.time', 'OS.status',
                                             'DRFS.time', 'DRFS.status',
                                             'DFS.time', 'DFS.status')]

        if (erstatus=='ern') {
            test[test$GRADE!='G3','GRADE'] <- 'G1.G2'
            test[!test$ONCOTREE_CODE!='IDC','ONCOTREE_CODE'] <- 'BREAST'
        }

        test <- test %>% mutate_if(is.character, factor)
    }

    graph.steps <- loadGraph(paste0('out/mfcimgmStEPS.metabric.rna.', erstatus ,'.rank.', cvString, '.txt'))
    graph.stars <- loadGraph(paste0('out/mfcimgmStARS.metabric.rna.', erstatus ,'.rank.', cvString, '.txt'))
    graph.ebic <- loadGraph(paste0('out/mfcimgmEBIC.metabric.rna.', erstatus ,'.rank.', cvString, '.txt'))
    graph.bic <- loadGraph(paste0('out/mfcimgmBIC.metabric.rna.', erstatus ,'.rank.', cvString, '.txt'))

    return(list(train=train, test=test,
                graph=list(StEPS=graph.steps,
                           StARS=graph.stars,
                           EBIC=graph.ebic,
                           BIC=graph.bic)))
}


loadExtData <- function(ext.ds, trainList, feats) {

    print(ext.ds)
    
    extData <- read.csv(paste0('../external_validation/',
                               ifelse(ext.ds=='GSE96058','scan-b',ext.ds),
                               '/',ext.ds, '_combat.csv'), check.names=T, row.names=1)

    exp.surv.feats <- c('DSS.time', 'DSS.status',
                        'OS.time', 'OS.status',
                        'DRFS.time', 'DRFS.status',
                        'DFS.time', 'DFS.status')

    feats <- unique(c('ER_STATUS', feats))

    cleanFeats <- gsub('_grp', '', feats)

    extData <- na.omit(extData[,intersect(c(cleanFeats, exp.surv.feats), colnames(extData))])
    
    ## extData <- extData[rownames(na.omit(extData[,!colnames(extData) %in% exp.surv.feats])),]

    ## ext.genes <- sapply(strsplit(colnames(extData), ' /// '), function(x) x[1])

    ## genes <- gsub('_grp', '', colnames(exData))
    ## feats <- gsub('HLA[.]', 'HLA-', feats)

    if (any(c('DSS.time', 'DSS.status') %in% colnames(extData))) {
        extData$DSS <- Surv(extData$DSS.time, extData$DSS.status)
    }
    if (any(c('OS.time', 'OS.status') %in% colnames(extData))) {
        extData$OS <- Surv(extData$OS.time, extData$OS.status)
    }
    if (any(c('DRFS.time', 'DRFS.status') %in% colnames(extData))) {
        extData$DRFS <- Surv(extData$DRFS.time, extData$DRFS.status)
    }
    if (any(c('DFS.time', 'DFS.status') %in% colnames(extData))) {
        extData$DFS <- Surv(extData$DFS.time, extData$DFS.status)
    }

    extData <- extData[,!colnames(extData) %in% exp.surv.feats]

    numeric.mask <- sapply(extData, function(x) !is.Surv(x) & is.numeric(x))
    ## numeric.mask[c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')] <- FALSE

    extData$ER_STATUS <- factor(extData$ER_STATUS, levels=c('Negative', 'Positive'))
    extData$LYMPH_NODE_STATUS <- factor(extData$LYMPH_NODE_STATUS, levels=c('Negative', 'Positive'))

    extData.npn <- extData

    par(mfrow=c(5, 5))
    for (col in colnames(extData)[numeric.mask]) {
        if (col %in% feats) {
            feat <- col
        } else {
            feat <- paste0(col, '_grp')
        }
        
        npnfun <- approxfun(x=trainList$train[,feat], y=trainList$train.npn[,feat], method='constant',
                            yleft=min(trainList$train.npn[,feat]),
                            yright=max(trainList$train.npn[,feat]))
        plot(x=sort(trainList$train[,feat]), y=sort(trainList$train.npn[,feat]), type='l', main=col,
             xlab='Original', ylab='NPN')
        points(x=extData[,col], y=npnfun(extData[,col]), col='red')
        extData.npn[,col] <- npnfun(extData[,col])
    }

    ## extData.npn[,numeric.mask] <- huge::huge.npn(as.matrix(extData[,numeric.mask]))
    ## extData.npn[,'TUMOR_SIZE'] <- log10(extData[,'TUMOR_SIZE'])

    ## if (ext.ds %in% c('GSE11121', 'GSE19615', 'GSE42568', 'GSE6532', 'GSE7390', 'GSE9195')) {
    ##     extData.npn[,'TUMOR_SIZE']  <- 10 * extData.npn[,'TUMOR_SIZE']
    ## }


    ## extData.npn <- extData %>%
    ##     mutate_if(function(x) !is.Surv(x) & is.numeric(x),
    ##               function(x) as.vector(huge::huge.npn(as.matrix(x))))
    
    for (feat in feats[grepl('_grp', feats)]) {
        colnames(extData)[grepl(gsub('_grp', '', feat), colnames(extData))] <- feat
        colnames(extData.npn)[grepl(gsub('_grp', '', feat), colnames(extData.npn))] <- feat
    }

    surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')
    
    print("Missing Features")
    print(feats[!feats %in% union(colnames(extData.npn), surv.feats)])

    return(list(original=extData, npn=extData.npn))
}

getNeighbors <- function(graph, target) {
    splitedges <- strsplit(graph$edges[grep(target, graph$edges)], ' ')
    neighbors <- setdiff(unique(as.vector(unlist(sapply(splitedges,
                                                        function(x) {
                                                            if (target %in% x) {
                                                                return(x[c(1,3)])
                                                            }
                                                        })))), target)
    return(neighbors)
}


evaluateModels <- function(model.list, ext.data.list, times, outcome) {
    eval.df <- data.frame()
    for (ds.name in names(ext.data.list)) {
        test <- ext.data.list[[ds.name]]
        if (!'AGE_AT_DIAGNOSIS' %in% colnames(test)) {
            test$AGE_AT_DIAGNOSIS <- 0
        }

        if (ds.name=='GSE11121') {
            test[,'DRFS'][,1] <- 30.437 * test[,'DRFS'][,1]
        }

        test$LYMPH_NODE_STATUS <- factor(test$LYMPH_NODE_STATUS,
                                         levels=c('Negative','Positive'))

        test$ER_STATUS <- factor(test$ER_STATUS,
                                 levels=c('Negative','Positive'))

        erPos <- any(test$ER_STATUS=='Positive')
        erNeg <- any(test$ER_STATUS=='Negative')

        if (outcome %in% colnames(test)) {
            print(paste(ds.name, 'Max time:', max(test[,outcome][,1], na.rm=T)/365.25, 'Years'))
            for (t in times) {
                conc.uno <- concordance(model.list[['all']],
                                        newdata=test[!is.na(test[,outcome]),],
                                        ymax=365.25*t,
                                        timewt='n/G2')
                conc.harrell <- concordance(model.list[['all']],
                                            newdata=test[!is.na(test[,outcome]),],
                                            ymax=365.25*t)
                eval.df <- rbind(eval.df,
                                 data.frame(`Uno Concordance` = conc.uno$concordance,
                                            `Harrell Concordance` = conc.harrell$concordance,
                                            `Uno SE` = sqrt(conc.uno$var),
                                            `Harrell SE` = sqrt(conc.harrell$var),
                                            Model = 'Stratified Causal',
                                            `ER Status` = 'All',
                                            Time = paste(t, 'Years'),
                                            Outcome = outcome,
                                            Dataset = ds.name,
                                            check.names=F))

                if (erPos) {
                    conc.uno <- concordance(model.list[['all']],
                                            newdata=test[!is.na(test[,outcome]) &
                                                         test$ER_STATUS=='Positive',],
                                            ymax=365.25*t,
                                            timewt='n/G2')
                    conc.harrell <- concordance(model.list[['all']],
                                                newdata=test[!is.na(test[,outcome]) &
                                                             test$ER_STATUS=='Positive',],
                                                ymax=365.25*t)
                    eval.df <- rbind(eval.df,
                                     data.frame(`Uno Concordance` = conc.uno$concordance,
                                                `Harrell Concordance` = conc.harrell$concordance,
                                                `Uno SE` = sqrt(conc.uno$var),
                                                `Harrell SE` = sqrt(conc.harrell$var),
                                                Model = 'Stratified Causal',
                                                `ER Status` = 'ER+',
                                                Time = paste(t, 'Years'),
                                                Outcome = outcome,
                                                Dataset = ds.name,
                                                check.names=F))
                }

                if (erNeg) {
                    conc.uno <- concordance(model.list[['all']],
                                            newdata=test[!is.na(test[,outcome]) &
                                                         test$ER_STATUS=='Negative',],
                                            ymax=365.25*t,
                                            timewt='n/G2')
                    conc.harrell <- concordance(model.list[['all']],
                                                newdata=test[!is.na(test[,outcome]) &
                                                             test$ER_STATUS=='Negative',],
                                                ymax=365.25*t)
                    eval.df <- rbind(eval.df,
                                     data.frame(`Uno Concordance` = conc.uno$concordance,
                                                `Harrell Concordance` = conc.harrell$concordance,
                                                `Uno SE` = sqrt(conc.uno$var),
                                                `Harrell SE` = sqrt(conc.harrell$var),
                                                Model = 'Stratified Causal',
                                                `ER Status` = 'ER-',
                                                Time = paste(t, 'Years'),
                                                Outcome = outcome,
                                                Dataset = ds.name,
                                                check.names=F))
                }

                risk <- rep(0, nrow(test))
                if (sum(test$ER_STATUS=='Negative') > 0) {
                    risk[test$ER_STATUS=='Negative'] <- -predict(model.list[['ern']],
                                                                 newdata=test[test$ER_STATUS=='Negative',])
                }
                if (sum(test$ER_STATUS=='Positive') > 0) {
                    risk[test$ER_STATUS=='Positive'] <- -predict(model.list[['erp']],
                                                                 newdata=test[test$ER_STATUS=='Positive',])
                }
                
                conc.uno <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t, timewt='n/G2')
                conc.harrell <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t)
                
                eval.df <- rbind(eval.df,
                                 data.frame(
                                     `Uno Concordance` = conc.uno$concordance,
                                     `Harrell Concordance` = conc.harrell$concordance,
                                     `Uno SE` = sqrt(conc.uno$var),
                                     `Harrell SE` = sqrt(conc.harrell$var),
                                     Model = 'Split Causal',
                                     `ER Status` = 'All',
                                     Time = paste(t, 'Years'),
                                     Outcome = outcome,
                                     Dataset = ds.name,
                                     check.names=F))
                if(erPos) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Split Causal',
                                         `ER Status` = 'ER+',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }

                if (erNeg) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Split Causal',
                                         `ER Status` = 'ER-',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }
                
            }
        }
    }

    return(eval.df)
}

evaluateModels <- function(model.list, ext.data.list, times, outcome) {
    eval.df <- data.frame()
    for (ds.name in names(ext.data.list)) {
        test <- ext.data.list[[ds.name]]
        if (!'AGE_AT_DIAGNOSIS' %in% colnames(test)) {
            test$AGE_AT_DIAGNOSIS <- 0
        }

        if (ds.name=='GSE11121') {
            test[,'DRFS'][,1] <- 30.437 * test[,'DRFS'][,1]
        }

        if (ds.name=='GSE96058') {
            test[,'HLA.DRB4'] <- 0
        }

        test$LYMPH_NODE_STATUS <- factor(test$LYMPH_NODE_STATUS,
                                         levels=c('Negative','Positive'))

        test$ER_STATUS <- factor(test$ER_STATUS,
                                 levels=c('Negative','Positive'))

        erPos <- any(test$ER_STATUS=='Positive')
        erNeg <- any(test$ER_STATUS=='Negative')
        
        if (outcome %in% colnames(test)) {
            print(paste(ds.name, 'Max time:', max(test[,outcome][,1], na.rm=T)/365.25, 'Years'))
            for (t in times) {
                conc.uno <- concordance(model.list[['nb.all']],
                                        newdata=test[!is.na(test[,outcome]),],
                                        ymax=365.25*t,
                                        timewt='n/G2')
                conc.harrell <- concordance(model.list[['nb.all']],
                                            newdata=test[!is.na(test[,outcome]),],
                                            ymax=365.25*t)
                eval.df <- rbind(eval.df,
                                 data.frame(`Uno Concordance` = conc.uno$concordance,
                                            `Harrell Concordance` = conc.harrell$concordance,
                                            `Uno SE` = sqrt(conc.uno$var),
                                            `Harrell SE` = sqrt(conc.harrell$var),
                                            Model = 'Stratified Causal NB',
                                            `ER Status` = 'All',
                                            Time = paste(t, 'Years'),
                                            Outcome = outcome,
                                            Dataset = ds.name,
                                            check.names=F))

                if (erPos) {
                    conc.uno <- concordance(model.list[['nb.all']],
                                            newdata=test[!is.na(test[,outcome]) &
                                                         test$ER_STATUS=='Positive',],
                                            ymax=365.25*t,
                                            timewt='n/G2')
                    conc.harrell <- concordance(model.list[['nb.all']],
                                                newdata=test[!is.na(test[,outcome]) &
                                                             test$ER_STATUS=='Positive',],
                                                ymax=365.25*t)
                    eval.df <- rbind(eval.df,
                                     data.frame(`Uno Concordance` = conc.uno$concordance,
                                                `Harrell Concordance` = conc.harrell$concordance,
                                                `Uno SE` = sqrt(conc.uno$var),
                                                `Harrell SE` = sqrt(conc.harrell$var),
                                                Model = 'Stratified Causal NB',
                                                `ER Status` = 'ER+',
                                                Time = paste(t, 'Years'),
                                                Outcome = outcome,
                                                Dataset = ds.name,
                                                check.names=F))
                }

                if (erNeg) {
                    conc.uno <- concordance(model.list[['nb.all']],
                                            newdata=test[!is.na(test[,outcome]) &
                                                         test$ER_STATUS=='Negative',],
                                            ymax=365.25*t,
                                            timewt='n/G2')
                    conc.harrell <- concordance(model.list[['nb.all']],
                                                newdata=test[!is.na(test[,outcome]) &
                                                             test$ER_STATUS=='Negative',],
                                                ymax=365.25*t)
                    eval.df <- rbind(eval.df,
                                     data.frame(`Uno Concordance` = conc.uno$concordance,
                                                `Harrell Concordance` = conc.harrell$concordance,
                                                `Uno SE` = sqrt(conc.uno$var),
                                                `Harrell SE` = sqrt(conc.harrell$var),
                                                Model = 'Stratified Causal NB',
                                                `ER Status` = 'ER-',
                                                Time = paste(t, 'Years'),
                                                Outcome = outcome,
                                                Dataset = ds.name,
                                                check.names=F))
                }
                
                conc.uno <- concordance(model.list[['mb.bic.all']],
                                        newdata=test[!is.na(test[,outcome]),],
                                        ymax=365.25*t,
                                        timewt='n/G2')
                conc.harrell <- concordance(model.list[['mb.bic.all']],
                                            newdata=test[!is.na(test[,outcome]),],
                                            ymax=365.25*t)
                eval.df <- rbind(eval.df,
                                 data.frame(`Uno Concordance` = conc.uno$concordance,
                                            `Harrell Concordance` = conc.harrell$concordance,
                                            `Uno SE` = sqrt(conc.uno$var),
                                            `Harrell SE` = sqrt(conc.harrell$var),
                                            Model = 'Stratified Causal BIC MB',
                                            `ER Status` = 'All',
                                            Time = paste(t, 'Years'),
                                            Outcome = outcome,
                                            Dataset = ds.name,
                                            check.names=F))

                if(erPos) {
                    conc.uno <- concordance(model.list[['mb.bic.all']],
                                            newdata=test[!is.na(test[,outcome]) &
                                                         test$ER_STATUS=='Positive',],
                                            ymax=365.25*t,
                                            timewt='n/G2')
                    conc.harrell <- concordance(model.list[['mb.bic.all']],
                                                newdata=test[!is.na(test[,outcome]) &
                                                             test$ER_STATUS=='Positive',],
                                                ymax=365.25*t)
                    eval.df <- rbind(eval.df,
                                     data.frame(`Uno Concordance` = conc.uno$concordance,
                                                `Harrell Concordance` = conc.harrell$concordance,
                                                `Uno SE` = sqrt(conc.uno$var),
                                                `Harrell SE` = sqrt(conc.harrell$var),
                                                Model = 'Stratified Causal BIC MB',
                                                `ER Status` = 'ER+',
                                                Time = paste(t, 'Years'),
                                                Outcome = outcome,
                                                Dataset = ds.name,
                                                check.names=F))
                }
                
                if(erNeg) {
                    conc.uno <- concordance(model.list[['mb.bic.all']],
                                            newdata=test[!is.na(test[,outcome]) &
                                                         test$ER_STATUS=='Negative',],
                                            ymax=365.25*t,
                                            timewt='n/G2')
                    conc.harrell <- concordance(model.list[['mb.bic.all']],
                                                newdata=test[!is.na(test[,outcome]) &
                                                             test$ER_STATUS=='Negative',],
                                                ymax=365.25*t)
                    eval.df <- rbind(eval.df,
                                     data.frame(`Uno Concordance` = conc.uno$concordance,
                                                `Harrell Concordance` = conc.harrell$concordance,
                                                `Uno SE` = sqrt(conc.uno$var),
                                                `Harrell SE` = sqrt(conc.harrell$var),
                                                Model = 'Stratified Causal BIC MB',
                                                `ER Status` = 'ER-',
                                                Time = paste(t, 'Years'),
                                                Outcome = outcome,
                                                Dataset = ds.name,
                                                check.names=F))
                }

                feats <- dimnames(coef(model.list[['lasso.all']]))[[1]]
                mod.mat <- data.frame(model.matrix(~., test[,!colnames(test) %in% c('ER_STATUS', 'DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS') & colSums(is.na(test))==0]))
                mod.mat[,feats[!feats %in% colnames(mod.mat)]] <- 0
                
                risk <- -predict(model.list[['lasso.all']], newx=as.matrix(mod.mat[,feats]), s=model.list[['lasso.all']]$lambda.1se)
                
                conc.uno <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t, timewt='n/G2')
                conc.harrell <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t)
                
                eval.df <- rbind(eval.df,
                                 data.frame(
                                     `Uno Concordance` = conc.uno$concordance,
                                     `Harrell Concordance` = conc.harrell$concordance,
                                     `Uno SE` = sqrt(conc.uno$var),
                                     `Harrell SE` = sqrt(conc.harrell$var),
                                     Model = 'Stratified LASSO 1SE',
                                     `ER Status` = 'All',
                                     Time = paste(t, 'Years'),
                                     Outcome = outcome,
                                     Dataset = ds.name,
                                     check.names=F))

                if(erPos) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Stratified LASSO 1SE',
                                         `ER Status` = 'ER+',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }

                if(erNeg) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Stratified LASSO 1SE',
                                         `ER Status` = 'ER-',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }
                
                ## Split Models
                
                risk <- rep(0, nrow(test))
                if (sum(test$ER_STATUS=='Negative') > 0) {
                    risk[test$ER_STATUS=='Negative'] <- -predict(model.list[['ern']],
                                                                 newdata=test[test$ER_STATUS=='Negative',])
                }
                if (sum(test$ER_STATUS=='Positive') > 0) {
                    risk[test$ER_STATUS=='Positive'] <- -predict(model.list[['erp']],
                                                                 newdata=test[test$ER_STATUS=='Positive',])
                }
                
                conc.uno <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t, timewt='n/G2')
                conc.harrell <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t)
                
                eval.df <- rbind(eval.df,
                                 data.frame(
                                     `Uno Concordance` = conc.uno$concordance,
                                     `Harrell Concordance` = conc.harrell$concordance,
                                     `Uno SE` = sqrt(conc.uno$var),
                                     `Harrell SE` = sqrt(conc.harrell$var),
                                     Model = 'Split Causal NB',
                                     `ER Status` = 'All',
                                     Time = paste(t, 'Years'),
                                     Outcome = outcome,
                                     Dataset = ds.name,
                                     check.names=F))
                if(erPos) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Split Causal NB',
                                         `ER Status` = 'ER+',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }

                if (erNeg) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Split Causal NB',
                                         `ER Status` = 'ER-',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }

                risk <- rep(0, nrow(test))
                if (sum(test$ER_STATUS=='Negative') > 0) {
                    risk[test$ER_STATUS=='Negative'] <- -predict(model.list[['mb.bic.ern']],
                                                                 newdata=test[test$ER_STATUS=='Negative',])
                }
                if (sum(test$ER_STATUS=='Positive') > 0) {
                    risk[test$ER_STATUS=='Positive'] <- -predict(model.list[['mb.bic.erp']],
                                                                 newdata=test[test$ER_STATUS=='Positive',])
                }
                
                conc.uno <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t, timewt='n/G2')
                conc.harrell <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t)
                
                eval.df <- rbind(eval.df,
                                 data.frame(
                                     `Uno Concordance` = conc.uno$concordance,
                                     `Harrell Concordance` = conc.harrell$concordance,
                                     `Uno SE` = sqrt(conc.uno$var),
                                     `Harrell SE` = sqrt(conc.harrell$var),
                                     Model = 'Split Causal BIC MB',
                                     `ER Status` = 'All',
                                     Time = paste(t, 'Years'),
                                     Outcome = outcome,
                                     Dataset = ds.name,
                                     check.names=F))

                if(erPos) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Split Causal BIC MB',
                                         `ER Status` = 'ER+',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }

                if(erNeg) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Split Causal BIC MB',
                                         `ER Status` = 'ER-',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }

                risk <- rep(0, nrow(test))
                if (sum(test$ER_STATUS=='Negative') > 0) {
                    risk[test$ER_STATUS=='Negative'] <- -predict(model.list[['lasso.ern']],
                                                                 newx=as.matrix(mod.mat[test$ER_STATUS=='Negative',feats]),
                                                                 s=model.list[['lasso.ern']]$lambda.1se)
                }
                if (sum(test$ER_STATUS=='Positive') > 0) {
                    risk[test$ER_STATUS=='Positive'] <- -predict(model.list[['lasso.erp']],
                                                                 newx=as.matrix(mod.mat[test$ER_STATUS=='Positive',feats]),
                                                                 s=model.list[['lasso.erp']]$lambda.1se)
                }
                
                conc.uno <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t, timewt='n/G2')
                conc.harrell <- concordance(
                    test[,outcome] ~ strata(test$ER_STATUS) + risk,
                    ymax=365.25*t)
                
                eval.df <- rbind(eval.df,
                                 data.frame(
                                     `Uno Concordance` = conc.uno$concordance,
                                     `Harrell Concordance` = conc.harrell$concordance,
                                     `Uno SE` = sqrt(conc.uno$var),
                                     `Harrell SE` = sqrt(conc.harrell$var),
                                     Model = 'Split LASSO 1SE',
                                     `ER Status` = 'All',
                                     Time = paste(t, 'Years'),
                                     Outcome = outcome,
                                     Dataset = ds.name,
                                     check.names=F))

                if(erPos) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Positive',outcome] ~ risk[test$ER_STATUS=='Positive'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Split LASSO 1SE',
                                         `ER Status` = 'ER+',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }

                if(erNeg) {
                    conc.uno <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t, timewt='n/G2')
                    conc.harrell <- concordance(
                        test[test$ER_STATUS=='Negative',outcome] ~ risk[test$ER_STATUS=='Negative'],
                        ymax=365.25*t)
                    
                    eval.df <- rbind(eval.df,
                                     data.frame(
                                         `Uno Concordance` = conc.uno$concordance,
                                         `Harrell Concordance` = conc.harrell$concordance,
                                         `Uno SE` = sqrt(conc.uno$var),
                                         `Harrell SE` = sqrt(conc.harrell$var),
                                         Model = 'Split LASSO 1SE',
                                         `ER Status` = 'ER-',
                                         Time = paste(t, 'Years'),
                                         Outcome = outcome,
                                         Dataset = ds.name,
                                         check.names=F))
                }
                
            }
        }
    }

    return(eval.df)
}


surv.feats <- c('DSS', 'OD', 'DistantRelapse', 'LocalRelapse', 'OS', 'DRFS', 'DFS')

ext.ds.names <- c('GSE11121', 'GSE19615', 'GSE3494', 'GSE42568', 'GSE45255', 'GSE6532',
                  'GSE7390', 'GSE9195', 'GSE96058')

ext.ds.names <- c('GSE19615', 'GSE3494', 'GSE42568', 'GSE45255', 'GSE6532',
                  'GSE7390', 'GSE9195', 'GSE96058')


fullData.list <- lapply(list(all='all', erp='erp', ern='ern'), loadFullData)
fullGraph.list <- lapply(list(all='all', erp='erp', ern='ern'), loadFullGraph)

fullGraph.list$ern$lambda

ern.feats <- unique(unlist(fullGraph.list$ern$neighbors[c('DSS', 'DistantRelapse', 'LocalRelapse', 'OD')]))
erp.feats <- unique(unlist(fullGraph.list$erp$neighbors[c('DSS', 'DistantRelapse', 'LocalRelapse', 'OD')]))
feats <- union(erp.feats, ern.feats)

feats <- setdiff(colnames(fullData.list$all$train.npn), c('HER2_STATUS', 'PR_STATUS', surv.feats))

extData.list <- lapply(ext.ds.names, loadExtData, fullData.list$all, feats)
names(extData.list) <- ext.ds.names

par(mfrow=c(3,3))
for (ds.name in ext.ds.names) {
    hist(extData.list[[ds.name]]$npn$TUMOR_SIZE, main=ds.name)
}

par(mfrow=c(3,3))
for (ds.name in ext.ds.names) {
    hist(extData.list[[ds.name]]$npn$AGE_AT_DIAGNOSIS, main=ds.name)
}

par(mfrow=c(3,3))
for (ds.name in ext.ds.names) {
    plot(extData.list[[ds.name]]$npn$ER_STATUS, main=ds.name)
}


shared.feats <- setdiff(colnames(extData.list[[1]]$npn), surv.feats)


pdf('rf_model_train_plots.pdf', width=8, height=8)
rfModels <- list()
for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {

    temp <- fullData.list$all$train.npn[,shared.feats]
    temp$time <- fullData.list$all$train.npn[,sf][,1]
    temp$status <- fullData.list$all$train.npn[,sf][,2]

    rfModels[[sf]] <- rfsrc(Surv(time, status) ~ ., temp,
                            ntree=5000, importance="anti")
    
    print(rfModels[[sf]])
    print(plot(rfModels[[sf]]))

}
dev.off()


pdf('lasso_model_train_plots.pdf', width=8, height=8)
lassoModels <- list()
for (erstatus in c('nostratify', 'all', 'ern', 'erp')) {
    lassoModels[[erstatus]] <- list()

    temp <- fullData.list$all$train.npn[,shared.feats]

    if (erstatus=='ern') {
        temp <- temp %>% filter(ER_STATUS=='Negative') %>% dplyr::select(-c('ER_STATUS'))
    }
    if (erstatus=='erp') {
        temp <- temp %>% filter(ER_STATUS=='Positive') %>% dplyr::select(-c('ER_STATUS'))
    }
    if (erstatus=='all') {
        temp <- temp %>% dplyr::select(-c('ER_STATUS'))
    }

    for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
        if (erstatus!='nostratify') {
            lassoModels[[erstatus]][[sf]] <- cv.glmnet(x=model.matrix(~., temp)[,-1],
                                                       y=fullData.list[[erstatus]]$train.npn[,sf],
                                                       family='cox', lambda.min.ratio=0.01)
        } else {
            lassoModels[[erstatus]][[sf]] <- cv.glmnet(x=model.matrix(~., temp)[,-1],
                                                       y=Surv(fullData.list$all$train.npn[,sf][,1],
                                                              fullData.list$all$train.npn[,sf][,2]),
                                                       family='cox', lambda.min.ratio=0.01)
        }
        print(lassoModels[[erstatus]][[sf]])
        print(plot(lassoModels[[erstatus]][[sf]], main=paste('LASSO',sf)))
    }
}
dev.off()

saveRDS(lassoModels, 'coxlasso_final_models.rds')
saveRDS(rfModels, 'rsf_final_models.rds')

lassoModels <- readRDS('coxlasso_final_models.rds')
rfModels <- readRDS('rsf_final_models.rds')

library(mstate)

prepTrainMSM <- function(data, graphs, tmat) {
    ## data <- rbind.data.frame(train, test)

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

    ## mstrain <- msdata[msdata$id %in% rownames(train),]
    ## mstest <- msdata[msdata$id %in% rownames(test),]

    return(msdata)
}


prepTestMSM <- function(data, graphs, tmat) {
    ## data <- rbind.data.frame(train, test)

    msm.surv.feats <- colnames(tmat)[-1]
    for (sf.idx in 1:4) {
        sf <- msm.surv.feats[sf.idx]
        data[,paste(sf, 'status', sep='.')] <- 0
        ## for (sf.jdx in (sf.idx+1):4) {
        ##     sf2 <- msm.surv.feats[sf.jdx]
        ##     if (is.na(sf2)) break
        ##     equal.vec <- apply(data[,sf][,1:2]==data[,sf2][,1:2], 1, all)
        ##     data[equal.vec & data[,sf][,2],sf2][,1] <- data[equal.vec & data[,sf][,2],sf2][,1] + 1
        ## }
        data[,paste(sf, 'time', sep='.')] <- runif(nrow(data)) + sf.idx
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

    ## mstrain <- msdata[msdata$id %in% rownames(train),]
    ## mstest <- msdata[msdata$id %in% rownames(test),]

    return(msdata)
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
    pt <- probtrans(msf1, predt=0, direction='forward')

    nTime <- nrow(pt[[1]])

    chf.df <- data.frame(Time=pt[[1]][1:nTime,'time'],
                         DSS=-log(1-pt[[1]][1:nTime,'pstate4']),
                         OS=-log(1-pt[[1]][1:nTime,'pstate4']-pt[[1]][1:nTime,'pstate5']),
                         DRFS=-log(1-pt[[1]][1:nTime,'pstate3']-pt[[1]][1:nTime,'pstate4']-pt[[1]][1:nTime,'pstate5']),
                         DFS=-log(pt[[1]][1:nTime,'pstate1']))

    chf.df <- na.omit(chf.df)

    chf.fun <- list()
    for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
        chf.fun[[sf]] <- approxfun(chf.df[,'Time'], chf.df[,sf],
                                   yleft=0, yright=max(chf.df[,sf]), method='constant')
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
                                                     ER_STATUS=mstest[mstest$id==samp,'ER_STATUS'][1],
                                                     tail(na.omit(chf.out$chf.df),1))
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
                      }, mc.cores=2))

    return(riskdf)
}

msmPredictMortality <- function(model, train, mstest, tmat) {
    require(parallel)
    
    samps <- unique(mstest$id)

    mortdf <- do.call(rbind.data.frame,
                      mclapply(samps, function(samp) {
                          chf.out <- patient.chf(model, mstest[mstest$id==samp,], tmat)

                          single <- data.frame(id=samp,
                                               DSS=0,
                                               OS=0,
                                               DRFS=0,
                                               DFS=0)
                          
                          for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
                              eventTimes <- unique(train[,sf][train[,sf][,2]==1,1])
                              for (tt in eventTimes) {
                                  single[,sf] <- single[,sf] + chf.out$chf.fun[[sf]](tt)
                              }
                          }
                          
                          return(single)
                      }, mc.cores=2))

    return(mortdf)
}



tmat <- transMat(list(c(2,3,4,5), c(3,4,5), c(4, 5), c(), c()),
                 c('Surgery', 'LocalRelapse', 'DistantRelapse', 'DSS', 'OD'))

fullNPN <- fullData.list$all$train.npn

## for (sf in surv.feats) {
##     fullNPN[,sf][fullNPN[,sf][,1] > 25 * 365.25,1] <- 25 * 365.25
##     fullNPN[,sf][fullNPN[,sf][,1] == 25 * 365.25,2] <- 0
## }

mstrain <- prepTrainMSM(fullNPN, fullGraph.list, tmat)

model <- msmTrain(fullNPN, mstrain, fullGraph.list, tmat)

msmconcdf <- data.frame()
for (ext.ds in names(extData.list)) {

    print(ext.ds)
    print(dim(extData.list[[ext.ds]]$npn))

    mstest <- prepTestMSM(extData.list[[ext.ds]]$npn, fullGraph.list, tmat)

    riskdf <- msmPredict(model, mstest, tmat, c(3,5,8))

    for (tt in c(3,5,8)) {
        temp <- riskdf %>% filter(Time==paste(tt,'Years'))

        ids <- as.character(temp$id)

        rownames(temp) <- ids

        temp <- cbind(extData.list[[ext.ds]]$npn[ids,], risk=-log(temp[ids,c('DSS','OS','DRFS','DFS')]))

        for (sf in intersect(colnames(extData.list[[ext.ds]]$npn), c('DSS','OS','DRFS','DFS'))) {

            f <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + risk.', sf))
            
            conc.harrell <- concordance(f, temp, ymax=365.25*tt)

            single <- data.frame(Dataset=ext.ds,
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='Stratified',
                                 `Concordance` = conc.harrell$concordance,
                                 `Concordance SE` = sqrt(conc.harrell$var),
                                 check.names=F)

            msmconcdf <- rbind(msmconcdf, single)

            f <- as.formula(paste0(sf, ' ~ risk.', sf))

            conc.harrell <- concordance(f, temp, ymax=365.25*tt)

            single <- data.frame(Dataset=ext.ds,
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='All',
                                 `Concordance` = conc.harrell$concordance,
                                 `Concordance SE` = sqrt(conc.harrell$var),
                                 check.names=F)

            msmconcdf <- rbind(msmconcdf, single)

            if (any(temp$ER_STATUS=='Positive')) {
                conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Positive'), ymax=365.25*tt)

                single <- data.frame(Dataset=ext.ds,
                                     Outcome=sf,
                                     Time=paste(tt,'Years'),
                                     `ER Status`='Positive',
                                     `Concordance` = conc.harrell$concordance,
                                     `Concordance SE` = sqrt(conc.harrell$var),
                                     check.names=F)

                msmconcdf <- rbind(msmconcdf, single)
            }

            if (any(temp$ER_STATUS=='Negative')) {
                conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Negative'), ymax=365.25*tt)

                single <- data.frame(Dataset=ext.ds,
                                     Outcome=sf,
                                     Time=paste(tt,'Years'),
                                     `ER Status`='Negative',
                                     `Concordance` = conc.harrell$concordance,
                                     `Concordance SE` = sqrt(conc.harrell$var),
                                     check.names=F)
                
                msmconcdf <- rbind(msmconcdf, single)
            }
        }
    }
}

## write.csv(msmconcdf, 'external_msm_chf_conc.csv')
## write.csv(msmconcdf, 'external_msm_mort_conc.csv')

msmconcdf <- read.csv('external_msm_chf_conc.csv', row.names=1, check.names=F)

head(msmconcdf)

ggplot(msmconcdf, ## %>% filter(`ER Status`!='All'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset)) +
    geom_pointrange(position=position_dodge(width=0.8)) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    theme_bw()


rfconcdf <- data.frame()
for (ext.ds in names(extData.list)) {

    print(ext.ds)
    print(dim(extData.list[[ext.ds]]$npn))

    temp <- extData.list[[ext.ds]]$npn

    for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
        pred <- predict(rfModels[[sf]], newdata=extData.list[[ext.ds]]$npn[,shared.feats])
        temp[,paste0('risk.',sf)] <- -pred$predicted
    }

    for (tt in c(3,5,8)) {

        for (sf in intersect(colnames(extData.list[[ext.ds]]$npn), c('DSS','OS','DRFS','DFS'))) {
    
            f <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + risk.', sf))
            
            conc.harrell <- concordance(f, temp, ymax=365.25*tt)

            single <- data.frame(Dataset=ext.ds,
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='Stratified',
                                 `Concordance` = conc.harrell$concordance,
                                 `Concordance SE` = sqrt(conc.harrell$var),
                                 check.names=F)

            rfconcdf <- rbind(rfconcdf, single)

            f <- as.formula(paste0(sf, ' ~ risk.', sf))

            conc.harrell <- concordance(f, temp, ymax=365.25*tt)

            single <- data.frame(Dataset=ext.ds,
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`='All',
                                 `Concordance` = conc.harrell$concordance,
                                 `Concordance SE` = sqrt(conc.harrell$var),
                                 check.names=F)

            rfconcdf <- rbind(rfconcdf, single)

            if (any(temp$ER_STATUS=='Positive')) {
                conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Positive'), ymax=365.25*tt)

                single <- data.frame(Dataset=ext.ds,
                                     Outcome=sf,
                                     Time=paste(tt,'Years'),
                                     `ER Status`='Positive',
                                     `Concordance` = conc.harrell$concordance,
                                     `Concordance SE` = sqrt(conc.harrell$var),
                                     check.names=F)

                rfconcdf <- rbind(rfconcdf, single)
            }

            if (any(temp$ER_STATUS=='Negative')) {
                conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Negative'), ymax=365.25*tt)

                single <- data.frame(Dataset=ext.ds,
                                     Outcome=sf,
                                     Time=paste(tt,'Years'),
                                     `ER Status`='Negative',
                                     `Concordance` = conc.harrell$concordance,
                                     `Concordance SE` = sqrt(conc.harrell$var),
                                     check.names=F)
                
                rfconcdf <- rbind(rfconcdf, single)
            }
        }
    }
}

head(rfconcdf)


## lassoPredict <- function(cv.out, target, shared.feats, train, test, times) {
##     require(pec)
##     eventTimes <- sort(unique(train[train[,target][,2]==1,target][,1]))

##     risk <- predict(
##         cv.out,
##         newx=as.matrix(model.matrix(~.,train[,setdiff(shared.feats, 'ER_STATUS')])[,-1]),
##         s=cv.out$lambda.min, type='response')

##     temp <- data.frame(train[,c(target, 'ER_STATUS')], risk=risk[,1])

##     f <- as.formula(paste(target, '~ strata(ER_STATUS) + risk'))
    
## }

lassoconcdf <- data.frame()
for (ext.ds in names(extData.list)) {

    print(ext.ds)
    print(dim(extData.list[[ext.ds]]$npn))

    for (erstatus in c('nostratify', 'all', 'erp', 'ern')) {

        temp <- extData.list[[ext.ds]]$npn

        for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
            if (erstatus != 'nostratify') {
                temp[,paste0('risk.',sf)] <- -as.vector(predict(lassoModels[[erstatus]][[sf]], newx=as.matrix(model.matrix(~.,extData.list[[ext.ds]]$npn[,setdiff(shared.feats, 'ER_STATUS')])[,-1]), s=lassoModels[[erstatus]][[sf]]$lambda.min, type='response'))
            } else {
                temp[,paste0('risk.',sf)] <- -as.vector(predict(lassoModels[[erstatus]][[sf]], newx=as.matrix(model.matrix(~.,extData.list[[ext.ds]]$npn[,shared.feats])[,-1]), s=lassoModels[[erstatus]][[sf]]$lambda.min, type='response'))
            }
        }

        if (erstatus=='ern') {
            temp <- temp %>% filter(ER_STATUS=='Negative') ## %>% dplyr::select(-c('ER_STATUS'))
        }
        if (erstatus=='erp') {
            temp <- temp %>% filter(ER_STATUS=='Positive') ## %>% dplyr::select(-c('ER_STATUS'))
        }
        ## if (erstatus=='all') {
        ##     temp <- temp ## %>% dplyr::select(-c('ER_STATUS'))
        ## }


        for (tt in c(3,5,8)) {

            for (sf in intersect(colnames(extData.list[[ext.ds]]$npn), c('DSS','OS','DRFS','DFS'))) {

                if (erstatus=='all') {
                    f <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + risk.', sf))
                    
                    conc.harrell <- concordance(f, temp, ymax=365.25*tt)

                    single <- data.frame(Dataset=ext.ds,
                                         Outcome=sf,
                                         Time=paste(tt,'Years'),
                                         `ER Status`='Stratified',
                                         `Concordance` = conc.harrell$concordance,
                                         `Concordance SE` = sqrt(conc.harrell$var),
                                         check.names=F)

                    lassoconcdf <- rbind(lassoconcdf, single)
                } else {
                    f <- as.formula(paste0(sf, ' ~ risk.', sf))

                    if (erstatus=='nostratify') {
                        conc.harrell <- concordance(f, temp, ymax=365.25*tt)
                    
                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='All',
                                             `Concordance` = conc.harrell$concordance,
                                             `Concordance SE` = sqrt(conc.harrell$var),
                                             check.names=F)

                        lassoconcdf <- rbind(lassoconcdf, single)
                    }
                    if (erstatus=='erp' && any(temp$ER_STATUS=='Positive')) {
                        conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Positive'), ymax=365.25*tt)

                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='Positive',
                                             `Concordance` = conc.harrell$concordance,
                                             `Concordance SE` = sqrt(conc.harrell$var),
                                             check.names=F)

                        lassoconcdf <- rbind(lassoconcdf, single)
                    }

                    if (erstatus=='ern' && any(temp$ER_STATUS=='Negative')) {
                        conc.harrell <- concordance(f, temp %>% filter(ER_STATUS=='Negative'), ymax=365.25*tt)

                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='Negative',
                                             `Concordance` = conc.harrell$concordance,
                                             `Concordance SE` = sqrt(conc.harrell$var),
                                             check.names=F)
                
                        lassoconcdf <- rbind(lassoconcdf, single)
                    }
                }
            }
        }
    }
}

head(lassoconcdf)


head(msmconcdf)

library(survcomp)
msmmetaconcdf <- data.frame()
for (sf in c('DSS','OS','DRFS','DFS')) {
    for (tt in c(5)) {
        for (erstatus in c('All', 'Stratified', 'Positive', 'Negative')) {
            temp <- msmconcdf %>% filter(Outcome==sf, Time==paste(tt,'Years'), `ER Status`==erstatus)

            meta.conc <- unlist(combine.est(temp$Concordance, temp$`Concordance SE`))

            single <- data.frame(Dataset='Meta Cohort',
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`=erstatus,
                                 `Concordance` = meta.conc['estimate'],
                                 `Concordance SE` = meta.conc['se'],
                                 check.names=F)

            msmmetaconcdf <- rbind(msmmetaconcdf, single)
                
        }
    }
}

library(survcomp)
rfmetaconcdf <- data.frame()
for (sf in c('DSS','OS','DRFS','DFS')) {
    for (tt in c(5)) {
        for (erstatus in c('All', 'Stratified', 'Positive', 'Negative')) {
            temp <- rfconcdf %>% filter(Outcome==sf, Time==paste(tt,'Years'), `ER Status`==erstatus)

            meta.conc <- unlist(combine.est(temp$Concordance, temp$`Concordance SE`))

            single <- data.frame(Dataset='Meta Cohort',
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`=erstatus,
                                 `Concordance` = meta.conc['estimate'],
                                 `Concordance SE` = meta.conc['se'],
                                 check.names=F)

            rfmetaconcdf <- rbind(rfmetaconcdf, single)
                
        }
    }
}

library(survcomp)
lassometaconcdf <- data.frame()
for (sf in c('DSS','OS','DRFS','DFS')) {
    for (tt in c(5)) {
        for (erstatus in c('All', 'Stratified', 'Positive', 'Negative')) {
            temp <- lassoconcdf %>% filter(Outcome==sf, Time==paste(tt,'Years'), `ER Status`==erstatus)

            meta.conc <- unlist(combine.est(temp$Concordance, temp$`Concordance SE`))

            single <- data.frame(Dataset='Meta Cohort',
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`=erstatus,
                                 `Concordance` = meta.conc['estimate'],
                                 `Concordance SE` = meta.conc['se'],
                                 check.names=F)

            lassometaconcdf <- rbind(lassometaconcdf, single)
                
        }
    }
}


library(survcomp)
lassometaconcdiffdf <- data.frame()
for (sf in c('DSS','OS','DRFS','DFS')) {
    for (tt in c(5)) {
        for (erstatus in c('All', 'Stratified', 'Positive', 'Negative')) {
            temp <- lassoconcdiffdf %>% filter(Outcome==sf, Time==paste(tt,'Years'), `ER Status`==erstatus)

            meta.conc <- unlist(combine.est(temp$`Delta Concordance`, temp$`Delta Concordance SE`))

            single <- data.frame(Dataset='Meta Cohort',
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`=erstatus,
                                 `Delta Concordance` = meta.conc['estimate'],
                                 `Delta Concordance SE` = meta.conc['se'],
                                 check.names=F)

            lassometaconcdiffdf <- rbind(lassometaconcdiffdf, single)
                
        }
    }
}

lassometaconcdiffdf$stat <- lassometaconcdiffdf[,5] / lassometaconcdiffdf[,6]
lassometaconcdiffdf$pval <- 2*pnorm(-abs(lassometaconcdiffdf$stat))


library(survcomp)
rfmetaconcdiffdf <- data.frame()
for (sf in c('DSS','OS','DRFS','DFS')) {
    for (tt in c(5)) {
        for (erstatus in c('All', 'Stratified', 'Positive', 'Negative')) {
            temp <- rfconcdiffdf %>% filter(Outcome==sf, Time==paste(tt,'Years'), `ER Status`==erstatus)

            meta.conc <- unlist(combine.est(temp$`Delta Concordance`, temp$`Delta Concordance SE`))

            single <- data.frame(Dataset='Meta Cohort',
                                 Outcome=sf,
                                 Time=paste(tt,'Years'),
                                 `ER Status`=erstatus,
                                 `Delta Concordance` = meta.conc['estimate'],
                                 `Delta Concordance SE` = meta.conc['se'],
                                 check.names=F)

            rfmetaconcdiffdf <- rbind(rfmetaconcdiffdf, single)
                
        }
    }
}

rfmetaconcdiffdf$stat <- rfmetaconcdiffdf[,5] / rfmetaconcdiffdf[,6]
rfmetaconcdiffdf$pval <- 2*pnorm(-abs(rfmetaconcdiffdf$stat))


lassozstat <- (msmtemp$Concordance - lassotemp$Concordance) / sqrt(msmtemp$`Concordance SE`^2 + lassotemp$`Concordance SE`^2)

2*pnorm(-abs(lassozstat))

rfzstat <- (msmtemp$Concordance - rftemp$Concordance) / sqrt(msmtemp$`Concordance SE`^2 + rftemp$`Concordance SE`^2)

2*pnorm(-abs(rfzstat))


set.seed(43)
lassoconcdiffdf <- data.frame()
B <- 500
for (ext.ds in ext.ds.names) {
    print(ext.ds)
    print(dim(extData.list[[ext.ds]]$npn))

    mstest <- prepTestMSM(extData.list[[ext.ds]]$npn, fullGraph.list, tmat)

    msmriskdf <- msmPredict(model, mstest, tmat, c(5))

    for (tt in c(5)) {

        msmtemp <- msmriskdf %>% filter(Time==paste(tt,'Years'))

        ids <- as.character(msmtemp$id)

        rownames(msmtemp) <- ids

        msmtemp <- cbind(extData.list[[ext.ds]]$npn[ids,], risk=-log(msmtemp[ids,c('DSS','OS','DRFS','DFS')]))


        for (erstatus in c('nostratify', 'all', 'erp', 'ern')) {

            lassotemp <- extData.list[[ext.ds]]$npn[ids,]

            for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
                if (erstatus != 'nostratify') {
                    lassotemp[,paste0('risk.',sf)] <- -as.vector(predict(lassoModels[[erstatus]][[sf]], newx=as.matrix(model.matrix(~.,extData.list[[ext.ds]]$npn[,setdiff(shared.feats, 'ER_STATUS')])[,-1]), s=lassoModels[[erstatus]][[sf]]$lambda.min, type='response'))
                } else {
                    lassotemp[,paste0('risk.',sf)] <- -as.vector(predict(lassoModels[[erstatus]][[sf]], newx=as.matrix(model.matrix(~.,extData.list[[ext.ds]]$npn[,shared.feats])[,-1]), s=lassoModels[[erstatus]][[sf]]$lambda.min, type='response'))
                }
            }

            if (erstatus=='ern') {
                lassotemp <- lassotemp %>% filter(ER_STATUS=='Negative') ## %>% dplyr::select(-c('ER_STATUS'))
            }
            if (erstatus=='erp') {
                lassotemp <- lassotemp %>% filter(ER_STATUS=='Positive') ## %>% dplyr::select(-c('ER_STATUS'))
            }

            for (sf in intersect(colnames(extData.list[[ext.ds]]$npn), c('DSS','OS','DRFS','DFS'))) {

                if (erstatus=='all') {
                    f <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + risk.', sf))
                    
                    deltavec <- rep(0, B)

                    for (b in 1:B) {
                        while (TRUE) {
                            ss <- sample(1:nrow(lassotemp), nrow(lassotemp), replace=T)
                            lasso.conc.harrell <- concordance(f, lassotemp[ss,], ymax=365.25*tt)

                            msm.conc.harrell <- concordance(f, msmtemp[ss,], ymax=365.25*tt)

                            deltavec[b] <- msm.conc.harrell$concordance - lasso.conc.harrell$concordance
                            if (!is.na(deltavec[b])) {
                                break
                            }
                        }
                    }

                    lasso.conc.harrell <- concordance(f, lassotemp, ymax=365.25*tt)

                    msm.conc.harrell <- concordance(f, msmtemp, ymax=365.25*tt)


                    single <- data.frame(Dataset=ext.ds,
                                         Outcome=sf,
                                         Time=paste(tt,'Years'),
                                         `ER Status`='Stratified',
                                         `Delta Concordance` = msm.conc.harrell$concordance - lasso.conc.harrell$concordance,
                                         `Delta Concordance SE` = sd(deltavec),
                                         check.names=F)

                    lassoconcdiffdf <- rbind(lassoconcdiffdf, single)
                    
                } else {
                    f <- as.formula(paste0(sf, ' ~ risk.', sf))

                    if (erstatus=='nostratify') {
                        deltavec <- rep(0, B)

                        for (b in 1:B) {
                            while (TRUE) {
                                ss <- sample(1:nrow(lassotemp), nrow(lassotemp), replace=T)
                                lasso.conc.harrell <- concordance(f, lassotemp[ss,], ymax=365.25*tt)

                                msm.conc.harrell <- concordance(f, msmtemp[ss,], ymax=365.25*tt)

                                deltavec[b] <- msm.conc.harrell$concordance - lasso.conc.harrell$concordance
                                if (!is.na(deltavec[b])) {
                                    break
                                }
                            }
                        }

                        lasso.conc.harrell <- concordance(f, lassotemp, ymax=365.25*tt)

                        msm.conc.harrell <- concordance(f, msmtemp, ymax=365.25*tt)


                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='All',
                                             `Delta Concordance` = msm.conc.harrell$concordance - lasso.conc.harrell$concordance,
                                             `Delta Concordance SE` = sd(deltavec),
                                             check.names=F)

                        lassoconcdiffdf <- rbind(lassoconcdiffdf, single)
                    }
                    if (erstatus=='erp' && any(lassotemp$ER_STATUS=='Positive')) {

                        deltavec <- rep(0, B)

                        for (b in 1:B) {
                            while (TRUE) {
                                ss <- sample(1:sum(lassotemp$ER_STATUS=='Positive'), sum(lassotemp$ER_STATUS=='Positive'), replace=T)
                                lasso.conc.harrell <- concordance(f, (lassotemp %>% filter(ER_STATUS=='Positive'))[ss,], ymax=365.25*tt)

                                msm.conc.harrell <- concordance(f, (msmtemp %>% filter(ER_STATUS=='Positive'))[ss,], ymax=365.25*tt)

                                deltavec[b] <- msm.conc.harrell$concordance - lasso.conc.harrell$concordance
                                if (!is.na(deltavec[b])) {
                                    break
                                }
                            }
                        }

                        lasso.conc.harrell <- concordance(f, (lassotemp %>% filter(ER_STATUS=='Positive')), ymax=365.25*tt)

                        msm.conc.harrell <- concordance(f, (msmtemp %>% filter(ER_STATUS=='Positive')), ymax=365.25*tt)


                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='Positive',
                                             `Delta Concordance` = msm.conc.harrell$concordance - lasso.conc.harrell$concordance,
                                             `Delta Concordance SE` = sd(deltavec),
                                             check.names=F)

                        lassoconcdiffdf <- rbind(lassoconcdiffdf, single)
                    }

                    if (erstatus=='ern' && any(lassotemp$ER_STATUS=='Negative')) {
                        deltavec <- rep(0, B)

                        for (b in 1:B) {
                            while (TRUE) {
                                ss <- sample(1:sum(lassotemp$ER_STATUS=='Negative'), sum(lassotemp$ER_STATUS=='Negative'), replace=T)
                                lasso.conc.harrell <- concordance(f, (lassotemp %>% filter(ER_STATUS=='Negative'))[ss,], ymax=365.25*tt)

                                msm.conc.harrell <- concordance(f, (msmtemp %>% filter(ER_STATUS=='Negative'))[ss,], ymax=365.25*tt)

                                deltavec[b] <- msm.conc.harrell$concordance - lasso.conc.harrell$concordance
                                if (!is.na(deltavec[b])) {
                                    break
                                }
                            }
                        }

                        lasso.conc.harrell <- concordance(f, (lassotemp %>% filter(ER_STATUS=='Negative')), ymax=365.25*tt)

                        msm.conc.harrell <- concordance(f, (msmtemp %>% filter(ER_STATUS=='Negative')), ymax=365.25*tt)


                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='Negative',
                                             `Delta Concordance` = msm.conc.harrell$concordance - lasso.conc.harrell$concordance,
                                             `Delta Concordance SE` = sd(deltavec),
                                             check.names=F)

                        lassoconcdiffdf <- rbind(lassoconcdiffdf, single)
                    }
                }
            }
        }
    }
}

lassoconcdiffdf

write.csv(lassoconcdiffdf, "msm_lasso_delta_conc.csv")


set.seed(43)
rfconcdiffdf <- data.frame()
B <- 500
for (ext.ds in ext.ds.names) {
    print(ext.ds)
    print(dim(extData.list[[ext.ds]]$npn))

    mstest <- prepTestMSM(extData.list[[ext.ds]]$npn, fullGraph.list, tmat)

    msmriskdf <- msmPredict(model, mstest, tmat, c(5))

    gc()

    for (tt in c(5)) {

        msmtemp <- msmriskdf %>% filter(Time==paste(tt,'Years'))

        ids <- as.character(msmtemp$id)

        rownames(msmtemp) <- ids

        msmtemp <- cbind(extData.list[[ext.ds]]$npn[ids,], risk=-log(msmtemp[ids,c('DSS','OS','DRFS','DFS')]))

        rftemp <- extData.list[[ext.ds]]$npn[ids,]

        for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
            pred <- predict(rfModels[[sf]], newdata=extData.list[[ext.ds]]$npn[,shared.feats])
            rftemp[,paste0('risk.',sf)] <- -pred$predicted
        }


        for (erstatus in c('nostratify', 'all', 'erp', 'ern')) {

            ## if (erstatus=='ern') {
            ##     rftemp <- rftemp %>% filter(ER_STATUS=='Negative') ## %>% dplyr::select(-c('ER_STATUS'))
            ## }
            ## if (erstatus=='erp') {
            ##     rftemp <- rftemp %>% filter(ER_STATUS=='Positive') ## %>% dplyr::select(-c('ER_STATUS'))
            ## }

            for (sf in intersect(colnames(extData.list[[ext.ds]]$npn), c('DSS','OS','DRFS','DFS'))) {

                if (erstatus=='all') {
                    f <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + risk.', sf))
                    
                    deltavec <- rep(0, B)

                    for (b in 1:B) {
                        while (TRUE) {
                            ss <- sample(1:nrow(rftemp), nrow(rftemp), replace=T)
                            rf.conc.harrell <- concordance(f, rftemp[ss,], ymax=365.25*tt)

                            msm.conc.harrell <- concordance(f, msmtemp[ss,], ymax=365.25*tt)

                            deltavec[b] <- msm.conc.harrell$concordance - rf.conc.harrell$concordance
                            if (!is.na(deltavec[b])) {
                                break
                            }
                        }
                    }

                    rf.conc.harrell <- concordance(f, rftemp, ymax=365.25*tt)

                    msm.conc.harrell <- concordance(f, msmtemp, ymax=365.25*tt)


                    single <- data.frame(Dataset=ext.ds,
                                         Outcome=sf,
                                         Time=paste(tt,'Years'),
                                         `ER Status`='Stratified',
                                         `Delta Concordance` = msm.conc.harrell$concordance - rf.conc.harrell$concordance,
                                         `Delta Concordance SE` = sd(deltavec),
                                         check.names=F)

                    rfconcdiffdf <- rbind(rfconcdiffdf, single)
                    
                } else {
                    f <- as.formula(paste0(sf, ' ~ risk.', sf))

                    if (erstatus=='nostratify') {
                        deltavec <- rep(0, B)

                        for (b in 1:B) {
                            while (TRUE) {
                                ss <- sample(1:nrow(rftemp), nrow(rftemp), replace=T)
                                rf.conc.harrell <- concordance(f, rftemp[ss,], ymax=365.25*tt)

                                msm.conc.harrell <- concordance(f, msmtemp[ss,], ymax=365.25*tt)

                                deltavec[b] <- msm.conc.harrell$concordance - rf.conc.harrell$concordance
                                if (!is.na(deltavec[b])) {
                                    break
                                }
                            }
                        }

                        rf.conc.harrell <- concordance(f, rftemp, ymax=365.25*tt)

                        msm.conc.harrell <- concordance(f, msmtemp, ymax=365.25*tt)


                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='All',
                                             `Delta Concordance` = msm.conc.harrell$concordance - rf.conc.harrell$concordance,
                                             `Delta Concordance SE` = sd(deltavec),
                                             check.names=F)

                        rfconcdiffdf <- rbind(rfconcdiffdf, single)
                    }
                    if (erstatus=='erp' && any(rftemp$ER_STATUS=='Positive')) {

                        deltavec <- rep(0, B)

                        for (b in 1:B) {
                            while (TRUE) {
                                ss <- sample(1:sum(rftemp$ER_STATUS=='Positive'), sum(rftemp$ER_STATUS=='Positive'), replace=T)
                                rf.conc.harrell <- concordance(f, (rftemp %>% filter(ER_STATUS=='Positive'))[ss,], ymax=365.25*tt)

                                msm.conc.harrell <- concordance(f, (msmtemp %>% filter(ER_STATUS=='Positive'))[ss,], ymax=365.25*tt)

                                deltavec[b] <- msm.conc.harrell$concordance - rf.conc.harrell$concordance
                                if (!is.na(deltavec[b])) {
                                    break
                                }
                            }
                        }

                        rf.conc.harrell <- concordance(f, (rftemp %>% filter(ER_STATUS=='Positive')), ymax=365.25*tt)

                        msm.conc.harrell <- concordance(f, (msmtemp %>% filter(ER_STATUS=='Positive')), ymax=365.25*tt)


                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='Positive',
                                             `Delta Concordance` = msm.conc.harrell$concordance - rf.conc.harrell$concordance,
                                             `Delta Concordance SE` = sd(deltavec),
                                             check.names=F)

                        rfconcdiffdf <- rbind(rfconcdiffdf, single)
                    }

                    if (erstatus=='ern' && any(rftemp$ER_STATUS=='Negative')) {
                        deltavec <- rep(0, B)

                        for (b in 1:B) {
                            while (TRUE) {
                                ss <- sample(1:sum(rftemp$ER_STATUS=='Negative'), sum(rftemp$ER_STATUS=='Negative'), replace=T)
                                rf.conc.harrell <- concordance(f, (rftemp %>% filter(ER_STATUS=='Negative'))[ss,], ymax=365.25*tt)

                                msm.conc.harrell <- concordance(f, (msmtemp %>% filter(ER_STATUS=='Negative'))[ss,], ymax=365.25*tt)

                                deltavec[b] <- msm.conc.harrell$concordance - rf.conc.harrell$concordance
                                if (!is.na(deltavec[b])) {
                                    break
                                }
                            }
                        }

                        rf.conc.harrell <- concordance(f, (rftemp %>% filter(ER_STATUS=='Negative')), ymax=365.25*tt)

                        msm.conc.harrell <- concordance(f, (msmtemp %>% filter(ER_STATUS=='Negative')), ymax=365.25*tt)


                        single <- data.frame(Dataset=ext.ds,
                                             Outcome=sf,
                                             Time=paste(tt,'Years'),
                                             `ER Status`='Negative',
                                             `Delta Concordance` = msm.conc.harrell$concordance - rf.conc.harrell$concordance,
                                             `Delta Concordance SE` = sd(deltavec),
                                             check.names=F)

                        rfconcdiffdf <- rbind(rfconcdiffdf, single)
                    }
                }
            }
        }
    }
}

rfconcdiffdf

write.csv(rfconcdiffdf, "msm_rf_delta_conc.csv")

## msmmetaconcdf.v1 <- msmmetaconcdf

cvmsmconcdf <- read.csv('cv_msm_conc.csv', row.names=1, check.names=FALSE) %>% filter(Time %in% paste(c(5), 'Years'))
cvmsmconcdf$Dataset <- 'METABRIC'
cvmsmconcdf$Method <- 'Causal MSM'

cvrfconcdf <- read.csv('cv_rf_conc.csv', row.names=1, check.names=FALSE) %>% filter(Time %in% paste(c(5), 'Years'))
cvrfconcdf$Dataset <- 'METABRIC'
cvrfconcdf$Method <- 'RSF'

cvlassoconcdf <- read.csv('cv_lasso_conc.csv', row.names=1, check.names=FALSE) %>% filter(Time %in% paste(c(5), 'Years'))
cvlassoconcdf$Dataset <- 'METABRIC'
cvlassoconcdf$Method <- 'LASSO Cox'

jointcvdf <- rbind(cvmsmconcdf, cvrfconcdf, cvlassoconcdf)

harrell.cv.wald.df <- jointcvdf %>%
    group_by(Fold, Dataset, `ER Status`, Outcome, Time) %>%
    mutate(Delta = Concordance - Concordance[1]) %>%
    filter(Method != 'Causal MSM') %>%
    group_by(Method, Dataset, `ER Status`, Outcome, Time) %>%
    summarize(mean=mean(Delta, na.rm=T),
              se=sqrt(1/10 + 1/9) * ifelse(sd(Delta, na.rm=T)==0, 1, sd(Delta, na.rm=T))) %>%
    mutate(stat=mean/se,
           pval=2*pt(-abs(stat), 9)) %>%
    as.data.frame

harrell.meta.wald.df <- rbind(cbind(Method='LASSO Cox', lassometaconcdiffdf),
                              cbind(Method='RSF', rfmetaconcdiffdf))
colnames(harrell.meta.wald.df)[6:7] <- c('mean', 'se')

harrell.joint.wald.df <- rbind(harrell.cv.wald.df, harrell.meta.wald.df) %>% arrange(Dataset, Method, `ER Status`, Outcome)

harrell.joint.wald.df <- harrell.joint.wald.df %>%
    filter(`ER Status` != 'All') %>%
    mutate(padj=p.adjust(pval, method='fdr', n=n()),
           signif=ifelse(padj<0.05,
                  ifelse(padj<0.01,
                  ifelse(padj<0.001, '***', '**'), '*'), '')) %>%
    as.data.frame

harrell.joint.wald.df


cvmsmconcdf <- read.csv('cv_msm_conc_sum.csv', row.names=1, check.names=FALSE)

cvmsmconcdf$Dataset <- 'METABRIC'

cvmsmconcdf <- cvmsmconcdf[,colnames(msmconcdf)] %>% filter(Time %in% paste(c(5), 'Years'))

jointmsmconcdf <- rbind(msmmetaconcdf, cvmsmconcdf, msmconcdf)

jointmsmconcdf$Outcome <- factor(jointmsmconcdf$Outcome, c('DSS', 'OS', 'DRFS', 'DFS'))
jointmsmconcdf$`ER Status` <- factor(jointmsmconcdf$`ER Status`, c('All', 'Stratified', 'Positive', 'Negative'))
jointmsmconcdf$Dataset <- factor(jointmsmconcdf$Dataset, c('METABRIC', 'Meta Cohort', ext.ds.names))

jointmsmconcdf$Method <- 'Causal MSM'

ggplot(jointmsmconcdf, ## %>% filter(`ER Status`!='All'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=2) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    theme_bw()

ggsave('causal_msm_chf_conc_full.png', height=8, width=10, dpi=400)

ggplot(jointmsmconcdf %>% filter(Dataset %in% c('METABRIC', 'Meta Cohort')),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=2) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    theme_bw()

ggsave('causal_msm_chf_conc_reduced.png', height=8, width=10, dpi=400)


cvrfconcdf <- read.csv('cv_rf_conc_sum.csv', row.names=1, check.names=FALSE)

cvrfconcdf$Dataset <- 'METABRIC'

cvrfconcdf <- cvrfconcdf[,colnames(rfconcdf)] %>% filter(Time %in% paste(c(5), 'Years'))

jointrfconcdf <- rbind(rfmetaconcdf, rfconcdf, cvrfconcdf)

jointrfconcdf$Outcome <- factor(jointrfconcdf$Outcome, c('DSS', 'OS', 'DRFS', 'DFS'))
jointrfconcdf$`ER Status` <- factor(jointrfconcdf$`ER Status`, c('All', 'Stratified', 'Positive', 'Negative'))
jointrfconcdf$Dataset <- factor(jointrfconcdf$Dataset, c('METABRIC', 'Meta Cohort', ext.ds.names))

jointrfconcdf$Method <- 'RSF'

cvlassoconcdf <- read.csv('cv_lasso_conc_sum.csv', row.names=1, check.names=FALSE)

cvlassoconcdf$Dataset <- 'METABRIC'

cvlassoconcdf <- cvlassoconcdf[,colnames(lassoconcdf)] %>% filter(Time %in% paste(c(5), 'Years'))

jointlassoconcdf <- rbind(lassometaconcdf, lassoconcdf, cvlassoconcdf)

jointlassoconcdf$Outcome <- factor(jointlassoconcdf$Outcome, c('DSS', 'OS', 'DRFS', 'DFS'))
jointlassoconcdf$`ER Status` <- factor(jointlassoconcdf$`ER Status`, c('All', 'Stratified', 'Positive', 'Negative'))
jointlassoconcdf$Dataset <- factor(jointlassoconcdf$Dataset, c('METABRIC', 'Meta Cohort', ext.ds.names))

jointlassoconcdf$Method <- 'LASSO Cox'



ggplot(rbind(jointmsmconcdf, jointrfconcdf, jointlassoconcdf) %>% filter(`ER Status`!='All'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset, linetype=Method, shape=Method)) +
    geom_pointrange(position=position_dodge(width=0.85), fatten=2) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    theme_bw()

ggsave('all_method_conc_full.png', height=8, width=18, dpi=400)


rbind(jointmsmconcdf, jointrfconcdf, jointlassoconcdf) %>% filter(Dataset %in% c('Meta Cohort'))

0.7597843 - 0.8050620
0.7880550 - 0.7968438

ggplot(rbind(jointmsmconcdf, jointrfconcdf, jointlassoconcdf) %>% filter(Dataset %in% c('METABRIC', 'Meta Cohort'), `ER Status`!='All'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset, linetype=Method, shape=Method)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=2) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    theme_bw()

ggsave('all_method_conc_reduced.png', height=8, width=10, dpi=400)


ggplot(rbind(jointmsmconcdf, jointrfconcdf, jointlassoconcdf) %>%
       filter(`ER Status`!='All', Time=='5 Years'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset, linetype=Method, shape=Method)) +
    geom_pointrange(position=position_dodge(width=0.85), fatten=2) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    theme_bw()

ggsave('all_method_5yr_conc_full.png', height=6, width=9, dpi=400)

jointsummary <- rbind(jointmsmconcdf, jointrfconcdf, jointlassoconcdf) %>% filter(Dataset %in% c('METABRIC', 'Meta Cohort'), `ER Status`!='All', Time=='5 Years')%>% arrange(Dataset, Method, `ER Status`, Outcome) %>% as.data.frame

harrell.joint.wald.df$Outcome <- factor(harrell.joint.wald.df$Outcome, c('DSS', 'OS', 'DRFS', 'DFS'))
harrell.joint.wald.df$`ER Status` <- factor(harrell.joint.wald.df$`ER Status`, c('All', 'Stratified', 'Positive', 'Negative'))
harrell.joint.wald.df$Dataset <- factor(harrell.joint.wald.df$Dataset, c('METABRIC', 'Meta Cohort', ext.ds.names))


harrell.joint.wald.df <- harrell.joint.wald.df %>% filter(Dataset %in% c('METABRIC', 'Meta Cohort'), `ER Status`!='All', Time=='5 Years')%>% arrange(Dataset, Method, `ER Status`, Outcome) %>% as.data.frame

all(jointsummary[13:36,colnames(harrell.joint.wald.df)[1:5]] == harrell.joint.wald.df[1:24,1:5])

jointsummary[13:36,c('pval', 'padj', 'signif')] <- harrell.joint.wald.df[1:24,c('pval', 'padj', 'signif')]

all(jointsummary[49:72,colnames(harrell.joint.wald.df)[1:5]] == harrell.joint.wald.df[25:48,1:5])

jointsummary[49:72,c('pval', 'padj', 'signif')] <- harrell.joint.wald.df[25:48,c('pval', 'padj', 'signif')]

jointsummary$Method <- factor(jointsummary$Method, levels=c('Causal MSM', 'LASSO Cox', 'RSF'))


ggplot(jointsummary %>% filter(Dataset %in% c('METABRIC', 'Meta Cohort'), `ER Status`!='All', Time=='5 Years'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset, linetype=Method, shape=Method)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=2) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    ylim(0.45, 1.05) +
    geom_bracket(aes(xmax=1.07 + 0.14 * (as.numeric(Method)-1), xmin=1.07, y.position=Concordance + 3 * `Concordance SE`+0.05, color='black', label=ifelse(padj < 1e-4, 'p < 1e-04', paste('p =', round(padj,4))), color='black'),
             jointsummary %>% filter(padj < 0.05) %>% group_by(Outcome, `ER Status`) %>% mutate_at(c('Concordance', 'Concordance SE'), function(x) max(x)), step.increase=0.015, label.size=1.8, tip.length=0.01) + 
    theme_bw(base_size=7)

ggsave('all_method_5yr_conc_reduced_signif.png', height=4.5, width=5.25, dpi=400)


ggplot(jointsummary %>% filter(Dataset %in% c('METABRIC', 'Meta Cohort'), `ER Status`=='Stratified', Time=='5 Years'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset, shape=Method)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=4) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    ylim(0.45, 1.05) +
    geom_bracket(aes(xmax=1.07 + 0.14 * (as.numeric(Method)-1), xmin=1.07, y.position=Concordance + 3 * `Concordance SE`+0.05, color='black', label=ifelse(padj < 1e-4, 'p < 1e-04', paste('p =', round(padj,4))), color='black'),
             jointsummary %>% filter(padj < 0.05, `ER Status`=='Stratified') %>% group_by(Outcome, `ER Status`) %>% mutate_at(c('Concordance', 'Concordance SE'), function(x) max(x)), step.increase=0.025, label.size=1.8, tip.length=0.01) + 
    theme_bw(base_size=7)

ggsave('all_method_5yr_conc_reduced_stratified_signif.png', height=3, width=5, dpi=600)


ggplot(jointsummary %>% filter(Dataset %in% c('METABRIC', 'Meta Cohort'), `ER Status` %in% c('Positive', 'Negative'), Time=='5 Years'),
       aes(x=Time, y=`Concordance`,
           ymin=Concordance-1.96*`Concordance SE`,
           ymax=Concordance+1.96*`Concordance SE`,
           color=Dataset, shape=Method)) +
    geom_pointrange(position=position_dodge(width=0.8), fatten=4) +
    geom_hline(yintercept=0.5, color='darkgrey', linetype=5) +
    facet_grid(cols=vars(Outcome), rows=vars(`ER Status`)) +
    ylim(0.45, 1.05) +
    geom_bracket(aes(xmax=1.07 + 0.14 * (as.numeric(Method)-1), xmin=1.07, y.position=Concordance + 3 * `Concordance SE`+0.05, color='black', label=ifelse(padj < 1e-4, 'p < 1e-04', paste('p =', round(padj,4))), color='black'),
             jointsummary %>% filter(padj < 0.05, `ER Status` %in% c('Positive', 'Negative')) %>% group_by(Outcome, `ER Status`) %>% mutate_at(c('Concordance', 'Concordance SE'), function(x) max(x)), step.increase=0.025, label.size=1.8, tip.length=0.01) + 
    theme_bw(base_size=7)

ggsave('all_method_5yr_conc_reduced_ersplit_signif.png', height=4, width=5, dpi=600)


for (ds.name in ext.ds.names) {
    if ('DRFS' %in% colnames(extData.list[[ds.name]])) {
        temp <- extData.list[[ds.name]][,c('DRFS', 'ER_STATUS', 'VCAM1')]
        print(prop.table(table(temp$ER_STATUS)))
        if (prop.table(table(temp$ER_STATUS))['Positive'] < 0.9) {
            temp <- temp %>% filter(ER_STATUS=='Negative')
            temp$DRFS.time <- temp$DRFS[,1]
            temp$DRFS.status <- temp$DRFS[,2]
            res.cut <- surv_cutpoint(temp, 'DRFS.time', 'DRFS.status', variables='VCAM1')
            print(ds.name)
            print(plot(res.cut))
            ## Sys.sleep(2)
            res.cat <- surv_categorize(res.cut)
            png(paste0('VCAM1_DRFS_', ds.name, '_ern.png'), height=5, width=5, units='in', res=400)
            print(ggsurvplot(survfit(Surv(DRFS.time, DRFS.status) ~ VCAM1, res.cat), pval=T, conf.int=T))
            dev.off()
            ## Sys.sleep(2)
        }
    }
}


for (ds.name in ext.ds.names) {
    if ('OS' %in% colnames(extData.list[[ds.name]])) {
        temp <- extData.list[[ds.name]][,c('OS', 'ER_STATUS', 'VCAM1')]
        if (prop.table(table(temp$ER_STATUS))['Positive'] < 0.9) {
            temp <- temp %>% filter(ER_STATUS=='Negative')
            temp$OS.time <- temp$OS[,1]
            temp$OS.status <- temp$OS[,2]
            res.cut <- surv_cutpoint(temp, 'OS.time', 'OS.status', variables='VCAM1')
            print(ds.name)
            print(plot(res.cut))
            ## Sys.sleep(2)
            res.cat <- surv_categorize(res.cut)
            png(paste0('VCAM1_OS_', ds.name, '_ern.png'), height=5, width=5, units='in', res=400)
            print(ggsurvplot(survfit(Surv(OS.time, OS.status) ~ VCAM1, res.cat), pval=T, conf.int=T))
            dev.off()
            ## Sys.sleep(2)
        }
    }
}


for (ds.name in ext.ds.names) {
    if ('DSS' %in% colnames(extData.list[[ds.name]])) {
        temp <- extData.list[[ds.name]][,c('DSS', 'ER_STATUS', 'VCAM1')]
        if (prop.table(table(temp$ER_STATUS))['Positive'] < 0.9) {
            temp <- temp %>% filter(ER_STATUS=='Negative')
            temp$DSS.time <- temp$DSS[,1]
            temp$DSS.status <- temp$DSS[,2]
            res.cut <- surv_cutpoint(temp, 'DSS.time', 'DSS.status', variables='VCAM1')
            print(ds.name)
            print(plot(res.cut))
            ## Sys.sleep(2)
            res.cat <- surv_categorize(res.cut)
            png(paste0('VCAM1_DSS_', ds.name, '_ern.png'), height=5, width=5, units='in', res=400)
            print(ggsurvplot(survfit(Surv(DSS.time, DSS.status) ~ VCAM1, res.cat), pval=T, conf.int=T))
            dev.off()
            ## Sys.sleep(2)
        }
    }
}


for (ds.name in ext.ds.names) {
    if ('DFS' %in% colnames(extData.list[[ds.name]])) {
        temp <- extData.list[[ds.name]][,c('DFS', 'ER_STATUS', 'VCAM1')]
        print(prop.table(table(temp$ER_STATUS)))
        if (prop.table(table(temp$ER_STATUS))['Positive'] < 0.9) {
            temp <- temp %>% filter(ER_STATUS=='Negative')
            temp$DFS.time <- temp$DFS[,1]
            temp$DFS.status <- temp$DFS[,2]
            res.cut <- surv_cutpoint(temp, 'DFS.time', 'DFS.status', variables='VCAM1')
            print(ds.name)
            print(plot(res.cut))
            ## Sys.sleep(2)
            res.cat <- surv_categorize(res.cut)
            png(paste0('VCAM1_DFS_', ds.name, '_ern.png'), height=5, width=5, units='in', res=400)
            print(ggsurvplot(survfit(Surv(DFS.time, DFS.status) ~ VCAM1, res.cat), pval=T, conf.int=T))
            dev.off()
            ## Sys.sleep(2)
        }
    }
}


joint.clin.df <- data.frame(
    Dataset='Metabric',
    fullData.list[['all']]$train[,c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')])
## joint.clin.df <- rbind(joint.clin.df,
##                        data.frame(Dataset="GSE11121",
##                                   AGE_AT_DIAGNOSIS=NA,
##                                   TUMOR_SIZE=extData.list[["GSE11121"]][,c('TUMOR_SIZE')],
##                                   row.names=rownames(extData.list[["GSE11121"]])))
for (ds.name in ext.ds.names) {
    joint.clin.df <- rbind(
        joint.clin.df,
        data.frame(Dataset=ds.name,
                   extData.list[[ds.name]][,c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')]))
}
joint.clin.df <- joint.clin.df %>% mutate_at(c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE'), as.numeric)

ggplot(joint.clin.df, aes(y=AGE_AT_DIAGNOSIS, x=Dataset, fill=Dataset)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    theme_classic()

ggplot(joint.clin.df, aes(y=TUMOR_SIZE, x=Dataset, fill=Dataset)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    theme_classic()

joint.clin.df$AGE_AT_DIAGNOSIS[!is.na(joint.clin.df$AGE_AT_DIAGNOSIS)] <- huge::huge.npn(as.matrix(joint.clin.df$AGE_AT_DIAGNOSIS[!is.na(joint.clin.df$AGE_AT_DIAGNOSIS)]))
joint.clin.df$TUMOR_SIZE <- huge::huge.npn(as.matrix(joint.clin.df$TUMOR_SIZE))

ggplot(joint.clin.df, aes(y=AGE_AT_DIAGNOSIS, x=Dataset, fill=Dataset)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    stat_compare_means() +
    theme_classic()

ggplot(joint.clin.df, aes(y=TUMOR_SIZE, x=Dataset, fill=Dataset)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    stat_compare_means() +
    theme_classic()

fullData.list[['all']]$train.npn[,c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')] <- joint.clin.df[rownames(fullData.list[['all']]$train.npn),c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')]

fullData.list[['erp']]$train.npn[,c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')] <- joint.clin.df[rownames(fullData.list[['erp']]$train.npn),c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')]

fullData.list[['ern']]$train.npn[,c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')] <- joint.clin.df[rownames(fullData.list[['ern']]$train.npn),c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')]

for (ds.name in ext.ds.names) {
    extData.list[[ds.name]][,c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')] <- joint.clin.df[rownames(extData.list[[ds.name]]),c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')]
    if (ds.name=="GSE11121") {
        extData.list[[ds.name]]$AGE_AT_DIAGNOSIS <- 0
    }
}

fullModels.list <- lapply(list(all='all', erp='erp', ern='ern'),
                          function(x) trainModels(fullData.list[[x]]$train.npn,
                                                  fullGraph.list[[x]], x))

lapply(fullModels.list, FUN=lapply, summary)

ext.feats <- lapply(extData.list, colnames)
shared.feats <- colnames(fullData.list$all$train)
for (idx in 1:length(ext.feats)) {
    shared.feats <- intersect(shared.feats, ext.feats[[idx]])
}
shared.feats

lassoModels.list <- lapply(list(all='all', erp='erp', ern='ern'),
                           function(x) trainLassoModels(fullData.list[[x]]$train.npn,
                                                        shared.feats, x))

erp.risk <- matrix(NA, 1815, 7)
ern.risk <- matrix(NA, 1815, 7)
colnames(erp.risk) <- paste0('risk.',surv.feats)
colnames(ern.risk) <- paste0('risk.',surv.feats)
rownames(erp.risk) <- rownames(fullData.list$all$train.npn)
rownames(ern.risk) <- rownames(fullData.list$all$train.npn)
for (sf in surv.feats) {
    if (sf %in% surv.feats[1:4]) {
    erp.risk[,paste0('risk.',sf)] <- predict(fullModels.list[['erp']][[sf]],
                                             newdata=fullData.list$all$train.npn)
    ern.risk[,paste0('risk.',sf)] <- predict(fullModels.list[['ern']][[sf]],
                                             newdata=fullData.list$all$train.npn)
    } else {
        erp.risk[,paste0('risk.',sf)] <- predict(fullModels.list[['erp']][[sf]],
                                             newdata=data.frame(erp.risk))
        ern.risk[,paste0('risk.',sf)] <- predict(fullModels.list[['ern']][[sf]],
                                                 newdata=data.frame(ern.risk))
    
    }
}

erp.risk.erp <- erp.risk[rownames(fullData.list[['erp']]$train),]
ern.risk.ern <- ern.risk[rownames(fullData.list[['ern']]$train),]

library(CCA)

gene.mask <- sapply(fullData.list[['erp']]$train, function(x) !is.Surv(x) & is.numeric(x))
gene.mask[c('AGE_AT_DIAGNOSIS', 'TUMOR_SIZE')] <- FALSE

cv.rcc <- estim.regul(risk.scores, scale(fullData.list[['all']]$train[,gene.mask]),
                      10^seq(-2, 0, length.out=5), 10^seq(-2, 0, length.out=5), plt=TRUE)

risk.scores <- data.frame(erp.risk[,1:4], ern.risk[,1:4])
colnames(risk.scores) <- c(paste0('erp.risk.',surv.feats[1:4]), paste0('ern.risk.',surv.feats[1:4]))

res <- rcc(as.matrix(risk.scores[rownames(fullData.list[['erp']]$train),1:4]), scale(fullData.list[['erp']]$train[,gene.mask]), 0.25, 0.25)

plt.cc(res, var.label=T)

barplot(res$cor, xlab = "Dimension", ylab = "Canonical correlations", names.arg = 1:8, ylim = c(0,1))

pred.erp <- as.matrix(scale(fullData.list[['erp']]$train[,gene.mask])) %*% res$ycoef %*% diag(res$cor) %*% t(res$xcoef)

cor(erp.risk.erp, pred.erp)

plot(scale(fullData.list[['erp']]$train[,gene.mask]) %*% res$ycoef %*% sqrt(diag(1/res$cor)),
     erp.risk.erp[,1:4] %*% res$xcoef %*% sqrt(diag(1/res$cor)))

plot(erp.risk.erp[,1],
     (scale(fullData.list[['erp']]$train[,gene.mask]) %*% res$ycoef %*% t(res$xcoef))[,1])

cov(scale(fullData.list[['erp']]$train[,gene.mask]) %*% res$ycoef %*% diag(1/sqrt(res$cor)),
    erp.risk.erp[,1:4] %*% res$xcoef %*% diag(1/sqrt(res$cor)))

pred.ern <- scale(fullData.list[['ern']]$train[,gene.mask]) %*% res$ycoef %*% diag(1/(res$cor)) %*% t(res$xcoef)

cov(risk.scores[rownames(fullData.list[['ern']]$train),1:4], pred.ern)

for (i in 1:4) {
    plot(risk.scores[rownames(fullData.list[['ern']]$train),i],
         pred.ern[,i], main=surv.feats[i], xlab='True', ylab='Ridge Predict')
    print(cor.test(risk.scores[rownames(fullData.list[['ern']]$train),i],
                   pred.ern[,i]))
}

par(mfcol=c(4,2))
for (i in 1:4) {
    cv.out <- cv.glmnet(scale(fullData.list[['erp']]$train[,gene.mask]),
                        erp.risk.erp[,i],
                        family='gaussian', alpha=0)
    ## plot(cv.out)
    plot(erp.risk[rownames(fullData.list[['ern']]$train),i],
         predict(cv.out, scale(fullData.list[['ern']]$train[,gene.mask]),
                 s=cv.out$lambda.min), main=surv.feats[i], xlab='True', ylab='Ridge Predict')
    print(cor.test(erp.risk[rownames(fullData.list[['ern']]$train),i],
                   predict(cv.out, scale(fullData.list[['ern']]$train[,gene.mask]),
                           s=cv.out$lambda.min)))
}

par(mfcol=c(2,7))
for (i in 1:7) {
    cv.out <- cv.glmnet(scale(fullData.list[['ern']]$train[,gene.mask]),
                        ern.risk.ern[,i],
                        family='gaussian', alpha=0)
    plot(cv.out)
    plot(ern.risk[rownames(fullData.list[['erp']]$train),i],
         predict(cv.out, scale(fullData.list[['erp']]$train[,gene.mask]),
                 s=cv.out$lambda.min), main=surv.feats[i], xlab='True', ylab='Ridge Predict')
    print(cor.test(ern.risk[rownames(fullData.list[['erp']]$train),i],
                   predict(cv.out, scale(fullData.list[['erp']]$train[,gene.mask]),
                           s=cv.out$lambda.min)))
}


cv.out <- cv.glmnet(scale(fullData.list[['erp']]$train[,gene.mask]),
                    erp.risk.erp[,1],
                    family='gaussian', alpha=0)

tierList <- list(
    c('AGE_AT_DIAGNOSIS'), c('INFERRED_MENOPAUSAL_STATE'),
    colnames(fullData.list[['ern']]$train)[!colnames(fullData.list[['ern']]$train) %in%
                                           c('AGE_AT_DIAGNOSIS',
                                             'INFERRED_MENOPAUSAL_STATE',
                                             surv.feats)],
    surv.feats)

knowledge.ern <- list(tiers=tierList,
                      forbiddenWithinTier=c(FALSE, FALSE, FALSE, TRUE),
                      ## forbidden=list(c('HER2_STATUS', 'ERBB2'),
                      ##                c('PR_STATUS', 'PGR')),
                      required=list(c('ERBB2', 'HER2_STATUS'),
                                    c('PGR', 'PR_STATUS')))

saveRDS(knowledge.ern, 'knowledge_subtype.rds')

ig.path.ern <- coxmgmPath(fullData.list[['ern']]$train, nLambda=20, rank=T, verbose=T)

plot(ig.path.ern)

gamma <- 0.5

penalty <- 1 + gamma * log(ncol(fullData.list[['ern']]$train)) / log(nrow(fullData.list[['ern']]$train))

penalty

ebic <- -2 * ig.path.ern$loglik + penalty * log(nrow(fullData.list[['ern']]$train)) * ig.path.ern$nParams

points(log10(ig.path.ern$lambdas), ebic / (2 * nrow(fullData.list[['ern']]$train)), col='purple', pch=19)

ig.path.ern$graph.bic

ig.path.ern$graph.ebic <- ig.path.ern$graphs[[which.min(ebic)]]

ig.path.ern$graph.ebic

ig.path.ern$graph.ebic$markov.blankets[surv.feats]

g.ern.v2 <- fciStable(fullData.list[['ern']]$train,
                      initialGraph=ig.path.ern$graph.ebic,
                      knowledge=knowledge.ern, orientRule='maxp',
                      alpha=0.1, fdr=T, rank=T, verbose=T)

g.ern

g.ern.v2

g.ern.v3

g.ern.v4

g.ern.v6

g.ern.v5

g.ern$markov.blankets[surv.feats]

lapply(surv.feats, function(x) g.ern$edges[grepl(paste0('^',x,' | ',x,'$'), g.ern$edges)])

lapply(surv.feats, function(x) g.ern.v2$edges[grepl(paste0('^',x,' | ',x,'$'), g.ern.v2$edges)])

g.ern.v2$edges[grepl('CODE', g.ern.v2$edges)]

lapply(surv.feats, function(x) g.ern.v3$edges[grepl(paste0('^',x,' | ',x,'$'), g.ern.v3$edges)])

g.ern.v3$edges[grepl('CODE', g.ern.v3$edges)]

unique(unlist(lapply(surv.feats, getNeighbors, graph=g.ern.v3)))

plot.graph(g.ern, unique(c(surv.feats, unlist(lapply(surv.feats, getNeighbors, graph=g.ern)))), list(fontsize=48))


tierList <- list(
    c('AGE_AT_DIAGNOSIS'), c('INFERRED_MENOPAUSAL_STATE'),
    c('ONCOTREE_CODE'),
    colnames(fullData.list[['erp']]$train)[!colnames(fullData.list[['erp']]$train) %in%
                                           c('AGE_AT_DIAGNOSIS',
                                             'INFERRED_MENOPAUSAL_STATE',
                                             'ONCOTREE_CODE',
                                             surv.feats)],
    surv.feats)

knowledge.erp <- list(tiers=tierList,
                      forbiddenWithinTier=c(FALSE, FALSE, FALSE, FALSE, TRUE),
                      forbidden=list(c('HER2_STATUS', 'ERBB2'),
                                     c('PR_STATUS', 'PGR')),
                      required=list(c('ERBB2', 'HER2_STATUS'),
                                    c('PGR', 'PR_STATUS')))

ig.path.erp <- coxmgmPath(fullData.list[['erp']]$train, nLambda=20, rank=T, verbose=T)

plot(ig.path.erp)

ig.path.erp$graph.bic

ig.path.erp$graph.bic$markov.blankets[surv.feats]

g.erp.v3 <- fciStable(fullData.list[['erp']]$train,
                      initialGraph=ig.path.erp$graph.bic,
                      knowledge=knowledge.erp, orientRule='maxp',
                      alpha=0.1, fdr=T, rank=T, verbose=T)

g.erp

g.erp.v3

g.erp.v2$ambiguous_triples

g.erp$markov.blankets[surv.feats]

g.erp.v2$markov.blankets[surv.feats]

lapply(surv.feats, function(x) g.erp$edges[grepl(paste0('^',x,' | ',x,'$'), g.erp$edges)])

lapply(surv.feats, function(x) g.erp.v3$edges[grepl(paste0('^',x,' | ',x,'$'), g.erp.v3$edges)])

unique(unlist(lapply(surv.feats, getNeighbors, graph=g.erp)))

plot.graph(g.erp, unique(c(surv.feats, unlist(lapply(surv.feats, getNeighbors, graph=g.erp)))), list(fontsize=48))

tierList <- list(
    c('AGE_AT_DIAGNOSIS'), c('INFERRED_MENOPAUSAL_STATE'),
    colnames(fullData.list[['all']]$train)[!colnames(fullData.list[['all']]$train) %in%
                                           c('AGE_AT_DIAGNOSIS',
                                             'INFERRED_MENOPAUSAL_STATE',
                                             surv.feats)],
    surv.feats)

knowledge.all <- list(tiers=tierList,
                      forbiddenWithinTier=c(FALSE, FALSE, FALSE, TRUE),
                      forbidden=list(c('ER_STATUS', 'DSS'),
                                     c('ER_STATUS', 'OD'),
                                     c('ER_STATUS', 'DistantRelapse'),
                                     c('ER_STATUS', 'LocalRelapse'),
                                     c('ER_STATUS', 'OS'),
                                     c('ER_STATUS', 'DRFS'),
                                     c('ER_STATUS', 'DFS')),
                      required=list(c('ESR1', 'ER_STATUS'),
                                    c('ERBB2', 'HER2_STATUS'),
                                    c('PGR', 'PR_STATUS')))

saveRDS(knowledge.all, 'knowledge_full.rds')

ig.path.all <- coxmgmPath(fullData.list[['all']]$train, nLambda=20, rank=T, verbose=T)

plot(ig.path.all)

ig.path.all$graph.bic

ig.path.all$graph.bic$markov.blankets[surv.feats]

plot(ig.path.all)

gamma <- 0.5

penalty <- 1 + gamma * log(ncol(fullData.list[['all']]$train)) / log(nrow(fullData.list[['all']]$train))

penalty

ebic <- -2 * ig.path.all$loglik + penalty * log(nrow(fullData.list[['all']]$train)) * ig.path.all$nParams

points(log10(ig.path.all$lambdas), ebic / (2 * nrow(fullData.list[['all']]$train)), col='purple', pch=19)

ig.path.all$graph.bic

ig.path.all$graph.ebic <- ig.path.all$graphs[[which.min(ebic)]]

ig.path.all$graph.ebic

ig.path.all$graph.ebic$markov.blankets[surv.feats]

cph.ebic <- lapply(ig.path.all$graph.ebic$markov.blankets[surv.feats],
                   function(x) coxph(train.npn$DSS ~ strata(train.npn$ER_STATUS) + .,
                                     train.npn[,x]))

cph.ebic.dir <- lapply(cph.ebic, function(x) ifelse(coef(x)>0, 'UP', 'DOWN'))

paste(names(cph.ebic.dir$DSS)[cph.ebic.dir$DSS=="UP"],collapse=' ')

paste(names(cph.ebic.dir$DSS)[cph.ebic.dir$DSS=="DOWN"],collapse=' ')

g.all <- fciStable(fullData.list[['all']]$train,
                   initialGraph=ig.path.all$graph.ebic,
                   knowledge=knowledge.all, orientRule='maxp',
                   alpha=0.1, fdr=T, rank=T, verbose=T)

g.all

g.all$markov.blankets[surv.feats]

lapply(surv.feats, function(x) g.all$edges[grepl(paste0('^',x,' | ',x,'$'), g.all$edges)])

g.all.v3

g.all.v3$markov.blankets[surv.feats]

lapply(surv.feats, function(x) g.all.v3$edges[grepl(paste0('^',x,' | ',x,'$'), g.all.v3$edges)])

lapply(surv.feats, function(x) summary(coxph(as.formula(paste(x, '~ strata(train.npn$ER_STATUS) + .')), train.npn[,c(x,sort(g.all$markov.blankets[[x]]))])))


prog.feats.all <- unique(unlist(lapply(surv.feats, getNeighbors, graph=g.all.v3)))

plot.graph(g.all, unique(c(surv.feats, unlist(lapply(surv.feats, getNeighbors, graph=g.all)))), list(fontsize=72))


fullData.list <- lapply(c('all', 'ern', 'erp'), loadGraphAndData, k=-1)
names(fullData.list) <- c('all', 'ern', 'erp')

paste(unique(unlist(lapply('DRFS', getNeighbors, graph=fullData.list[['all']]$graph$BIC))), collapse=' + ')

paste(unique(unlist(lapply(surv.feats[1:3], getNeighbors, graph=fullData.list[['all']]$graph$BIC))), collapse=' + ')

summary(coxph(DRFS ~ strata(ER_STATUS) + AGE_AT_DIAGNOSIS + IGFBP5 + INFERRED_MENOPAUSAL_STATE + NQO1 + PRAME + SOX11 + LYMPH_NODE_STATUS + MYBPC1 + TUMOR_SIZE, train.npn))

paste(coef(coxph(DRFS ~ strata(ER_STATUS) + IGFBP5 + NQO1 + S100P + SERPINA1 + UBE2C_grp + LYMPH_NODE_STATUS + TUMOR_SIZE + AGE_AT_DIAGNOSIS + MMP11 + MYBPC1, train.npn)), collapse=', ')

ext.data.list <- lapply(ext.ds.names, loadExtData, feats=colnames(fullData$train))
names(ext.data.list) <- ext.ds.names

model.list <- list()

train.npn <- fullData.list[['all']]$train %>%
    mutate_if(function(x) !is.Surv(x) & is.numeric(x),
              function(x) as.vector(huge::huge.npn(as.matrix(x))))

train.npn$LYMPH_NODE_STATUS <- ifelse(train.npn$LYMPH_NODE_STATUS=='NodeNegative',
                                      'Negative', 'Positive')

for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
    model.list[[sf]] <- list()
    if (sf %in% surv.feats[5:7]) {
        nb.all <- setdiff(
            sort(
                unique(
                    unlist(
                        lapply(surv.feats[1:(which(surv.feats==sf)-3)], getNeighbors, graph=fullData.list[['all']]$graph[['StEPS']])
                    )
                )
            ),
            surv.feats)
        
        mb.bic.all <- setdiff(
            sort(
                unique(
                    unlist(
                        lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                               function(x) growShrinkMB(fullData.list[['all']]$train, x, graph=fullData.list[['all']]$graph[['StEPS']], rank=T, verbose=T))
                    )
                )
            ),
            surv.feats)

        nb.ern <- setdiff(
            sort(
                unique(
                    unlist(
                        lapply(surv.feats[1:(which(surv.feats==sf)-3)], getNeighbors, graph=fullData.list[['ern']]$graph[['StARS']])
                    )
                )
            ),
            surv.feats)
        
        mb.bic.ern <- setdiff(
            sort(
                unique(
                    unlist(
                        lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                               function(x) growShrinkMB(fullData.list[['ern']]$train, x, graph=fullData.list[['ern']]$graph[['StARS']], rank=T, verbose=T))
                    )
                )
            ),
            surv.feats)

        nb.erp <- setdiff(
            sort(
                unique(
                    unlist(
                        lapply(surv.feats[1:(which(surv.feats==sf)-3)], getNeighbors, graph=fullData.list[['erp']]$graph[['StARS']])
                    )
                )
            ),
            surv.feats)
        
        mb.bic.erp <- setdiff(
            sort(
                unique(
                    unlist(
                        lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                               function(x) growShrinkMB(fullData.list[['erp']]$train, x, graph=fullData.list[['erp']]$graph[['StARS']], rank=T, verbose=T))
                    )
                )
            ),
            surv.feats)

        
    } else {
        nb.all <- sort(getNeighbors(fullData.list[['all']]$graph[['StEPS']], sf))
        mb.bic.all <- sort(growShrinkMB(fullData.list[['all']]$train, sf, fullData.list[['all']]$graph[['StEPS']], rank=T, verbose=T))

        nb.ern <- sort(getNeighbors(fullData.list[['ern']]$graph[['StARS']], sf))
        mb.bic.ern <- sort(growShrinkMB(fullData.list[['ern']]$train, sf, fullData.list[['ern']]$graph[['StARS']], rank=T, verbose=T))

        nb.erp <- sort(getNeighbors(fullData.list[['erp']]$graph[['StARS']], sf))
        mb.bic.erp <- sort(growShrinkMB(fullData.list[['erp']]$train, sf, fullData.list[['erp']]$graph[['StARS']], rank=T, verbose=T))
    }

    f.nb.all <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(nb.all, collapse=' + ')))
    f.mb.bic.all <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(mb.bic.all, collapse=' + ')))

    f.nb.ern <- as.formula(paste0(sf, ' ~ ', paste(nb.ern, collapse=' + ')))
    f.mb.bic.ern <- as.formula(paste0(sf, ' ~ ', paste(mb.bic.ern, collapse=' + ')))

    f.nb.erp <- as.formula(paste0(sf, ' ~ ', paste(nb.erp, collapse=' + ')))
    f.mb.bic.erp <- as.formula(paste0(sf, ' ~ ', paste(mb.bic.erp, collapse=' + ')))

    cph.nb.all <- coxph(f.nb.all, train.npn)
    cph.mb.bic.all <- coxph(f.mb.bic.all, train.npn)

    model.list[[sf]][['nb.all']] <- cph.nb.all
    model.list[[sf]][['mb.bic.all']] <- cph.mb.bic.all

    cph.nb.ern <- coxph(f.nb.ern, train.npn[train.npn$ER_STATUS=='Negative',])
    cph.mb.bic.ern <- coxph(f.mb.bic.ern, train.npn[train.npn$ER_STATUS=='Negative',])

    model.list[[sf]][['nb.ern']] <- cph.nb.ern
    model.list[[sf]][['mb.bic.ern']] <- cph.mb.bic.ern

    cph.nb.erp <- coxph(f.nb.erp, train.npn[train.npn$ER_STATUS=='Positive',])
    cph.mb.bic.erp <- coxph(f.mb.bic.erp, train.npn[train.npn$ER_STATUS=='Positive',])

    model.list[[sf]][['nb.erp']] <- cph.nb.erp
    model.list[[sf]][['mb.bic.erp']] <- cph.mb.bic.erp

    cv.lasso.all <- glmnet::cv.glmnet(model.matrix(~., train.npn[,!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], train.npn[,sf], family='cox')

    plot(cv.lasso.all)

    model.list[[sf]][['lasso.all']] <- cv.lasso.all

    cv.lasso.ern <- glmnet::cv.glmnet(model.matrix(~., train.npn[train.npn$ER_STATUS=='Negative',!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], train.npn[train.npn$ER_STATUS=='Negative',sf], family='cox')

    plot(cv.lasso.ern)

    model.list[[sf]][['lasso.ern']] <- cv.lasso.ern

    cv.lasso.erp <- glmnet::cv.glmnet(model.matrix(~., train.npn[train.npn$ER_STATUS=='Positive',!colnames(train.npn) %in% c('ER_STATUS', surv.feats)])[,-1], train.npn[train.npn$ER_STATUS=='Positive',sf], family='cox')

    plot(cv.lasso.erp)

    model.list[[sf]][['lasso.erp']] <- cv.lasso.erp
    
}

library(survcomp)
ext.conc.df <- data.frame()
for (sf in c('DSS', 'OS', 'DFS', 'DRFS')) {
    outcome.eval.df <- evaluateModels(temp, extData.list, c(2, 5, 8), sf)

    meta.eval.df <- data.frame()
    for (t in c(2, 5, 8)) {
        for (mod in unique(outcome.eval.df$Model)) {
            for (erstatus in c('All', 'ER+', 'ER-')) {
                single.eval.df <- outcome.eval.df %>% filter(Time==paste(t, 'Years'),
                                                             Model==mod,
                                                             `ER Status`==erstatus)
                single.eval.df <- na.omit(single.eval.df)
                meta.uno <- combine.est(single.eval.df$`Uno Concordance`,
                                        ifelse(single.eval.df$`Uno SE`==0, 0.15,
                                               single.eval.df$`Uno SE`))
                meta.harrell <- combine.est(single.eval.df$`Harrell Concordance`,
                                            ifelse(single.eval.df$`Harrell SE`==0, 0.15,
                                                   single.eval.df$`Harrell SE`))

                meta.eval.df <- rbind(meta.eval.df,
                                      data.frame(`Uno Concordance` = meta.uno$estimate,
                                                 `Harrell Concordance` = meta.harrell$estimate,
                                                 `Uno SE` = meta.uno$se,
                                                 `Harrell SE` = meta.harrell$se,
                                                 Model = mod,
                                                 `ER Status` = erstatus,
                                                 Time = paste(t, 'Years'),
                                                 Outcome = sf,
                                                 Dataset = 'Meta Cohort',
                                                 check.names=F))
            }
        }
    }

    ext.conc.df <- rbind(ext.conc.df,
                         outcome.eval.df,
                         meta.eval.df)
}

ext.conc.df <- evaluateCausalModels(fullModels.list, extData.list, times=c(3, 5, 8))

ext.conc.df <- rbind(
    ext.conc.df,
    evaluateLassoModels(lassoModels.list, extData.list, shared.feats, c(3, 5, 8), 'Min'),
    evaluateLassoModels(lassoModels.list, extData.list, shared.feats, c(3, 5, 8), '1SE'))

library(survcomp)
meta.eval.df <- data.frame()
for (sf in c('DSS', 'OS', 'DRFS', 'DFS')) {
    for (t in c(3, 5, 8)) {
        for (mod in unique(ext.conc.df$Model)) {
            for (erstatus in c('All', 'ER+', 'ER-')) {
                single.eval.df <- ext.conc.df %>% filter(Time==paste(t, 'Years'),
                                                         Model==mod,
                                                         `ER Status`==erstatus,
                                                         Outcome==sf)
                single.eval.df <- na.omit(single.eval.df)
                meta.uno <- combine.est(single.eval.df$`Uno Concordance`,
                                        ifelse(single.eval.df$`Uno SE`==0, 0.01,
                                               single.eval.df$`Uno SE`))
                meta.harrell <- combine.est(single.eval.df$`Harrell Concordance`,
                                            ifelse(single.eval.df$`Harrell SE`==0, 0.01,
                                                   single.eval.df$`Harrell SE`))

                meta.eval.df <- rbind(meta.eval.df,
                                      data.frame(
                                          `Uno Concordance` = meta.uno$estimate,
                                          `Harrell Concordance` = meta.harrell$estimate,
                                          `Uno SE` = meta.uno$se,
                                          `Harrell SE` = meta.harrell$se,
                                          Model = mod,
                                          `ER Status` = erstatus,
                                          Time = paste(t, 'Years'),
                                          Outcome = sf,
                                          Dataset = 'Meta Cohort',
                                          `Number of Features`=mean(single.eval.df$`Number of Features`),
                                          check.names=F))
            }
        }
    }
}

ext.conc.df <- rbind(ext.conc.df, meta.eval.df)

ext.conc.df$Dataset <- factor(ext.conc.df$Dataset, c('Meta Cohort', ext.ds.names))
ext.conc.df$Outcome <- factor(ext.conc.df$Outcome, levels=c('DSS', 'OS', 'DRFS', 'DFS'))
ext.conc.df$Time <- factor(ext.conc.df$Time, levels=paste(c(3,5,8), 'Years'))
ext.conc.df$Model <- factor(ext.conc.df$Model)
ext.conc.df$`ER Status` <- factor(ext.conc.df$`ER Status`, levels=c('All', 'ER+', 'ER-'))

write.csv(ext.conc.df, 'ext_concordance_all_numFeats.csv')

library(ggsci)

gg8 <- ggplot(ext.conc.df %>% filter(grepl('Split', Model), !grepl('1SE', Model)),
       aes(x=Time, y=`Uno Concordance`,
           ymin=`Uno Concordance` - 1.96*`Uno SE`,
           ymax=pmin(1,`Uno Concordance` + 1.96*`Uno SE`),
           color=Dataset, lty=Model, shape=Model)) +
    geom_point(position=position_dodge(width=0.8)) +
    geom_errorbar(position='dodge', width=0.8) +
    facet_grid(rows=vars(Outcome), cols=vars(`ER Status`), scales='free') +
    scale_color_d3() +
    scale_linetype_manual(values=c(1, 2), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    scale_shape_manual(values=c(19, 17), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    geom_hline(yintercept=0.5, color='darkgrey', lty=5) +
    ylab("Uno's Concordance") +
    ylim(0.1,1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

gg8

ggsave('external_concordance_all_ds_min.png', gg8, width=12, height=12, dpi=400)

gg9 <- ggplot(ext.conc.df %>%
              filter(grepl('Split', Model), ## !grepl('1SE', Model),
                     Dataset=='Meta Cohort'),
       aes(x=Time, y=`Uno Concordance`,
           ymin=`Uno Concordance` - 1.96*`Uno SE`,
           ymax=pmin(1,`Uno Concordance` + 1.96*`Uno SE`),
           color=`ER Status`, lty=Model, shape=Model)) +
    geom_point(position=position_dodge(width=0.9), size=1.5) +
    geom_errorbar(position=position_dodge(width=0.9), width=0.8) +
    facet_grid(cols=vars(Outcome)) +
    scale_color_brewer(palette='Set1') +
    scale_linetype_manual(values=c(1, 2, 4), labels=c('CausalCoxMGM', 'LASSO Cox 1SE', 'LASSO Cox Min')) +
    scale_shape_manual(values=c(19, 17, 15), labels=c('CausalCoxMGM', 'LASSO Cox 1SE', 'LASSO Cox Min')) +
    geom_hline(yintercept=0.5, color='darkgrey', lty=5) +
    ylab("Uno's Concordance") +
    ## ylim(0.1,1) +
    theme_bw() ## + 
    ## theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

gg9

ggsave('external_concordance_meta_all.png', gg9, width=12, height=4, dpi=400)

gg10 <- ggplot(ext.conc.df %>%
               filter(grepl('Split', Model),
                      Dataset=='Meta Cohort'),
               aes(x=`ER Status`, y=`Number of Features`,
                   fill=Model)) +
    geom_bar(stat='identity', position=position_dodge(width=0.9), width=0.8) +
    facet_grid(cols=vars(Outcome)) +
    scale_fill_brewer(palette='Set2', labels=c('CausalCoxMGM', 'LASSO Cox 1SE', 'LASSO Cox Min')) +
    ## scale_fill_manual(labels=c('CausalCoxMGM', 'LASSO Cox')) +
    ylab("Number of Features") +
    ## ylim(0.1,1) +
    theme_bw() ## + 
    ## theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

gg10

ggsave('final_selected_numFeatures_all.png', gg10, width=12, height=4, dpi=400)

library(ggsci)

gg8 <- ggplot(ext.conc.df %>% filter(grepl('Split', Model), !grepl('Min', Model)),
       aes(x=Time, y=`Uno Concordance`,
           ymin=`Uno Concordance` - 1.96*`Uno SE`,
           ymax=pmin(1,`Uno Concordance` + 1.96*`Uno SE`),
           color=Dataset, lty=Model, shape=Model)) +
    geom_point(position=position_dodge(width=0.8)) +
    geom_errorbar(position='dodge', width=0.8) +
    facet_grid(rows=vars(Outcome), cols=vars(`ER Status`), scales='free') +
    scale_color_d3() +
    scale_linetype_manual(values=c(1, 2), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    scale_shape_manual(values=c(19, 17), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    geom_hline(yintercept=0.5, color='darkgrey', lty=5) +
    ylab("Uno's Concordance") +
    ylim(0.1,1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

gg8

ggsave('external_concordance_all_ds_1se.png', gg8, width=12, height=12, dpi=400)

gg9 <- ggplot(ext.conc.df %>%
              filter(grepl('Split', Model), !grepl('Min', Model),
                     Dataset=='Meta Cohort'),
       aes(x=Time, y=`Uno Concordance`,
           ymin=`Uno Concordance` - 1.96*`Uno SE`,
           ymax=pmin(1,`Uno Concordance` + 1.96*`Uno SE`),
           color=`ER Status`, lty=Model, shape=Model)) +
    geom_point(position=position_dodge(width=0.8)) +
    geom_errorbar(position='dodge', width=0.8) +
    facet_grid(cols=vars(Outcome)) +
    scale_color_brewer(palette='Set1') +
    scale_linetype_manual(values=c(1, 2), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    scale_shape_manual(values=c(19, 17), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    geom_hline(yintercept=0.5, color='darkgrey', lty=5) +
    ylab("Uno's Concordance") +
    ## ylim(0.1,1) +
    theme_bw() ## + 
    ## theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

gg9

ggsave('external_concordance_meta_1se.png', gg9, width=10, height=4, dpi=400)





gg2 <- ggplot(ext.conc.df,
       aes(x=Time, y=`Uno Concordance`,
           ymin=`Uno Concordance` - 1.96*`Uno SE`,
           ymax=pmin(`Uno Concordance` + 1.96*`Uno SE`, 1),
           color=Dataset, lty=Model)) +
    geom_point(position=position_dodge(width=1)) +
    geom_errorbar(position='dodge', width=1) +
    facet_grid(rows=vars(Outcome), cols=vars(`ER Status`)) +
    scale_color_brewer(palette='Paired', direction=-1) +
    ## scale_linetype_manual(values=c(1, 2), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    geom_hline(yintercept=0.5, color='red', lty=5) +
    ylab("Uno's Concordance") +
    ylim(0.2,1.1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

gg2

library(cowplot)

gg <- plot_grid(gg1, gg4, gg5, nrow=1)
gg

ggsave('uno_concordance_ext_both_plot.pdf', gg, height=15, width=15)
ggsave('uno_concordance_ext_both_plot.png', gg, height=15, width=15, dpi=400)


gg <- ggplot(ext.conc.df %>%
             filter(Model %in% c('Split LASSO 1SE', 'Split Causal NB'),
                    Dataset=='Meta Cohort'),
       aes(x=Time, y=`Uno Concordance`,
           ymin=`Uno Concordance` - 1.96*`Uno SE`,
           ymax=`Uno Concordance` + 1.96*`Uno SE`,
           color=`ER Status`, lty=Model)) +
    geom_point(position=position_dodge(width=1)) +
    geom_errorbar(position='dodge', width=1) +
    facet_grid(cols=vars(Outcome)) +
    scale_color_brewer(palette='Set1', direction=1) +
    scale_linetype_manual(values=c(1, 2), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    geom_hline(yintercept=0.5, color='black', lty=5) +
    ylab("Uno's Concordance") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

gg

ggsave('uno_concordance_ext_split_plot_wide.pdf', gg, height=4, width=10)
ggsave('uno_concordance_ext_split_plot_wide.png', gg, height=4, width=10, dpi=400)


gg <- ggplot(ext.conc.df %>%
             filter(Model %in% c('Stratified LASSO 1SE', 'Stratified Causal NB'),
                    Dataset=='Meta Cohort'),
       aes(x=Time, y=`Uno Concordance`,
           ymin=`Uno Concordance` - 1.96*`Uno SE`,
           ymax=`Uno Concordance` + 1.96*`Uno SE`,
           color=`ER Status`, lty=Model)) +
    geom_point(position=position_dodge(width=1)) +
    geom_errorbar(position='dodge', width=1) +
    facet_grid(cols=vars(Outcome)) +
    scale_color_brewer(palette='Set1', direction=1) +
    scale_linetype_manual(values=c(1, 2), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    geom_hline(yintercept=0.5, color='black', lty=5) +
    ylab("Uno's Concordance") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust=NULL, hjust=1))

gg

ggsave('uno_concordance_ext_stratified_plot_wide.pdf', gg, height=4, width=10)
ggsave('uno_concordance_ext_stratified_plot_wide.png', gg, height=4, width=10, dpi=400)



## ext.conc.df$`ER Status` <- factor(ext.conc.df$`ER Status`)

gg <- ggplot(ext.conc.df %>%
             filter(Model %in% c('Split LASSO 1SE', 'Split Causal NB'),
                    Dataset=='Meta Cohort'),
       aes(y=Time, x=`Uno Concordance`,
           xmin=`Uno Concordance` - 1.96*`Uno SE`,
           xmax=`Uno Concordance` + 1.96*`Uno SE`,
           color=`ER Status`, lty=Model)) +
    geom_point(position=position_dodge(width=1)) +
    geom_errorbar(position='dodge', width=1) +
    facet_grid(rows=vars(Outcome)) +
    scale_color_brewer(palette='Set1', direction=1) +
    scale_linetype_manual(values=c(1, 2), labels=c('CausalCoxMGM', 'LASSO Cox')) +
    geom_vline(xintercept=0.5, color='black', lty=5) +
    ylab("Uno's Concordance") +
    theme_bw()## + 
    ## coord_flip()

gg

ggsave('uno_concordance_ext_split_plot_tall.pdf', gg, width=6, height=15)
ggsave('uno_concordance_ext_split_plot_tall.png', gg, width=6, height=15, dpi=400)



hr.forest.df <- data.frame()
mod.names <- c(all='Stratified Full Model', ern='ER- Model', erp='ER+ Model')
for (mod in c('all', 'erp', 'ern')) {
    nb <- sort(getNeighbors(fullData.list[[mod]]$graph[[ifelse(mod=='all', 'StEPS', 'StARS')]], 'DSS'))
    f <- as.formula(paste0('DSS', ' ~ strata(ER_STATUS) + ', paste(nb, collapse=' + ')))
    for (erstatus in c('All', 'ER+', 'ER-')) {
        samp.idx <- 1:nrow(train.npn)
        if (erstatus=='ER+') samp.idx <- samp.idx[train.npn$ER_STATUS=='Positive']
        if (erstatus=='ER-') samp.idx <- samp.idx[train.npn$ER_STATUS=='Negative']
        cph <- coxph(f, train.npn[samp.idx,])
        coefs <- cph$coefficients
        coefs.se <- sqrt(diag(cph$var))

        hr.forest.df <- rbind(hr.forest.df,
                              data.frame(Feature=ifelse(nb=='LYMPH_NODE_STATUS',
                                                        'LNS:Positive', nb),
                                         Coef=coefs,
                                         Coef.SE=coefs.se,
                                         HR=exp(coefs),
                                         HR.Lower=exp(coefs-1.96*coefs.se),
                                         HR.Upper=exp(coefs+1.96*coefs.se),
                                         Model=mod.names[mod],
                                         `ER Status`=erstatus, check.names=F))
    }
}

all.nb.feats <- union(union(nb.all, nb.erp), nb.ern)

hr.forest.df.v2 <- data.frame()

all.nb.feats <- c('LYMPH_NODE_STATUS', 'TUMOR_SIZE', 'UBE2C_grp', 'IGFBP5',
                  'S100P', 'NQO1', 'POSTN', 'FGFR4', 'HLA.DRB4', 'SERPINA1', 'SERPINA6')

f <- as.formula(paste0('DSS', ' ~ strata(ER_STATUS) + ', paste(all.nb.feats, collapse=' + ')))
for (erstatus in c('All', 'ER+', 'ER-')) {
    samp.idx <- 1:nrow(train.npn)
    if (erstatus=='ER+') samp.idx <- samp.idx[train.npn$ER_STATUS=='Positive']
    if (erstatus=='ER-') samp.idx <- samp.idx[train.npn$ER_STATUS=='Negative']
    cph <- coxph(f, train.npn[samp.idx,])
    coefs <- cph$coefficients
    coefs.se <- sqrt(diag(cph$var))

    hr.forest.df.v2 <- rbind(hr.forest.df.v2,
                          data.frame(Feature=ifelse(all.nb.feats=='LYMPH_NODE_STATUS',
                                                    'LNS:Positive', all.nb.feats),
                                     Coef=coefs,
                                     Coef.SE=coefs.se,
                                     HR=exp(coefs),
                                     HR.Lower=exp(coefs-1.96*coefs.se),
                                     HR.Upper=exp(coefs+1.96*coefs.se),
                                     Model=c('All models', 'Stratified & ER+',
                                             'Stratified & ER+', 'Stratified', 'Stratified',
                                             'Stratified', 'ER-', 'ER-', 'ER-',
                                             'ER-', 'ER-'),
                                     `ER Status`=erstatus, check.names=F))
}

hr.forest.df.v2$Feature <- factor(hr.forest.df.v2$Feature,
                                  levels=rev(ifelse(all.nb.feats=='LYMPH_NODE_STATUS',
                                                    'LNS:Positive', all.nb.feats)))

hr.forest.df.v2$`ER Status` <- factor(hr.forest.df.v2$`ER Status`,
                                      levels=c('All', 'ER+', 'ER-'))

hr.forest.df.v2$Model <- factor(hr.forest.df.v2$Model,
                                      levels=c('All models', 'Stratified & ER+',
                                               'Stratified', 'ER-'))

ggplot(hr.forest.df.v2, aes(HR, Feature, xmin=HR.Lower, xmax=HR.Upper,
                            shape=Model,
                            color=`ER Status`)) +
    geom_pointrange(position=position_dodge(width=1), fatten=4) +
    ## facet_grid(cols=vars(Model), rows=vars(`ER Status`)) +
    geom_vline(xintercept=1, col='red', lty=5) +
    scale_x_continuous(trans='log10') +
    scale_color_brewer(palette='Set1', direction=1) +
    xlab('Hazard Ratio') +
    theme_bw()

ggsave('dss_forestplot_ersplit.pdf', width=9, height=9)
ggsave('dss_forestplot_ersplit.png', width=6, height=6, dpi=400)


hr.forest.df <- data.frame()
mod.names <- c(all='Stratified Full Model', ern='ER- Model', erp='ER+ Model')
all.nb.feats <- c()
for (mod in c('all', 'erp', 'ern')) {
    nb <- sort(getNeighbors(fullData.list[[mod]]$graph[[ifelse(mod=='all', 'StEPS', 'StARS')]], 'DistantRelapse'))
    all.nb.feats <- union(all.nb.feats, nb)
    f <- as.formula(paste0('DistantRelapse', ' ~ strata(ER_STATUS) + ', paste(nb, collapse=' + ')))
    for (erstatus in c('All', 'ER+', 'ER-')) {
        samp.idx <- 1:nrow(train.npn)
        if (erstatus=='ER+') samp.idx <- samp.idx[train.npn$ER_STATUS=='Positive']
        if (erstatus=='ER-') samp.idx <- samp.idx[train.npn$ER_STATUS=='Negative']
        cph <- coxph(f, train.npn[samp.idx,])
        coefs <- cph$coefficients
        coefs.se <- sqrt(diag(cph$var))

        hr.forest.df <- rbind(hr.forest.df,
                              data.frame(Feature=ifelse(nb=='LYMPH_NODE_STATUS',
                                                        'LNS:Positive', nb),
                                         Coef=coefs,
                                         Coef.SE=coefs.se,
                                         HR=exp(coefs),
                                         HR.Lower=exp(coefs-1.96*coefs.se),
                                         HR.Upper=exp(coefs+1.96*coefs.se),
                                         Model=mod.names[mod],
                                         `ER Status`=erstatus, check.names=F))
    }
}


## all.nb.feats <- union(union(nb.all, nb.erp), nb.ern)
hr.forest.df$Feature <- factor(hr.forest.df$Feature,
                               levels=rev(ifelse(all.nb.feats=='LYMPH_NODE_STATUS',
                                                 'LNS:Positive', all.nb.feats)))

ggplot(hr.forest.df, aes(HR, Feature, xmin=HR.Lower, xmax=HR.Upper)) +
    geom_pointrange() +
    facet_grid(cols=vars(Model), rows=vars(`ER Status`)) +
    geom_vline(xintercept=1, col='red', lty=5) +
    scale_x_continuous(trans='log10') +
    theme_bw()

ggsave('dr_forestplot_ersplit.pdf', width=9, height=9)
ggsave('dr_forestplot_ersplit.png', width=9, height=9, dpi=400)



library(survminer)
for (sf in c('DSS', 'OS', 'DFS', 'DRFS')) {
    ggforest(model.list[[sf]][['nb.ern']])
    ggsave(paste0(sf,'_ern_forestplot.png'), width=9, height=5, dpi=400)
    ggforest(model.list[[sf]][['nb.erp']])
    ggsave(paste0(sf,'_erp_forestplot.png'), width=9, height=5, dpi=400)
}

for (sf in c('DSS', 'OS', 'DFS', 'DRFS')) {
    if (sf %in% surv.feats[1:4]) {
        nb <- sort(getNeighbors(fullData$graph, sf))
        mb <- sort(setdiff(fullData$graph$markov.blankets[[sf]], surv.feats))
    } else {
        ## nb <- paste0(surv.feats[1:(which(sf==surv.feats)-3)], '.risk.nb')
        ## mb <- paste0(surv.feats[1:(which(sf==surv.feats)-3)], '.risk.mb')
        
        nb <- sort(unique(unlist(lapply(surv.feats[1:(which(sf==surv.feats)-3)],
                                        function(x) getNeighbors(fullData$graph, x)))))
        mb <- setdiff(
            sort(
                unique(
                    unlist(
                        fullData$graph$markov.blankets[surv.feats[1:(which(sf==surv.feats)-3)]]
                    )
                )
            ),
            surv.feats)
    }

    f.nb <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(nb, collapse=' + ')))
    f.mb <- as.formula(paste0(sf, ' ~ strata(ER_STATUS) + ', paste(mb, collapse=' + ')))

    cph.nb <- coxph(f.nb, train.npn)
    cph.mb <- coxph(f.mb, train.npn)

    f.nb <- as.formula(paste0('~ ', paste(nb, collapse=' + ')))
    f.mb <- as.formula(paste0('~ ', paste(mb, collapse=' + ')))
        
    ridge.nb <- glmnet::glmnet(as.matrix(model.matrix(f.nb, train.npn)), train.npn[,sf], family='cox', lambda=0.5, alpha=0)
    ridge.mb <- glmnet::glmnet(as.matrix(model.matrix(f.mb, train.npn)), train.npn[,sf], family='cox', lambda=0.5, alpha=0)

    extData <- read.csv(paste0('../external_validation/',ext.ds,'/',ext.ds, '.csv'))
    extData <- na.omit(extData)

    if (any(c('DSS.time', 'DSS.status') %in% colnames(extData))) {
        extData$DSS <- Surv(extData$DSS.time, extData$DSS.status)
    }
    if (any(c('OS.time', 'OS.status') %in% colnames(extData))) {
        extData$OS <- Surv(extData$OS.time, extData$OS.status)
    }
    if (any(c('DRFS.time', 'DRFS.status') %in% colnames(extData))) {
        extData$DRFS <- Surv(extData$DRFS.time, extData$DRFS.status)
    }
    if (any(c('DFS.time', 'DFS.status') %in% colnames(extData))) {
        extData$DFS <- Surv(extData$DFS.time, extData$DFS.status)
    }

    extData.npn <- extData[,colnames(extData) %in% gsub('_grp', '', colnames(train.npn))] %>%
        mutate_if(function(x) !is.Surv(x) & is.numeric(x),
                  function(x) as.vector(huge::huge.npn(as.matrix(x))))

    colnames(train.npn)[!colnames(train.npn) %in% colnames(extData.npn)]
}
