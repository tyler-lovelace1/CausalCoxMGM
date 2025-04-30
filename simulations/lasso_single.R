library(rCausalMGM)
library(glmnet)
library(survival)
library(dplyr)
library(optparse)

option_list = list(
    make_option(c("-d", "--deg"), type="character", default="4",
                help="average graph degree", metavar="character"),
    make_option(c("-p", "--numFeatures"), type="numeric", default="110",
                help="number of variables", metavar="numeric"),
    make_option(c("-n", "--sampleSize"), type="character", default="00500",
                help="number of samples", metavar="character"),
    make_option(c("-t", "--graphType"), type="character", default="ER",
                help="graph type", metavar="character"),
    make_option(c("-i", "--index"), type="character", default="0",
                help="simulation index", metavar="character")
    );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

p <- opt$numFeatures
deg <- opt$deg
n <- opt$sampleSize
graphType <- opt$graphType
idx <- opt$index

lambdas <- round(10^seq(log10(0.1), log10(0.8), length.out=11), 2)

for (cr in c('0.3', '0.7')) {
    train <- read.csv(paste0('data/sdata_', graphType, '_deg', deg, '_p', p,
                             '_', idx, '_', cr, '_', n, '.csv'), header=T)

    for (i in 1:(p/11)) {
        train[,paste0('Survival', i)] <- Surv(train[,paste0('Survival', i)],
                                              train[,paste0('Censor', i-1)])
    }

    train <- train[,!grepl('Censor', colnames(train))]

    xmat <- model.matrix(~.,train[,!grepl('Survival', colnames(train))])[,-1]

    ## group.vec <- rep(0, ncol(xmat))
    ## group <- 1
    ## for (col in colnames(train)[!grepl('Survival', colnames(train))]) {
    ##     group.idx <- grep(paste0('^',col,'$'), colnames(xmat))
    ##     if (length(group.idx)==0) {
    ##         group.idx <- grep(paste0('^',col,'Category'), colnames(xmat))
    ##     }
    ##     group.vec[group.idx] <- group
    ##     group <- group + 1
    ## }

    ## sgl.results <- list()
    glmnet.results <- list()
    for (i in 1:(p/11)) {
        glmnet.results[[paste0('Survival', i)]] <- cv.glmnet(x=xmat,
                                                             y=train[,paste0('Survival', i)],
                                                             family='cox',
                                                             lambda=lambdas)
        
        ## sgl.results[[paste0('Survival', i)]] <- cvSGL(list(x=xmat,
        ##                                                    time=train[,paste0('Survival', i)][,1],
        ##                                                    status=train[,paste0('Survival', i)][,2]),
        ##                                               index=group.vec, type='cox', alpha=0,
        ##                                               nlam=30, min.frac=0.01)
    }

    ## saveRDS(sgl.results, paste0('out/coxSGL_', graphType, '_deg', deg, '_p', p,
    ##                             '_', idx, '_', cr, '_', n, '.rds'))

    saveRDS(glmnet.results, paste0('out/coxLASSO_', graphType, '_deg', deg, '_p', p,
                                   '_', idx, '_', cr, '_', n, '.rds'))
}
