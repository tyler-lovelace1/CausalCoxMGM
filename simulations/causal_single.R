library(rCausalMGM)
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
                help="simulation index", metavar="character"),
    make_option(c("-m", "--mgm"), action='store_true', default='FALSE',
                help="Learn an initial skeleton with MGM"),
    make_option(c("-P", "--pc"), action='store_true', default='FALSE',
                help="Learn causal graph with PC"),
    make_option(c("-f", "--fci"), action='store_true', default='FALSE',
                help="Learn causal graph with FCI"),
    make_option(c("-o", "--orientrule"), type="character", default="majority",
                help="Orientation rule: majority, maxp, or conservative"),
    make_option(c("-a", "--alpha"), type="numeric", default="0.01",
                help="Alpha for CI tests")
    );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

p <- opt$numFeatures
deg <- opt$deg
n <- opt$sampleSize
graphType <- opt$graphType
idx <- opt$index
mgmFlag <- opt$mgm
pcFlag <- opt$pc
fciFlag <- opt$fci
orientRule <- opt$orientrule
## alpha <- opt$alpha

lambdas <- round(10^seq(log10(0.1), log10(0.8), length.out=11), 2)

for (cr in c('0.3', '0.7')) {
    train <- read.csv(paste0('data/sdata_', graphType, '_deg', deg, '_p', p,
                             '_', idx, '_', cr, '_', n, '.csv'), header=T)

    for (i in 1:(p/11)) {
        train[,paste0('Survival', i)] <- Surv(train[,paste0('Survival', i)],
                                              train[,paste0('Censor', i-1)])
    }

    train <- train[,!grepl('Censor', colnames(train))]
    

    if (mgmFlag) {

        ig.steps <- steps(train, lambdas=lambdas, verbose=T)

        pdf(paste0('steps.instabs.', graphType, '_deg', deg, '_p', p, '_', idx,
                   '_', cr, '_', n, '.pdf'), width=5, height=5)
        
        log10params <- log10(ig.steps$lambdas)
        plot(x=log10params, y=ig.steps$instability[,6], col='black', pch=19,
             xlab=expression(log10(lambda)),
             ylab="Edge instability across subsamples",
             ylim=c(0,min(0.5, 2*max(ig.steps$instab, na.rm=T))),
             cex=1)

        points(x=log10params, y=ig.steps$instability[,1], col='red', pch=19, cex=1)
        points(x=log10params, y=ig.steps$instability[,2], col='dodgerblue', pch=19, cex=1)
        points(x=log10params, y=ig.steps$instability[,3], col='purple', pch=19, cex=1)
        points(x=log10params, y=ig.steps$instability[,4], col='orange', pch=19, cex=1)
        points(x=log10params, y=ig.steps$instability[,5], col='green', pch=19, cex=1)
        
        abline(h=ig.steps$gamma, lty=5, col='gray', lwd=3, cex=1)

        allIdx <- c(1, which(ig.steps$instability[1:which.max(ig.steps$instability[,6]),6]<ig.steps$gamma))
        allIdx <- allIdx[length(allIdx)]

        abline(v=log10params[allIdx], col='black',  lty=2, lwd=3)
        abline(v=log10(ig.steps$graph$lambda[1]), col='red',  lty=2, lwd=3)
        abline(v=log10(ig.steps$graph$lambda[2]), col='dodgerblue',  lty=2, lwd=3)
        abline(v=log10(ig.steps$graph$lambda[3]), col='purple',  lty=2, lwd=3)
        abline(v=log10(ig.steps$graph$lambda[4]), col='orange',  lty=2, lwd=3)
        abline(v=log10(ig.steps$graph$lambda[5]), col='green',  lty=2, lwd=3)
        
        legend(x = "topleft", title="Edge Type", 
               legend = c("All", "CC", "CD", "DD", "SC", "SD"), 
               col = c("black","red", "dodgerblue", "purple", "orange", "green"),
               pch = 19, cex=1)

        dev.off()

        ## saveGraph(ig.steps$graph,
        ##           paste0('out/mgmStEPS.metabric.rna.',splitstatus,'.',erstatus,
        ##                  ifelse(rankFlag, '.rank', '.linear'),
        ##                  '.full.txt'))
        ## saveGraph(ig.steps$graph,
        ##           paste0('graph/mgmStEPS.metabric.rna.',splitstatus,'.',erstatus,
        ##                  ifelse(rankFlag, '.rank', '.linear'),
        ##                  '.full.sif'))

        ## ig.steps$graph.stars <- coxmgm(train, lambda=ig.steps$lambdas[allIdx], verbose=T)

        saveRDS(ig.steps,
                paste0('out/coxmgmStEPS_', graphType, '_deg', deg, '_p', p,
                       '_', idx, '_', cr, '_', n, '_v2.rds'))

        ## saveGraph(ig.stars,
        ##           paste0('out/mgmStARS.metabric.rna.',splitstatus,'.',erstatus,
        ##                  ifelse(rankFlag, '.rank', '.linear'),
        ##                  '.full.txt'))
        ## saveGraph(ig.stars,
        ##           paste0('graph/mgmStARS.metabric.rna.',splitstatus,'.',erstatus,
        ##                  ifelse(rankFlag, '.rank', '.linear'),
        ##                  '.full.sif'))

        ig.path <- coxmgmPath(train, lambdas=lambdas, verbose=T)

        kappa <- log(ncol(train)) / log(nrow(train))
        ## gamma <- max(0, 1-1/(4*kappa)) / 2
        gamma <- 0.5
        penalty <- 1 + 4 * gamma * kappa
        ebic <- -2 *ig.path$loglik + penalty * log(nrow(train)) * ig.path$nParams

        pdf(paste0('solution.path.', graphType, '_deg', deg, '_p', p, '_', idx,
                   '_', cr, '_', n, '.full.pdf'), width=5, height=5)
        plot(ig.path)
        points(log10(ig.path$lambdas), ebic / (2 * nrow(train)), col='purple', pch=19)
        abline(v=log10(ig.path$lambdas)[which.min(ebic)], lty=2, col='purple', lwd=2)
        dev.off()

        ig.path$graph.ebic <- ig.path$graphs[[which.min(ebic)]]

        saveRDS(ig.path,
                paste0('out/coxmgmPath_', graphType, '_deg', deg, '_p', p,
                       '_', idx, '_', cr, '_', n, '_v2.rds'))

        ## saveGraph(ig.path$graph.bic,
        ##           paste0('out/mgmBIC.metabric.rna.',splitstatus,'.',erstatus,
        ##                  ifelse(rankFlag, '.rank', '.linear'),
        ##                  '.full.txt'))
        ## saveGraph(ig.path$graph.bic,
        ##           paste0('graph/mgmBIC.metabric.rna.',splitstatus,'.',erstatus,
        ##                  ifelse(rankFlag, '.rank', '.linear'),
        ##                  '.full.sif'))

        ## saveGraph(ig.path$graph.ebic,
        ##           paste0('out/mgmEBIC.metabric.rna.',splitstatus,'.',erstatus,
        ##                  ifelse(rankFlag, '.rank', '.linear'),
        ##                  '.full.txt'))
        ## saveGraph(ig.path$graph.ebic,
        ##           paste0('graph/mgmEBIC.metabric.rna.',splitstatus,'.',erstatus,
        ##                  ifelse(rankFlag, '.rank', '.linear'),
        ##                  '.full.sif'))
    }

    if (fciFlag) {

        for (alpha in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
            ig.steps <- readRDS(paste0('out/coxmgmStEPS_', graphType, '_deg', deg, '_p', p,
                                       '_', idx, '_', cr, '_', n, '.rds'))

            g <- fciStable(train, initialGraph=ig.steps$graph, alpha=alpha,
                           orientRule=orientRule, fdr=F, verbose=T)

            saveGraph(g,
                      paste0('out/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
                             'coxmgmStEPS_', graphType, '_deg', deg, '_p', p,
                             '_', idx, '_', cr, '_', n, '.a',
                             gsub('0[.]', '', as.character(alpha)),
                             '.txt'))

            g <- fciStable(train, initialGraph=ig.steps$graph.stars, alpha=alpha,
                           orientRule=orientRule, fdr=F, verbose=T)

            saveGraph(g,
                      paste0('out/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
                             'coxmgmStARS_', graphType, '_deg', deg, '_p', p,
                             '_', idx, '_', cr, '_', n, '.a',
                             gsub('0[.]', '', as.character(alpha)),
                             '.txt'))

            ig.path <- readRDS(paste0('out/coxmgmPath_', graphType, '_deg', deg, '_p', p,
                                      '_', idx, '_', cr, '_', n, '.rds'))

            g <- fciStable(train, initialGraph=ig.path$graph.bic, alpha=alpha,
                           orientRule=orientRule, fdr=F, verbose=T)

            saveGraph(g,
                      paste0('out/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
                             'coxmgmBIC_', graphType, '_deg', deg, '_p', p,
                             '_', idx, '_', cr, '_', n, '.a',
                             gsub('0[.]', '', as.character(alpha)),
                             '.txt'))
        }
    }

    if (pcFlag) {

        for (alpha in c(0.001, 0.005, 0.01, 0.05, 0.1)) {

            ig.steps <- readRDS(paste0('out/coxmgmStEPS_', graphType, '_deg', deg, '_p', p,
                                       '_', idx, '_', cr, '_', n, '_v2.rds'))

            g <- pcStable(train, initialGraph=ig.steps$graph.steps, alpha=alpha,
                          orientRule=c('majority', 'maxp'), fdr=F, verbose=T)

            saveGraph(g$majority,
                      paste0('out/mpccoxmgmStEPS_', graphType, '_deg', deg, '_p', p,
                             '_', idx, '_', cr, '_', n, '.a',
                             gsub('0[.]', '', as.character(alpha)),
                             '_v2.txt'))

            saveGraph(g$maxp,
                      paste0('out/pcmaxcoxmgmStEPS_', graphType, '_deg', deg, '_p', p,
                             '_', idx, '_', cr, '_', n, '.a',
                             gsub('0[.]', '', as.character(alpha)),
                             '_v2.txt'))


            ## g <- pcStable(train, initialGraph=ig.steps$graph.stars, alpha=alpha,
            ##               orientRule=orientRule, fdr=F, verbose=T)

            ## saveGraph(g,
            ##           paste0('out/', ifelse(orientRule=='majority', 'mpc', 'pcmax'),
            ##                  'coxmgmStARS_', graphType, '_deg', deg, '_p', p,
            ##                  '_', idx, '_', cr, '_', n, '.a',
            ##                  gsub('0[.]', '', as.character(alpha)),
            ##                  '.txt'))

            ## ig.path <- readRDS(paste0('out/coxmgmPath_', graphType, '_deg', deg, '_p', p,
            ##                           '_', idx, '_', cr, '_', n, '.rds'))
            
            ## g <- pcStable(train, initialGraph=ig.path$graph.bic, alpha=alpha,
            ##               orientRule=orientRule, fdr=F, verbose=T)
            
            ## saveGraph(g,
            ##           paste0('out/', ifelse(orientRule=='majority', 'mpc', 'pcmax'),
            ##                  'coxmgmBIC_', graphType, '_deg', deg, '_p', p,
            ##                  '_', idx, '_', cr, '_', n, '.a',
            ##                  gsub('0[.]', '', as.character(alpha)),
            ##                  '.txt'))
        }	
    }
}
