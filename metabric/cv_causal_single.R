library(rCausalMGM)
library(survival)
library(dplyr)
library(optparse)

option_list = list(
    make_option(c("-k", "--fold"), type="numeric", default="0",
                help="CV Fold (1 to 10)"),
    make_option(c("-e", "--erstatus"), type="character", default="all",
                help="ER+, ER- or all samples"),
    make_option(c("-m", "--mgm"), action='store_true', default='FALSE',
                help="Learn an initial skeleton with MGM"),
    make_option(c("-f", "--fci"), action='store_true', default='FALSE',
                help="Learn causal graph with FCI50"),
    make_option(c("-r", "--rank"), action='store_true', default='FALSE',
                help="Learn rank-based associations"),
    make_option(c("-s", "--split"), action='store_true', default='FALSE',
                help="Split only or include composite outcomes"),
    make_option(c("-o", "--orientrule"), type="character", default="majority",
                help="Orientation rule: majority, maxp, or conservative"),
    make_option(c("-a", "--alpha"), type="numeric", default="0.1",
                help="Alpha for FDR control")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

k <- opt$fold
erstatus <- opt$erstatus
mgmFlag <- opt$mgm
fciFlag <- opt$fci
rankFlag <- opt$rank
splitFlag <- opt$split
orientRule <- opt$orientrule
alpha <- opt$alpha

train <- read.csv(paste0('data/metabric.rna',ifelse(erstatus=='all','', paste0('.',erstatus)),'.train.cv', k, '.csv'), header=T, row.names=1)

if (any(c('DSS.time', 'DSS.status') %in% colnames(train))) {
   train$DSS <- Surv(train$DSS.time, train$DSS.status)
   if (erstatus=='all' | erstatus=='relapse') {
      attr(train$DSS, 'strata') <- train$ER_STATUS
   }
}
if (any(c('OD.time', 'OD.status') %in% colnames(train))) {
   train$OD <- Surv(train$OD.time, train$OD.status)
   ## if (erstatus=='all') {
   ##    attr(train$OD, 'strata') <- train$ER_STATUS
   ## }
}
if (any(c('DistantRelapse.time', 'DistantRelapse.status') %in% colnames(train))) {
   train$DistantRelapse <- Surv(train$DistantRelapse.time, train$DistantRelapse.status)
   if (erstatus=='all') {
      attr(train$DistantRelapse, 'strata') <- train$ER_STATUS
   }
}
if (any(c('LocalRelapse.time', 'LocalRelapse.status') %in% colnames(train))) {
   train$LocalRelapse <- Surv(train$LocalRelapse.time, train$LocalRelapse.status)
   if (erstatus=='all') {
      attr(train$LocalRelapse, 'strata') <- train$ER_STATUS
   }
}
if (!splitFlag) {
   if (any(c('OS.time', 'OS.status') %in% colnames(train))) {
      train$OS <- Surv(train$OS.time, train$OS.status)
      if (erstatus=='all' | erstatus=='relapse') {
      	 attr(train$OS, 'strata') <- train$ER_STATUS
      }
   }
   if (any(c('DRFS.time', 'DRFS.status') %in% colnames(train))) {
      train$DRFS <- Surv(train$DRFS.time, train$DRFS.status)
      if (erstatus=='all' | erstatus=='relapse') {
      	 attr(train$DRFS, 'strata') <- train$ER_STATUS
      }
   }
   if (any(c('DFS.time', 'DFS.status') %in% colnames(train))) {
      train$DFS <- Surv(train$DFS.time, train$DFS.status)
      if (erstatus=='all' | erstatus=='relapse') {
      	 attr(train$DFS, 'strata') <- train$ER_STATUS
      }
   }
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
   knowledge <- readRDS('knowledge_full.rds')
} else if (erstatus=='relapse') {
   knowledge <- readRDS('knowledge_relapse.rds')
} else {
   knowledge <- readRDS('knowledge_subtype.rds')
}

if (splitFlag) {
   knowledge$tiers[[4]] <- setdiff(knowledge$tiers[[4]], c('OS', 'DFS', 'DRFS'))
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
   train[train$ONCOTREE_CODE!='IDC','ONCOTREE_CODE'] <- 'BREAST'
}

ig <- NULL

splitstatus <- ifelse(splitFlag, 'split', 'composite')

if (mgmFlag) {

    ## ig.steps <- steps(train, rank=rankFlag, verbose=T)

    ## pdf(paste0('steps.instabs.metabric.rna.',splitstatus,'.',erstatus,'.cv', k, '.pdf'), width=5, height=5)
    ## log10params <- log10(ig.steps$lambdas)
    ## plot(x=log10params, y=ig.steps$instability[,6], col='black', pch=19,
    ## 		 xlab=expression(log10(lambda)),
    ##  		 ylab="Edge instability across subsamples",
    ##  		 ylim=c(0,min(0.5, 2*max(ig.steps$instab, na.rm=T))),
    ##  		 cex=1)

    ## points(x=log10params, y=ig.steps$instability[,1], col='red', pch=19, cex=1)
    ## points(x=log10params, y=ig.steps$instability[,2], col='dodgerblue', pch=19, cex=1)
    ## points(x=log10params, y=ig.steps$instability[,3], col='purple', pch=19, cex=1)
    ## points(x=log10params, y=ig.steps$instability[,4], col='orange', pch=19, cex=1)
    ## points(x=log10params, y=ig.steps$instability[,5], col='green', pch=19, cex=1)

    ## abline(h=ig.steps$gamma, lty=5, col='gray', lwd=3, cex=1)

    ## allIdx <- c(1, which(ig.steps$instability[1:which.max(ig.steps$instability[,6]),6]<ig.steps$gamma))
    ## allIdx <- allIdx[length(allIdx)]

    ## abline(v=log10params[allIdx], col='black',  lty=2, lwd=3)
    ## abline(v=log10(ig.steps$graph$lambda[1]), col='red',  lty=2, lwd=3)
    ## abline(v=log10(ig.steps$graph$lambda[2]), col='dodgerblue',  lty=2, lwd=3)
    ## abline(v=log10(ig.steps$graph$lambda[3]), col='purple',  lty=2, lwd=3)
    ## abline(v=log10(ig.steps$graph$lambda[4]), col='orange',  lty=2, lwd=3)
    ## abline(v=log10(ig.steps$graph$lambda[5]), col='green',  lty=2, lwd=3)

    ## legend(x = "topleft", title="Edge Type", 
    ## 	   legend = c("All", "CC", "CD", "DD", "SC", "SD"), 
    ##    	   col = c("black","red", "dodgerblue", "purple", "orange", "green"),
    ##    	   pch = 19, cex=1)

    ## dev.off()

    ## saveGraph(ig.steps$graph,
    ##           paste0('out/mgmStEPS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),
    ##                  '.cv', k, '.txt'))
    ## saveGraph(ig.steps$graph,
    ##           paste0('graph/mgmStEPS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),
    ##                  '.cv', k, '.sif'))

    ## ig.stars <- coxmgm(train, lambda=ig.steps$lambdas[allIdx], rank=rankFlag, verbose=T)

    ## saveGraph(ig.stars,
    ##           paste0('out/mgmStARS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),
    ##                  '.cv', k, '.txt'))
    ## saveGraph(ig.stars,
    ##           paste0('graph/mgmStARS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),
    ##                  '.cv', k, '.sif'))

    ig.path <- coxmgmPath(train, rank=rankFlag, verbose=T)
    
    p <- ncol(train)
    n <- nrow(train)
    kappa <- log(p) / log(n)
    gamma <- max(0, 1-1/(4*kappa)) / 2
    gamma <- 0.25
    penalty <- 1 + 4 * gamma * kappa
    ebic <- -2 *ig.path$loglik + penalty * log(nrow(train)) * ig.path$nParams

    pdf(paste0('solution.path.metabric.rna.',splitstatus,'.',erstatus,'.cv', k, '.pdf'), width=5, height=5)
    plot(ig.path)
    points(log10(ig.path$lambdas), ebic / (2 * nrow(train)), col='purple', pch=19)
    abline(v=log10(ig.path$lambdas)[which.min(ebic)], lty=2, col='purple', lwd=2)
    dev.off()

    ig.path$graph.ebic <- ig.path$graphs[[which.min(ebic)]]

    saveGraph(ig.path$graph.bic,
              paste0('out/mgmBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),
                     '.cv', k, '.txt'))
    saveGraph(ig.path$graph.bic,
              paste0('graph/mgmBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),
                     '.cv', k, '.sif'))

    saveGraph(ig.path$graph.ebic,
              paste0('out/mgmEBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),
                     '.cv', k, '.txt'))
    saveGraph(ig.path$graph.ebic,
              paste0('graph/mgmEBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),
                     '.cv', k, '.sif'))
}

if (fciFlag) {
    ## ig <- loadGraph(paste0('out/mgmStARS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),
    ##                  '.cv', k, '.txt'))

    ## g <- fciStable(train, initialGraph=ig, alpha=alpha, orientRule=orientRule,
    ##                knowledge=knowledge, fdr=T, rank=rankFlag, verbose=T)

    ## saveGraph(g,
    ##           paste0('out/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
    ## 	             'mgmStARS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),'.FDR',
    ## 		     gsub('0[.]', '', as.character(alpha)),
    ##                  '.cv', k, '.txt'))
    
    ## saveGraph(g,
    ##           paste0('graph/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
    ## 	      	     'mgmStARS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),'.FDR',
    ## 		     gsub('0[.]', '', as.character(alpha)),
    ##                  '.cv', k, '.sif'))

    ## ig <- loadGraph(paste0('out/mgmStEPS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),
    ##                  '.cv', k, '.txt'))

    ## g <- fciStable(train, initialGraph=ig, alpha=alpha, orientRule=orientRule,
    ##                knowledge=knowledge, fdr=T, rank=rankFlag, verbose=T)

    ## saveGraph(g,
    ##           paste0('out/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
    ## 	             'mgmStEPS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),'.FDR',
    ## 		     gsub('0[.]', '', as.character(alpha)),
    ##                  '.cv', k, '.txt'))
    
    ## saveGraph(g,
    ##           paste0('graph/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
    ## 	      	     'mgmStEPS.metabric.rna.',splitstatus,'.',erstatus,
    ##                  ifelse(rankFlag, '.rank', '.linear'),'.FDR',
    ## 		     gsub('0[.]', '', as.character(alpha)),
    ##                  '.cv', k, '.sif'))

    ig <- loadGraph(paste0('out/mgmEBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),
                     '.cv', k, '.txt'))

    g <- fciStable(train, initialGraph=ig, alpha=alpha, orientRule=orientRule,
                   knowledge=knowledge, fdr=T, rank=rankFlag, verbose=T)

    saveGraph(g,
              paste0('out/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
	             'mgmEBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),'.FDR',
		     gsub('0[.]', '', as.character(alpha)),
                     '.cv', k, '.txt'))
    
    saveGraph(g,
              paste0('graph/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
	      	     'mgmEBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),'.FDR',
		     gsub('0[.]', '', as.character(alpha)),
                     '.cv', k, '.sif'))

    ig <- loadGraph(paste0('out/mgmBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),
                     '.cv', k, '.txt'))

    g <- fciStable(train, initialGraph=ig, alpha=alpha, orientRule=orientRule,
                   knowledge=knowledge, fdr=T, rank=rankFlag, verbose=T)

    saveGraph(g,
              paste0('out/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
    	             'mgmBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),'.FDR',
    		     gsub('0[.]', '', as.character(alpha)),
                     '.cv', k, '.txt'))
    
    saveGraph(g,
              paste0('graph/', ifelse(orientRule=='majority', 'mfci', 'fcimax'),
    	      	     'mgmBIC.metabric.rna.',splitstatus,'.',erstatus,
                     ifelse(rankFlag, '.rank', '.linear'),'.FDR',
    		     gsub('0[.]', '', as.character(alpha)),
                     '.cv', k, '.sif'))


}