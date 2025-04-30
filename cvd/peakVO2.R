library(randomForestSRC)
library(rCausalMGM)
library(survival)
library(dplyr)
library(survminer)
library(ggplot2)

data(peakVO2)

peakVO2$survival <- Surv(peakVO2$ttodead, peakVO2$died)
peakVO2 <- peakVO2 %>% select(-c("ttodead","died"))
peakVO2 <- peakVO2 %>% mutate_at(c('betablok', 'dilver', 'nifed', 'acei', 'angioten.II',
                                   'anti.arrhy', 'anti.coag', 'aspirin', 'digoxin',
                                   'nitrates', 'vasodilator', 'diuretic.loop',
                                   'diuretic.thiazide', 'diuretic.potassium.spar',
                                   'lipidrx.statin', 'insulin', 'surgery.pacemaker',
                                   'surgery.cabg', 'surgery.pci', 'surgery.aicd.implant',
                                   'smknow', 'q.wave.mi', 'niddm', 'cad', 'male', 'black'),
                                 factor)

## Stratified 5-fold cross-validation based on all-cause mortality events, select for minimum total deviance of censored outcome

set.seed(37)

idx0 <- which(peakVO2$survival[,2]==0)
idx1 <- which(peakVO2$survival[,2]==1)

foldid <- rep(0, nrow(peakVO2))
foldid[idx0] <- sample(((1:length(idx0))-1) %% 5 + 1)
foldid[idx1] <- sample(((1:length(idx1))-1) %% 5 + 1)

table(foldid, peakVO2$survival[,2])

## lambdas <- runif(60, 0.05, 0.25)
## alphas <- runif(60, 0.001, 0.05)

lambdas <- exp(seq(log(0.25), log(0.05), length.out=20))

alphas <- c(0.01, 0.05, 0.1)

loglik <- matrix(0, 60, 5)
size <- matrix(0, 60, 5)
for (k in 1:5) {
    ig.path <- coxmgmPath(peakVO2[foldid!=k,], lambdas=lambdas, rank=F)
    idx <- 0
    for (ig in ig.path$graphs) {
        for (alph in alphas) {
            idx <- idx + 1
            g <- fciStable(peakVO2[foldid!=k,], initialGraph=ig,
                           alpha=alph, orientRule="maxp", fdr=T, verbose=T)
            mb <- g$markov.blankets$survival
            size[idx,k] <- length(mb)
            if (length(mb)==1) {
                mb <- c(1)
            }
            f <- as.formula(paste("survival ~", paste(mb, collapse=" + ")))
            res <- coxph(f, peakVO2[foldid!=k,])
            test.risk <- predict(res, newdata=peakVO2[foldid==k,])
            res.test <- coxph(survival ~ offset(test.risk), peakVO2[foldid==k,])
            loglik[idx,k] <- -as.numeric(logLik(res.test))
        }
    }
}

sizeMean <- rowMeans(size)
loglikMean <- rowMeans(loglik)
loglikSd <- apply(loglik, 1, sd)

plot(sizeMean, loglikMean, pch=19, col='red')



minIdx <- which.min(loglikMean)

o <- order(sizeMean)

seIdx <- 1
for (idx in 1:100) {
    print(paste("Size =", sizeMean[o[idx]], ": mean =", loglikMean[o[idx]], ": se =", loglikSd[o[idx]]))
    seIdx <- idx
    if (loglikMean[o[idx]] < loglikMean[minIdx] + loglikSd[minIdx]) {
        break
    }
}


## Model selection using BIC for CoxMGM and causal discovery with adjacency false discovery rate (FDR) controled at alpha < 0.1

lambdas <- exp(seq(log(0.25), log(0.05), length.out=50))

ig.path <- coxmgmPath(peakVO2, lambdas=lambdas, verbose=T)

plot(ig.path)

ig.path

plot(ig.path$graph.bic, "survival", list(fontsize=32))

g05 <- fciStable(peakVO2, initialGraph=ig.path$graph.bic,  orientRule="maxp", fdr=T, alpha=0.05, verbose=T)

hist(sapply(g05$neighbors, length))

g05

plot(g05)

plot(g05, "survival", list(fontsize=32))

mb <- g05$markov.blankets$survival

f <- as.formula(paste("survival ~", paste(mb, collapse=" + ")))
res <- coxph(f, peakVO2)

summary(res)


g1 <- fciStable(peakVO2, initialGraph=ig.path$graph.bic,  orientRule="maxp", fdr=T, alpha=0.1, verbose=T)

hist(sapply(g1$neighbors, length))

g1

plot(g1)

saveGraph(g1, "peakVO2_BIC_fdr1.txt")
saveGraph(g1, "peakVO2_BIC_fdr1.sif")

g1 <- loadGraph("peakVO2_BIC_fdr1.txt")

plot(g1, "survival", list(fontsize=32))

plot(g1, "digoxin", list(fontsize=32))

mb <- g1$markov.blankets$survival
f <- as.formula(paste("survival ~", paste(mb, collapse=" + ")))
res <- coxph(f, peakVO2 %>% mutate_at(c("bun", "interval", "lvef.metabl", "peak.vo2"), function(x) (x - mean(x))/sd(x) ))

mb <- g1$markov.blankets$survival[order(coef(res))]
f <- as.formula(paste("survival ~", paste(mb, collapse=" + ")))
res <- coxph(f, peakVO2 %>% mutate_at(c("bun", "interval", "lvef.metabl", "peak.vo2"), function(x) (x - mean(x))/sd(x) ))


summary(res)

ggforest(res)

ggsave("peakVO2_forestplot.png", width=7, height=4)
