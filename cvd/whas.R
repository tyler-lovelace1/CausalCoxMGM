library(rCausalMGM)
library(survival)
library(dplyr)
library(survminer)
library(ggplot2)

data <- read.csv('whas500.csv', row.names=1)

data$LengthOfStay <- Surv(data$los, 1-data$dstat)

data$Survival <- Surv(data$lenfol, data$fstat)

data <- data %>% select(-c("admitdate", "disdate", "fdate", "dstat", "lenfol", "fstat", "year", "los"))

data <- data %>% mutate_at(c("cvd", "afb", "sho", "chf", "av3", "gender", "miord", "mitype"), factor)

## Stratified 5-fold cross-validation based on all-cause mortality and hospital discharge events, select for minimum total deviance of censored outcomes

set.seed(43)

idx00 <- which(data$Survival[,2]==0 & data$LengthOfStay[,2]==0)
idx10 <- which(data$Survival[,2]==1 & data$LengthOfStay[,2]==0)
idx01 <- which(data$Survival[,2]==0 & data$LengthOfStay[,2]==1)
idx11 <- which(data$Survival[,2]==1 & data$LengthOfStay[,2]==1)

foldid <- rep(0, 500)
foldid[idx00] <- sample(((1:length(idx00))-1) %% 5 + 1)
foldid[idx10] <- sample(((1:length(idx10))-1) %% 5 + 1)
foldid[idx01] <- sample(((1:length(idx01))-1) %% 5 + 1)
foldid[idx11] <- sample(((1:length(idx11))-1) %% 5 + 1)

table(foldid, data$Survival[,2])
table(foldid, data$LengthOfStay[,2])

lambdas <- runif(100, 0.05, 0.5)

alphas <- runif(100, 0.01, 0.25)

loglik <- matrix(0, 100, 5)
size <- matrix(0, 100, 5)
for (k in 1:5) {
    ig.path <- coxmgmPath(data[foldid!=k,], lambdas=lambdas, rank=F)
    idx <- 0
    for (ig in ig.path$graphs) {
        idx <- idx + 1
        g <- fciStable(data[foldid!=k,], initialGraph=ig,
                       alpha=alphas[idx], orientRule="maxp", rank=F)
        mb <- g$markov.blankets$Survival
        size[idx,k] <- length(mb)
        if (length(mb)==1) {
            mb <- c(1)
        }
        f <- as.formula(paste("Survival ~", paste(mb, collapse=" + ")))
        res <- coxph(f, data[foldid!=k,])
        test.risk <- predict(res, newdata=data[foldid==k,])
        res.test <- coxph(Survival ~ offset(test.risk), data[foldid==k,])
        loglik[idx,k] <- -as.numeric(logLik(res.test))

        mb <- g$markov.blankets$LengthOfStay
        size[idx,k] <- size[idx,k] + length(mb)
        if (length(mb)==1) {
            mb <- c(1)
        }
        f <- as.formula(paste("LengthOfStay ~", paste(mb, collapse=" + ")))
        res <- coxph(f, data[foldid!=k,])
        test.risk <- predict(res, newdata=data[foldid==k,])
        res.test <- coxph(LengthOfStay ~ offset(test.risk), data[foldid==k,])
        loglik[idx,k] <- loglik[idx,k] + -as.numeric(logLik(res.test))
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


ig.path <- coxmgmPath(data, lambda=lambdas)

g <- fciStable(data, initialGraph=ig.path$graphs[[minIdx]],
               alpha=alphas[minIdx], verbose=T, orientRule="maxp")

g

plot(g)

mb <- g$markov.blankets$Survival
f <- as.formula(paste("Survival ~", paste(mb, collapse=" + ")))
res <- coxph(f, data %>% mutate_at(c("diasbp", "age", "hr", "bmi"), function(x) (x - mean(x))/sd(x)) %>% mutate_if(is.factor, as.numeric))

summary(res)

mb <- g$markov.blankets$Survival[order(coef(res))]
f <- as.formula(paste("Survival ~", paste(mb, collapse=" + ")))
res <- coxph(f, data %>% mutate_at(c("diasbp", "age", "hr", "bmi"), function(x) (x - mean(x))/sd(x)) %>% mutate_if(is.factor, as.numeric))


summary(res)

ggforest(res, data %>% mutate_at(c("diasbp", "age", "hr", "bmi"), function(x) (x - mean(x))/sd(x)) %>% mutate_if(is.factor, as.numeric))

ggsave("whas500_survival_forestplot.png", width=7, height=4)


mb <- g$markov.blankets$LengthOfStay
f <- as.formula(paste("LengthOfStay ~", paste(mb, collapse=" + ")))
res <- coxph(f, data)

mb <- g$markov.blankets$LengthOfStay[order(coef(res))]
f <- as.formula(paste("LengthOfStay ~", paste(mb, collapse=" + ")))
res <- coxph(f, data)


summary(res)

ggforest(res, data=data)

ggsave("whas500_lengthofstay_forestplot.png", width=7, height=4)


write.csv(loglik, "whas500_loglik.csv")
write.csv(size, "whas500_mbSize.csv")

saveGraph(g, "whas500_randCV.txt")
saveGraph(g, "whas500_randCV.sif")

g <- loadGraph("whas500_randCV.txt")


g.boot <- bootstrap(data, graph=g, numBoots=200, verbose=T)

plot(g.boot)

graphdf <- graphTable(g, g.boot$stabilities)
graphdf$Edge <- g$edges

graphdf <- graphdf[order(graphdf$adjFreq, graphdf$orientFreq),]
graphdf$Edge <- factor(graphdf$Edge, levels=graphdf$Edge)

ggplot(graphdf, aes(x=adjFreq, y=Edge)) +
    geom_bar(stat='identity') +
    geom_bar(aes(orientFreq, Edge), graphdf, stat='identity', fill='red') +
    theme_bw()

