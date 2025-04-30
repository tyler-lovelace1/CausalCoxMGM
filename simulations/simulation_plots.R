library(rCausalMGM)
library(stringr)
library(ggplot2)
library(pracma)
library(parallel)
library(dplyr)
library(cowplot)

subGraphEdgeType <- function(graph, edgetype='all') {
    require(stringr)
    edgetype <- toupper(edgetype)
    if (edgetype=='ALL') {
        return(graph)
    }
    cont.counts <- sapply(graph$edges, str_count, pattern='X')
    disc.counts <- sapply(graph$edges, str_count, pattern='Y')
    surv.counts <- sapply(graph$edges, str_count, pattern='Survival')
    ## cont.nodes <- graph$nodes[grepl('X', graph$nodes)]
    ## cont.nodes <- graph$nodes[grepl('Y', graph$nodes)]
    ## cont.nodes <- graph$nodes[grepl('X', graph$nodes)]
    subGraph <- graph

    if (edgetype=='CC') {
        subGraph$edges <- subGraph$edges[cont.counts==2]
    } else if (edgetype=='CD') {
        subGraph$edges <- subGraph$edges[cont.counts==1 & disc.counts==1]
    } else if (edgetype=='DD') {
        subGraph$edges <- subGraph$edges[disc.counts==2]
    } else if (edgetype=='SC') {
        subGraph$edges <- subGraph$edges[cont.counts==1 & surv.counts==1]
    } else if (edgetype=='SD') {
        subGraph$edges <- subGraph$edges[disc.counts==1 & surv.counts==1]
    } else if (edgetype=='SURV') {
        subGraph$edges <- subGraph$edges[surv.counts==1]
    } else {
        stop(paste("Unrecognized edge type:", edgetype))
    }
    return(subGraph)
}

prMetricsByEdgeTypeFast <- function(est, true) {
    edgetypes <- c('ALL', 'CC', 'CD', 'DD', 'SC', 'SD', 'SURV', 'NoSURV')
    tpVals <- rep(0, 8)
    fpVals <- rep(0, 8)
    fnVals <- rep(0, 8)
    estEdgeList <- sort(est$edges)
    trueEdgeList <- sort(skeleton(true)$edges)

    idx <- 1
    jdx <- 1
    
    while (idx <= length(estEdgeList) || jdx <= length(trueEdgeList)) {
        edge1 <- estEdgeList[idx]
        esplit1 <- strsplit(edge1, ' ')[[1]]
        
        survCount1 <- sum(grepl("Survival", esplit1))
        contCount1 <- sum(grepl("X", esplit1))
        discCount1 <- sum(grepl("Y", esplit1))

        edge2 <- trueEdgeList[jdx]
        esplit2 <- strsplit(edge2, ' ')[[1]]

        if (!is.na(edge1==edge2) && edge1==edge2) {
            kdx <- c(1)
            if (survCount1==0) {
                kdx <- c(kdx, 8)
                if (contCount1==2) {
                    kdx <- c(kdx, 2)
                } else {
                    if (discCount1==2) {
                        kdx <- c(kdx, 4)
                    } else {
                        kdx <- c(kdx, 3)
                    }
                }
            } else {
                kdx <- c(kdx, 7)
                if (contCount1==1) {
                    kdx <- c(kdx, 5)
                } else {
                    kdx <- c(kdx, 6)
                }
            }
            tpVals[kdx] <- tpVals[kdx] + 1
            idx <- idx + 1
            jdx <- jdx + 1
        } else {
            if (is.na(edge1)) {
                survCount2 <- sum(grepl("Survival", esplit2))
                contCount2 <- sum(grepl("X", esplit2))
                discCount2 <- sum(grepl("Y", esplit2))
                
                kdx <- c(1)
                if (survCount2==0) {
                    kdx <- c(kdx, 8)
                    if (contCount2==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount2==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount2==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fnVals[kdx] <- fnVals[kdx] + 1
                jdx <- jdx + 1
            } else if (is.na(edge2) || edge1 < edge2) {
                kdx <- c(1)
                if (survCount1==0) {
                    kdx <- c(kdx, 8)
                    if (contCount1==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount1==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount1==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fpVals[kdx] <- fpVals[kdx] + 1
                idx <- idx + 1
            } else {
                survCount2 <- sum(grepl("Survival", esplit2))
                contCount2 <- sum(grepl("X", esplit2))
                discCount2 <- sum(grepl("Y", esplit2))
                
                kdx <- c(1)
                if (survCount2==0) {
                    kdx <- c(kdx, 8)
                    if (contCount2==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount2==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount2==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fnVals[kdx] <- fnVals[kdx] + 1
                jdx <- jdx + 1
            }
        }
    }

    skelPrec <- tpVals / (tpVals + fpVals)
    skelRec <- tpVals / (tpVals + fnVals)

    tpVals <- rep(0, 8)
    fpVals <- rep(0, 8)
    fnVals <- rep(0, 8)
    estEdgeList <- sort(est$edges)
    trueEdgeList <- sort(moral(true)$edges)

    idx <- 1
    jdx <- 1
    
    while (idx <= length(estEdgeList) || jdx <= length(trueEdgeList)) {
        edge1 <- estEdgeList[idx]
        esplit1 <- strsplit(edge1, ' ')[[1]]
        
        survCount1 <- sum(grepl("Survival", esplit1))
        contCount1 <- sum(grepl("X", esplit1))
        discCount1 <- sum(grepl("Y", esplit1))

        edge2 <- trueEdgeList[jdx]
        esplit2 <- strsplit(edge2, ' ')[[1]]

        if (!is.na(edge1==edge2) && edge1==edge2) {
            kdx <- c(1)
            if (survCount1==0) {
                kdx <- c(kdx, 8)
                if (contCount1==2) {
                    kdx <- c(kdx, 2)
                } else {
                    if (discCount1==2) {
                        kdx <- c(kdx, 4)
                    } else {
                        kdx <- c(kdx, 3)
                    }
                }
            } else {
                kdx <- c(kdx, 7)
                if (contCount1==1) {
                    kdx <- c(kdx, 5)
                } else {
                    kdx <- c(kdx, 6)
                }
            }
            tpVals[kdx] <- tpVals[kdx] + 1
            idx <- idx + 1
            jdx <- jdx + 1
        } else {
            if (is.na(edge1)) {
                survCount2 <- sum(grepl("Survival", esplit2))
                contCount2 <- sum(grepl("X", esplit2))
                discCount2 <- sum(grepl("Y", esplit2))
                
                kdx <- c(1)
                if (survCount2==0) {
                    kdx <- c(kdx, 8)
                    if (contCount2==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount2==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount2==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fnVals[kdx] <- fnVals[kdx] + 1
                jdx <- jdx + 1
            } else if (is.na(edge2) || edge1 < edge2) {
                kdx <- c(1)
                if (survCount1==0) {
                    kdx <- c(kdx, 8)
                    if (contCount1==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount1==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount1==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fpVals[kdx] <- fpVals[kdx] + 1
                idx <- idx + 1
            } else {
                survCount2 <- sum(grepl("Survival", esplit2))
                contCount2 <- sum(grepl("X", esplit2))
                discCount2 <- sum(grepl("Y", esplit2))
                
                kdx <- c(1)
                if (survCount2==0) {
                    kdx <- c(kdx, 8)
                    if (contCount2==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount2==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount2==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fnVals[kdx] <- fnVals[kdx] + 1
                jdx <- jdx + 1
            }
        }
    }

    moralPrec <- tpVals / (tpVals + fpVals)
    moralRec <- tpVals / (tpVals + fnVals)

    data.frame(`Edge Type`=edgetypes,
               Lambda = c(mean(est$lambda), est$lambda,
                          mean(est$lambda[4:5]),
                          mean(est$lambda[1:3])),
               skeletonPrecision = skelPrec,
               skeletonRecall = skelRec,
               skeletonF1 = 2 * (skelPrec * skelRec) / (skelPrec + skelRec),
               moralPrecision = moralPrec,
               moralRecall = moralRec,
               moralF1 = 2 * (moralPrec * moralRec) / (moralPrec + moralRec),
               check.names=F)
}

allMetricsByEdgeTypeFast <- function(est, true) {
    edgetypes <- c('ALL', 'CC', 'CD', 'DD', 'SC', 'SD', 'SURV', 'NoSURV')
    tpVals <- rep(0, 8)
    fpVals <- rep(0, 8)
    fnVals <- rep(0, 8)
    tpValsOrient <- rep(0, 8)
    fpValsOrient <- rep(0, 8)
    fnValsOrient <- rep(0, 8)
    tpValsCausal <- rep(0, 8)
    fpValsCausal <- rep(0, 8)
    fnValsCausal <- rep(0, 8)
    shdVals <- rep(0, 8)

    estAdjList <- skeleton(est)$edges
    estEdgeList <- est$edges
    trueAdjList <- skeleton(true)$edges
    tureEdgeList <- true$edges
    estAdjO <- order(estAdjList)
    trueAdjO <- order(trueAdjList)

    idx <- 1
    jdx <- 1
    
    while (idx <= length(estEdgeList) || jdx <= length(trueEdgeList)) {
        edge1 <- estAdjList[estAdjO[idx]]
        esplit1 <- strsplit(edge1, ' ')[[1]]
        
        survCount1 <- sum(grepl("Survival", esplit1))
        contCount1 <- sum(grepl("X", esplit1))
        discCount1 <- sum(grepl("Y", esplit1))

        edge2 <- trueEdgeList[jdx]
        esplit2 <- strsplit(edge2, ' ')[[1]]

        if (!is.na(edge1==edge2) && edge1==edge2) {
            kdx <- c(1)
            if (survCount1==0) {
                kdx <- c(kdx, 8)
                if (contCount1==2) {
                    kdx <- c(kdx, 2)
                } else {
                    if (discCount1==2) {
                        kdx <- c(kdx, 4)
                    } else {
                        kdx <- c(kdx, 3)
                    }
                }
            } else {
                kdx <- c(kdx, 7)
                if (contCount1==1) {
                    kdx <- c(kdx, 5)
                } else {
                    kdx <- c(kdx, 6)
                }
            }
            tpVals[kdx] <- tpVals[kdx] + 1
            idx <- idx + 1
            jdx <- jdx + 1
        } else {
            if (is.na(edge1)) {
                survCount2 <- sum(grepl("Survival", esplit2))
                contCount2 <- sum(grepl("X", esplit2))
                discCount2 <- sum(grepl("Y", esplit2))
                
                kdx <- c(1)
                if (survCount2==0) {
                    kdx <- c(kdx, 8)
                    if (contCount2==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount2==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount2==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fnVals[kdx] <- fnVals[kdx] + 1
                jdx <- jdx + 1
            } else if (is.na(edge2) || edge1 < edge2) {
                kdx <- c(1)
                if (survCount1==0) {
                    kdx <- c(kdx, 8)
                    if (contCount1==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount1==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount1==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fpVals[kdx] <- fpVals[kdx] + 1
                idx <- idx + 1
            } else {
                survCount2 <- sum(grepl("Survival", esplit2))
                contCount2 <- sum(grepl("X", esplit2))
                discCount2 <- sum(grepl("Y", esplit2))
                
                kdx <- c(1)
                if (survCount2==0) {
                    kdx <- c(kdx, 8)
                    if (contCount2==2) {
                        kdx <- c(kdx, 2)
                    } else {
                        if (discCount2==2) {
                            kdx <- c(kdx, 4)
                        } else {
                            kdx <- c(kdx, 3)
                        }
                    }
                } else {
                    kdx <- c(kdx, 7)
                    if (contCount2==1) {
                        kdx <- c(kdx, 5)
                    } else {
                        kdx <- c(kdx, 6)
                    }
                }
                fnVals[kdx] <- fnVals[kdx] + 1
                jdx <- jdx + 1
            }
        }
    }

    skelPrec <- tpVals / (tpVals + fpVals)
    skelRec <- tpVals / (tpVals + fnVals)
}

prMetricsByEdgeType <- function(est, true) {
    edgetypes <- c('ALL', 'CC', 'CD', 'DD', 'SC', 'SD', 'SURV')
    graphtypes <- c('Skeleton', 'Moral')
    pr.metrics <- data.frame()
    gtype.pr.metrics <- list()
    for (gtype in graphtypes) {
        gtype.pr.metrics[[gtype]] <- data.frame()
        if (gtype=='Skeleton') {
            true.undir <- skeleton(true)
        } else if (gtype=='Moral') {
            true.undir <- moral(true)
        } else {
            stop("Invalid graph type")
        }

        for (etype in edgetypes) {
            pr.vec <- prMetricsAdjacency(subGraphEdgeType(est,etype),
                                         subGraphEdgeType(true.undir,etype))
            lam.idx <- 1:5
            if (etype=='SURV') {
                lam.idx <- 4:5
            } else if (etype %in% c('CC', 'CD', 'DD', 'SC', 'SD')) {
                lam.idx <- which(etype == edgetypes)-1
            }

            if (gtype=='Skeleton') {
                gtype.pr.metrics[[gtype]] <- rbind(
                    gtype.pr.metrics[[gtype]],
                    data.frame(`Edge Type`=etype,
                               Lambda=mean(est$lambda[lam.idx]),
                               skeletonPrecision=pr.vec[1],
                               skeletonRecall=pr.vec[2],
                               skeletonF1=pr.vec[3],
                               check.names=F))
            }

            if (gtype=='Moral') {
                gtype.pr.metrics[[gtype]] <- rbind(
                    gtype.pr.metrics[[gtype]],
                    data.frame(`Edge Type`=etype,
                               Lambda=mean(est$lambda[lam.idx]),
                               moralPrecision=pr.vec[1],
                               moralRecall=pr.vec[2],
                               moralF1=pr.vec[3],
                               check.names=F))
            }
        }
    }
    pr.metrics <- data.frame(gtype.pr.metrics[[1]],
                             gtype.pr.metrics[[2]][,3:5],
                             check.names=F)
    
    pr.metrics[is.na(pr.metrics[,colnames(pr.metrics)[grepl('Precision', colnames(pr.metrics))][1]]),grepl('Precision', colnames(pr.metrics))] <- 1
    pr.metrics[is.na(pr.metrics[,colnames(pr.metrics)[grepl('F1', colnames(pr.metrics))][1]]),grepl('F1', colnames(pr.metrics))] <- 0
    rownames(pr.metrics) <- NULL
    return(pr.metrics)
}

directionalPrMetric <- function(est, true) {
    pos.est <- length(est$edges)
    pos.true <- length(true$edges)
    tp <- sum(est$edges %in% true$edges)
    prec <- tp / pos.est
    rec <- tp / pos.true
    f1 <- 2 * prec * rec / (prec + rec)
    return(c(directionalPrecision=prec,
             directionalRecall=rec,
             directionalF1=f1))
}

allMetricsByEdgeType <- function(est, true) {
    edgetypes <- c('ALL', 'CC', 'CD', 'DD', 'SC', 'SD', 'SURV')
    ## graphtypes <- c('Skeleton', 'Moral')
    if (true$type=="") {
        true$type <- "directed acyclic graph"
    }
    true.cpdag <- cpdag(true)
    all.metrics <- data.frame()
    for (etype in edgetypes) {
        all.vec <- allMetrics(subGraphEdgeType(est,etype),
                              subGraphEdgeType(true.cpdag,etype))
        dir.vec <- prMetricsCausal(subGraphEdgeType(est,etype),
                                   subGraphEdgeType(true,etype))
        all.metrics <- rbind(all.metrics,
                             data.frame(`Edge Type`=etype,
                                        Lambda=ifelse(is.null(est$lambda[1]),
                                                      NA, est$lambda[1]),
                                        Alpha=est$alpha,
                                        SHD=all.vec['SHD'],
                                        skeletonPrecision=all.vec['adjPrecision'],
                                        skeletonRecall=all.vec['adjRecall'],
                                        skeletonF1=all.vec['adjF1'],
                                        skeletonMCC=all.vec['adjMCC'],
                                        orientationPrecision=all.vec['orientPrecision'],
                                        orientationRecall=all.vec['orientRecall'],
                                        orientationF1=all.vec['orientF1'],
                                        orientationMCC=all.vec['orientMCC'],
                                        causalPrecision=dir.vec['causalPrecision'],
                                        causalRecall=dir.vec['causalRecall'],
                                        causalF1=dir.vec['causalF1'],
                                        check.names=F))
    }
    all.metrics[is.na(all.metrics[,'skeletonPrecision']),'skeletonPrecision'] <- 1
    all.metrics[is.na(all.metrics[,'orientationPrecision']),'orientationPrecision'] <- 1
    all.metrics[is.na(all.metrics[,'causalPrecision']),'causalPrecision'] <- 1
    all.metrics[is.na(all.metrics[,'skeletonF1']),'skeletonF1'] <- 0
    all.metrics[is.na(all.metrics[,'orientationF1']),'orientationF1'] <- 0
    all.metrics[is.na(all.metrics[,'causalF1']),'causalF1'] <- 0
    all.metrics[is.na(all.metrics[,'skeletonRecall']),'skeletonRecall'] <- 0
    all.metrics[is.na(all.metrics[,'orientationRecall']),'orientationRecall'] <- 0
    all.metrics[is.na(all.metrics[,'causalRecall']),'causalRecall'] <- 0
    
    rownames(all.metrics) <- NULL
    return(all.metrics)
}

calcAuPR <- function(recall, precision) {
    mask <- !(is.na(recall) | is.na(precision))
    recall <- c(0, recall[mask], 1)
    precision <- c(1, precision[mask], 0)
    return(trapz(recall, precision))
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


lassoFeatSelectEval <- function(true, lassoRes) {
    nodes <- true$nodes
    
    group.vec <- c()
    group <- 1
    for (col in nodes[!grepl('Survival', nodes)]) {
        if (grepl('X', col)) group.vec <- c(group.vec, group)
        else if (grepl('Y', col)) group.vec <- c(group.vec, rep(group,2))
        group <- group + 1
    }

    p.true <- 0
    tp.min <- 0
    tp.1se <- 0
    predp.min <- 0
    predp.1se <- 0
    for (col in nodes[grepl('Survival', nodes)]) {
        nb.true <- getNeighbors(true, col)
        nb.min <- nodes[unique(group.vec[as.vector(coef(lassoRes[[col]], s=lassoRes[[col]]$lambda.min))!=0])]
        nb.1se <- nodes[unique(group.vec[as.vector(coef(lassoRes[[col]], s=lassoRes[[col]]$lambda.1se))!=0])]

        ## print(col)
        ## print(sort(nb.true))
        ## print(sort(nb.1se))
        ## print(sort(nb.min))
        ## print("")
        
        p.true <- p.true + length(nb.true)

        tp.min <- tp.min + length(intersect(nb.true, nb.min))
        tp.1se <- tp.1se + length(intersect(nb.true, nb.1se))

        predp.min <- predp.min + length(nb.min)
        predp.1se <- predp.1se + length(nb.1se)
    
    }

    prec.min <- tp.min/predp.min
    prec.1se <- tp.1se/predp.1se

    rec.min <- tp.min/p.true
    rec.1se <- tp.1se/p.true

    f1.min <- 2 * prec.min * rec.min / (prec.min + rec.min)
    f1.1se <- 2 * prec.1se * rec.1se / (prec.1se + rec.1se)

    eval.df <- data.frame(Method='LASSO',
                          `Feature Select`=c('Min', '1SE'),
                          Alpha=NA,
                          Precision=c(prec.min, prec.1se),
                          Recall=c(rec.min, rec.1se),
                          F1=c(f1.min, f1.1se),
                          check.names=F)

    return(eval.df)

}

causalFeatSelectEval <- function(true, est) {
    nodes <- true$nodes
    
    p.true <- 0
    tp.mb <- 0
    tp.nb <- 0
    predp.mb <- 0
    predp.nb <- 0
    for (col in nodes[grepl('Survival', nodes)]) {
        nb.true <- getNeighbors(true, col)
        nb <- getNeighbors(est, col)
        mb <- est$markov.blankets[[col]]

        ## print(col)
        ## print(sort(nb.true))
        ## print(sort(nb))
        ## print(sort(mb))
        ## print("")
        
        p.true <- p.true + length(nb.true)

        tp.mb <- tp.mb + length(intersect(nb.true, mb))
        tp.nb <- tp.nb + length(intersect(nb.true, nb))

        predp.mb <- predp.mb + length(mb)
        predp.nb <- predp.nb + length(nb)
    
    }

    prec.mb <- tp.mb/predp.mb
    prec.nb <- tp.nb/predp.nb

    rec.mb <- tp.mb/p.true
    rec.nb <- tp.nb/p.true

    f1.mb <- 2 * prec.mb * rec.mb / (prec.mb + rec.mb)
    f1.nb <- 2 * prec.nb * rec.nb / (prec.nb + rec.nb)

    eval.df <- data.frame(Method='CausalCoxMGM',
                          `Feature Select`=c('MB', 'NB'),
                          Alpha=est$alpha,
                          Precision=c(prec.mb, prec.nb),
                          Recall=c(rec.mb, rec.nb),
                          F1=c(f1.mb, f1.nb),
                          check.names=F)

    return(eval.df)

}

edgetypes <- c('ALL', 'CC', 'CD', 'DD', 'SC', 'SD', 'SURV')

graphTypes <- c('SF', 'ER')
sampSizes <- c('00100', '00250', '00500', '01000', '05000', '10000')
censorRates <- c('0.3', '0.7')
numFeats <- c(55, 110, 550)
degrees <- c(2, 4, 6)
lambdas <- 10^seq(log10(0.05), log10(0.8), length.out=20)

#### CoxMGM precision recall curves across sample sizes

p <- 110
d <- 4

mgmBySampSizeMetrics <- data.frame()
mgmBySampSizeSelectMetrics <- data.frame()
mgmBySampSizeSummaryMetrics <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (N in sampSizes) {
            for (idx in (0:19)) {
                ig.path <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                }, error= function(cond) { NULL })
                
                ig.steps <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmStEPS_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                }, error= function(cond) { NULL })

                if (is.null(ig.path)) {
                    next
                }

                if (is.null(ig.steps)) {
                    next
                }
                
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))
                singleMgmMetrics <- do.call(rbind.data.frame, mclapply(ig.path$graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmMetrics$F1 <- 2 * (singleMgmMetrics$moralPrecision * singleMgmMetrics$skeletonRecall) / (singleMgmMetrics$moralPrecision + singleMgmMetrics$skeletonRecall)
                singleMgmMetrics$`Graph Type` <- gt
                singleMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmMetrics$`Sample Size` <- as.numeric(N)
                singleMgmMetrics$Degree <- d
                singleMgmMetrics$`Number of Features` <- p
                singleMgmMetrics$`Index` <- idx
                ## mgmBySampSizeMetrics <- rbind(mgmBySampSizeMetrics, singleMgmMetrics)
                mgmBySampSizeMetrics[nrow(mgmBySampSizeMetrics)+1:nrow(singleMgmMetrics),colnames(singleMgmMetrics)] <- singleMgmMetrics
                
                singleMgmSummaryMetrics <- data.frame()
                for (edgetype in edgetypes) {
                    singleMgmSummaryMetrics <- rbind(
                        singleMgmSummaryMetrics,
                        data.frame(`Edge Type`=edgetype,
                                   AUPR=calcAuPR(
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'moralPrecision']),
                                   check.names=F))
                }
                singleMgmSummaryMetrics$`Graph Type` <- gt
                singleMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSummaryMetrics$Degree <- d
                singleMgmSummaryMetrics$`Number of Features` <- p
                singleMgmSummaryMetrics$`Index` <- idx
                ## mgmBySampSizeSummaryMetrics <- rbind(mgmBySampSizeSummaryMetrics,
                ##                                      singleMgmSummaryMetrics)

                mgmBySampSizeSummaryMetrics[nrow(mgmBySampSizeSummaryMetrics)+1:nrow(singleMgmSummaryMetrics),colnames(singleMgmSummaryMetrics)] <- singleMgmSummaryMetrics

                mgmSelect.graphs <- list(BIC=ig.path$graph.bic, StARS=ig.steps$graph.stars, StEPS=ig.steps$graph, Oracle=ig.path$graphs[[which.max(singleMgmMetrics$F1[singleMgmMetrics$`Edge Type`=='ALL'])]])
                    
                singleMgmSelectMetrics <- do.call(rbind.data.frame, mclapply(mgmSelect.graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmSelectMetrics$F1 <- 2 * (singleMgmSelectMetrics$moralPrecision * singleMgmSelectMetrics$skeletonRecall) / (singleMgmSelectMetrics$moralPrecision + singleMgmSelectMetrics$skeletonRecall)
                singleMgmSelectMetrics$`Graph Type` <- gt
                singleMgmSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSelectMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSelectMetrics$Degree <- d
                singleMgmSelectMetrics$`Number of Features` <- p
                singleMgmSelectMetrics$`Index` <- idx
                singleMgmSelectMetrics$`CoxMGM Select` <- gsub('[.][0-9]$', '', rownames(singleMgmSelectMetrics))
                ## mgmBySampSizeSelectMetrics <- rbind(mgmBySampSizeSelectMetrics, singleMgmSelectMetrics)
                mgmBySampSizeSelectMetrics[nrow(mgmBySampSizeSelectMetrics)+1:nrow(singleMgmSelectMetrics),colnames(singleMgmSelectMetrics)] <- singleMgmSelectMetrics
                

            }
        }
    }
}

#### CoxMGM precision recall curves across number of features

N <- '00500'
d <- 4

mgmByNumFeatsMetrics <- data.frame()
mgmByNumFeatsSelectMetrics <- data.frame()
mgmByNumFeatsSummaryMetrics <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (p in numFeats) {
            for (idx in (0:19)) {
                if (p == 550 && idx >= 10) next
                
                ig.path <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                }, error= function(cond) { NULL })
                
                ig.steps <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmStEPS_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                }, error= function(cond) { NULL })

                if (is.null(ig.path)) {
                    next
                }

                if (is.null(ig.steps)) {
                    next
                }
                
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))
                singleMgmMetrics <- do.call(rbind.data.frame, mclapply(ig.path$graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmMetrics$F1 <- 2 * (singleMgmMetrics$moralPrecision * singleMgmMetrics$skeletonRecall) / (singleMgmMetrics$moralPrecision + singleMgmMetrics$skeletonRecall)
                singleMgmMetrics$`Graph Type` <- gt
                singleMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmMetrics$`Sample Size` <- as.numeric(N)
                singleMgmMetrics$Degree <- d
                singleMgmMetrics$`Number of Features` <- p
                singleMgmMetrics$`Index` <- idx
                ## mgmByNumFeatsMetrics <- rbind(mgmByNumFeatsMetrics, singleMgmMetrics)
                mgmByNumFeatsMetrics[nrow(mgmByNumFeatsMetrics)+1:nrow(singleMgmMetrics),colnames(singleMgmMetrics)] <- singleMgmMetrics
                
                singleMgmSummaryMetrics <- data.frame()
                for (edgetype in edgetypes) {
                    singleMgmSummaryMetrics <- rbind(
                        singleMgmSummaryMetrics,
                        data.frame(`Edge Type`=edgetype,
                                   AUPR=calcAuPR(
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'moralPrecision']),
                                   check.names=F))
                }
                singleMgmSummaryMetrics$`Graph Type` <- gt
                singleMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSummaryMetrics$Degree <- d
                singleMgmSummaryMetrics$`Number of Features` <- p
                singleMgmSummaryMetrics$`Index` <- idx
                ## mgmByNumFeatsSummaryMetrics <- rbind(mgmByNumFeatsSummaryMetrics,
                ##                                      singleMgmSummaryMetrics)

                mgmByNumFeatsSummaryMetrics[nrow(mgmByNumFeatsSummaryMetrics)+1:nrow(singleMgmSummaryMetrics),colnames(singleMgmSummaryMetrics)] <- singleMgmSummaryMetrics

                mgmSelect.graphs <- list(BIC=ig.path$graph.bic, StARS=ig.steps$graph.stars, StEPS=ig.steps$graph, Oracle=ig.path$graphs[[which.max(singleMgmMetrics$F1[singleMgmMetrics$`Edge Type`=='ALL'])]])
                    
                singleMgmSelectMetrics <- do.call(rbind.data.frame, mclapply(mgmSelect.graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmSelectMetrics$F1 <- 2 * (singleMgmSelectMetrics$moralPrecision * singleMgmSelectMetrics$skeletonRecall) / (singleMgmSelectMetrics$moralPrecision + singleMgmSelectMetrics$skeletonRecall)
                singleMgmSelectMetrics$`Graph Type` <- gt
                singleMgmSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSelectMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSelectMetrics$Degree <- d
                singleMgmSelectMetrics$`Number of Features` <- p
                singleMgmSelectMetrics$`Index` <- idx
                singleMgmSelectMetrics$`CoxMGM Select` <- gsub('[.][0-9]$', '', rownames(singleMgmSelectMetrics))
                ## mgmByNumFeatsSelectMetrics <- rbind(mgmByNumFeatsSelectMetrics, singleMgmSelectMetrics)
                mgmByNumFeatsSelectMetrics[nrow(mgmByNumFeatsSelectMetrics)+1:nrow(singleMgmSelectMetrics),colnames(singleMgmSelectMetrics)] <- singleMgmSelectMetrics
                

            }
        }
    }
}


#### CoxMGM precision recall curves across average degree

N <- '00500'
p <- 110

mgmByDegreeMetrics <- data.frame()
mgmByDegreeSelectMetrics <- data.frame()
mgmByDegreeSummaryMetrics <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (d in degrees) {
            for (idx in (0:19)) {
                if (p == 550 && idx >= 10) next
                
                ig.path <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                }, error= function(cond) { NULL })
                
                ig.steps <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmStEPS_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                }, error= function(cond) { NULL })

                if (is.null(ig.path)) {
                    next
                }

                if (is.null(ig.steps)) {
                    next
                }
                
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))
                singleMgmMetrics <- do.call(rbind.data.frame, mclapply(ig.path$graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmMetrics$F1 <- 2 * (singleMgmMetrics$moralPrecision * singleMgmMetrics$skeletonRecall) / (singleMgmMetrics$moralPrecision + singleMgmMetrics$skeletonRecall)
                singleMgmMetrics$`Graph Type` <- gt
                singleMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmMetrics$`Sample Size` <- as.numeric(N)
                singleMgmMetrics$Degree <- d
                singleMgmMetrics$`Number of Features` <- p
                singleMgmMetrics$`Index` <- idx
                ## mgmByDegreeMetrics <- rbind(mgmByDegreeMetrics, singleMgmMetrics)
                mgmByDegreeMetrics[nrow(mgmByDegreeMetrics)+1:nrow(singleMgmMetrics),colnames(singleMgmMetrics)] <- singleMgmMetrics
                
                singleMgmSummaryMetrics <- data.frame()
                for (edgetype in edgetypes) {
                    singleMgmSummaryMetrics <- rbind(
                        singleMgmSummaryMetrics,
                        data.frame(`Edge Type`=edgetype,
                                   AUPR=calcAuPR(
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'moralPrecision']),
                                   check.names=F))
                }
                singleMgmSummaryMetrics$`Graph Type` <- gt
                singleMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSummaryMetrics$Degree <- d
                singleMgmSummaryMetrics$`Number of Features` <- p
                singleMgmSummaryMetrics$`Index` <- idx
                ## mgmByDegreeSummaryMetrics <- rbind(mgmByDegreeSummaryMetrics,
                ##                                      singleMgmSummaryMetrics)

                mgmByDegreeSummaryMetrics[nrow(mgmByDegreeSummaryMetrics)+1:nrow(singleMgmSummaryMetrics),colnames(singleMgmSummaryMetrics)] <- singleMgmSummaryMetrics

                mgmSelect.graphs <- list(BIC=ig.path$graph.bic, StARS=ig.steps$graph.stars, StEPS=ig.steps$graph, Oracle=ig.path$graphs[[which.max(singleMgmMetrics$F1[singleMgmMetrics$`Edge Type`=='ALL'])]])
                    
                singleMgmSelectMetrics <- do.call(rbind.data.frame, mclapply(mgmSelect.graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmSelectMetrics$F1 <- 2 * (singleMgmSelectMetrics$moralPrecision * singleMgmSelectMetrics$skeletonRecall) / (singleMgmSelectMetrics$moralPrecision + singleMgmSelectMetrics$skeletonRecall)
                singleMgmSelectMetrics$`Graph Type` <- gt
                singleMgmSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSelectMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSelectMetrics$Degree <- d
                singleMgmSelectMetrics$`Number of Features` <- p
                singleMgmSelectMetrics$`Index` <- idx
                singleMgmSelectMetrics$`CoxMGM Select` <- gsub('[.][0-9]$', '', rownames(singleMgmSelectMetrics))
                ## mgmByDegreeSelectMetrics <- rbind(mgmByDegreeSelectMetrics, singleMgmSelectMetrics)
                mgmByDegreeSelectMetrics[nrow(mgmByDegreeSelectMetrics)+1:nrow(singleMgmSelectMetrics),colnames(singleMgmSelectMetrics)] <- singleMgmSelectMetrics
                

            }
        }
    }
}


#### CoxMGM precision recall curves across number of features

N <- '00500'
d <- 4

mgmByNumFeatsMetricsV2 <- data.frame()
mgmByNumFeatsSelectMetricsV2 <- data.frame()
mgmByNumFeatsSummaryMetricsV2 <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (p in numFeats) {
            for (idx in (0:19)) {
                if (p == 550 && idx >= 10) next
                
                ig.path <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '_v2.rds'))
                }, error= function(cond) { NULL })
                
                ig.steps <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmStEPS_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '_v2.rds'))
                }, error= function(cond) { NULL })

                if (is.null(ig.path)) {
                    next
                }

                if (is.null(ig.steps)) {
                    next
                }
                
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))
                singleMgmMetrics <- do.call(rbind.data.frame, mclapply(ig.path$graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmMetrics$F1 <- 2 * (singleMgmMetrics$moralPrecision * singleMgmMetrics$skeletonRecall) / (singleMgmMetrics$moralPrecision + singleMgmMetrics$skeletonRecall)
                singleMgmMetrics$`Graph Type` <- gt
                singleMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmMetrics$`Sample Size` <- as.numeric(N)
                singleMgmMetrics$Degree <- d
                singleMgmMetrics$`Number of Features` <- p
                singleMgmMetrics$`Index` <- idx
                ## mgmByNumFeatsMetricsV2 <- rbind(mgmByNumFeatsMetricsV2, singleMgmMetrics)
                mgmByNumFeatsMetricsV2[nrow(mgmByNumFeatsMetricsV2)+1:nrow(singleMgmMetrics),colnames(singleMgmMetrics)] <- singleMgmMetrics
                
                singleMgmSummaryMetrics <- data.frame()
                for (edgetype in edgetypes) {
                    singleMgmSummaryMetrics <- rbind(
                        singleMgmSummaryMetrics,
                        data.frame(`Edge Type`=edgetype,
                                   AUPR=calcAuPR(
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'moralPrecision']),
                                   check.names=F))
                }
                singleMgmSummaryMetrics$`Graph Type` <- gt
                singleMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSummaryMetrics$Degree <- d
                singleMgmSummaryMetrics$`Number of Features` <- p
                singleMgmSummaryMetrics$`Index` <- idx
                ## mgmByNumFeatsSummaryMetricsV2 <- rbind(mgmByNumFeatsSummaryMetricsV2,
                ##                                      singleMgmSummaryMetrics)

                mgmByNumFeatsSummaryMetricsV2[nrow(mgmByNumFeatsSummaryMetricsV2)+1:nrow(singleMgmSummaryMetrics),colnames(singleMgmSummaryMetrics)] <- singleMgmSummaryMetrics

                mgmSelect.graphs <- list(BIC=ig.path$graph.bic, StARS=ig.steps$graph.stars, StEPS=ig.steps$graph.steps, Oracle=ig.path$graphs[[which.max(singleMgmMetrics$F1[singleMgmMetrics$`Edge Type`=='ALL'])]])
                    
                singleMgmSelectMetrics <- do.call(rbind.data.frame, mclapply(mgmSelect.graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmSelectMetrics$F1 <- 2 * (singleMgmSelectMetrics$moralPrecision * singleMgmSelectMetrics$skeletonRecall) / (singleMgmSelectMetrics$moralPrecision + singleMgmSelectMetrics$skeletonRecall)
                singleMgmSelectMetrics$`Graph Type` <- gt
                singleMgmSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSelectMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSelectMetrics$Degree <- d
                singleMgmSelectMetrics$`Number of Features` <- p
                singleMgmSelectMetrics$`Index` <- idx
                singleMgmSelectMetrics$`CoxMGM Select` <- gsub('[.][0-9]$', '', rownames(singleMgmSelectMetrics))
                ## mgmByNumFeatsSelectMetricsV2 <- rbind(mgmByNumFeatsSelectMetricsV2, singleMgmSelectMetrics)
                mgmByNumFeatsSelectMetricsV2[nrow(mgmByNumFeatsSelectMetricsV2)+1:nrow(singleMgmSelectMetrics),colnames(singleMgmSelectMetrics)] <- singleMgmSelectMetrics
                

            }
        }
    }
}


#### CoxMGM precision recall curves across average degree

N <- '00500'
p <- 110

mgmByDegreeMetricsV2 <- data.frame()
mgmByDegreeSelectMetricsV2 <- data.frame()
mgmByDegreeSummaryMetricsV2 <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (d in degrees) {
            for (idx in (0:19)) {
                if (p == 550 && idx >= 10) next
                
                ig.path <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '_v2.rds'))
                }, error= function(cond) { NULL })
                
                ig.steps <- tryCatch(expr={
                    readRDS(paste0('out/coxmgmStEPS_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '_v2.rds'))
                }, error= function(cond) { NULL })

                if (is.null(ig.path)) {
                    next
                }

                if (is.null(ig.steps)) {
                    next
                }
                
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))
                singleMgmMetrics <- do.call(rbind.data.frame, mclapply(ig.path$graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmMetrics$F1 <- 2 * (singleMgmMetrics$moralPrecision * singleMgmMetrics$skeletonRecall) / (singleMgmMetrics$moralPrecision + singleMgmMetrics$skeletonRecall)
                singleMgmMetrics$`Graph Type` <- gt
                singleMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmMetrics$`Sample Size` <- as.numeric(N)
                singleMgmMetrics$Degree <- d
                singleMgmMetrics$`Number of Features` <- p
                singleMgmMetrics$`Index` <- idx
                ## mgmByDegreeMetricsV2 <- rbind(mgmByDegreeMetricsV2, singleMgmMetrics)
                mgmByDegreeMetricsV2[nrow(mgmByDegreeMetricsV2)+1:nrow(singleMgmMetrics),colnames(singleMgmMetrics)] <- singleMgmMetrics
                
                singleMgmSummaryMetrics <- data.frame()
                for (edgetype in edgetypes) {
                    singleMgmSummaryMetrics <- rbind(
                        singleMgmSummaryMetrics,
                        data.frame(`Edge Type`=edgetype,
                                   AUPR=calcAuPR(
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'moralPrecision']),
                                   check.names=F))
                }
                singleMgmSummaryMetrics$`Graph Type` <- gt
                singleMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSummaryMetrics$Degree <- d
                singleMgmSummaryMetrics$`Number of Features` <- p
                singleMgmSummaryMetrics$`Index` <- idx
                ## mgmByDegreeSummaryMetricsV2 <- rbind(mgmByDegreeSummaryMetricsV2,
                ##                                      singleMgmSummaryMetrics)

                mgmByDegreeSummaryMetricsV2[nrow(mgmByDegreeSummaryMetricsV2)+1:nrow(singleMgmSummaryMetrics),colnames(singleMgmSummaryMetrics)] <- singleMgmSummaryMetrics

                mgmSelect.graphs <- list(BIC=ig.path$graph.bic, StARS=ig.steps$graph.stars, StEPS=ig.steps$graph.steps, Oracle=ig.path$graphs[[which.max(singleMgmMetrics$F1[singleMgmMetrics$`Edge Type`=='ALL'])]])
                    
                singleMgmSelectMetrics <- do.call(rbind.data.frame, mclapply(mgmSelect.graphs, prMetricsByEdgeTypeFast, true=true, mc.cores=16))
                singleMgmSelectMetrics$F1 <- 2 * (singleMgmSelectMetrics$moralPrecision * singleMgmSelectMetrics$skeletonRecall) / (singleMgmSelectMetrics$moralPrecision + singleMgmSelectMetrics$skeletonRecall)
                singleMgmSelectMetrics$`Graph Type` <- gt
                singleMgmSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSelectMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSelectMetrics$Degree <- d
                singleMgmSelectMetrics$`Number of Features` <- p
                singleMgmSelectMetrics$`Index` <- idx
                singleMgmSelectMetrics$`CoxMGM Select` <- gsub('[.][0-9]$', '', rownames(singleMgmSelectMetrics))
                ## mgmByDegreeSelectMetricsV2 <- rbind(mgmByDegreeSelectMetricsV2, singleMgmSelectMetrics)
                mgmByDegreeSelectMetricsV2[nrow(mgmByDegreeSelectMetricsV2)+1:nrow(singleMgmSelectMetrics),colnames(singleMgmSelectMetrics)] <- singleMgmSelectMetrics
                

            }
        }
    }
}



## mgmBySampSizeMetrics$F1 <- 2 * (mgmBySampSizeMetrics$moralPrecision * mgmBySampSizeMetrics$skeletonRecall) / (mgmBySampSizeMetrics$moralPrecision + mgmBySampSizeMetrics$skeletonRecall)

## mgmBySampSizeSelectMetrics$F1 <- 2 * (mgmBySampSizeSelectMetrics$moralPrecision * mgmBySampSizeSelectMetrics$skeletonRecall) / (mgmBySampSizeSelectMetrics$moralPrecision + mgmBySampSizeSelectMetrics$skeletonRecall)

library(dplyr)

ggplot(mgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
       group_by(Lambda, `Graph Type`, `Edge Type`, `Censoring Rate`, `Sample Size`) %>%
       summarize_at(c('skeletonRecall', 'moralPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=moralPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=moralPrecision_mean - 1.96 * moralPrecision_se,
           ymax=moralPrecision_mean + 1.96 * moralPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Precision') +
    xlab('Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('coxmgm_pr_curves_by_samplesize.png', width=10, height=6, dpi=400)

ggplot(mgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500, 1000)) %>%
       group_by(Lambda, `Graph Type`, `Edge Type`, `Censoring Rate`, `Sample Size`) %>%
       summarize_at(c('skeletonRecall', 'moralPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=moralPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=moralPrecision_mean - 1.96 * moralPrecision_se,
           ymax=moralPrecision_mean + 1.96 * moralPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Precision') +
    xlab('Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('coxmgm_pr_curves_reduced.png', width=7, height=6, dpi=400)


p1 <- ggplot(mgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(mgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(mgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.25,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

library(cowplot)
plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('coxmgm_aupr_curves_reduced.png', width=6, height=4, dpi=400)


p1 <- ggplot(mgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))) %>%
       mutate(Comparison='Sample Size'),
       aes(x=`Sample Size`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background.y = element_blank(),
          strip.text.y = element_blank())


p3 <- ggplot(mgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))) %>%
       mutate(Comparison='Degree'),
       aes(x=Degree, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(mgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))) %>%
       mutate(Comparison='Number of Features'),
       aes(x=`Number of Features`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.25,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background.y = element_blank(),
          strip.text.y = element_blank(), axis.ticks.y = element_blank())

library(cowplot)
plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('coxmgm_aupr_curves_all.png', width=6, height=3, dpi=400)


p1 <- ggplot(mgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))) %>%
       mutate(Comparison='Sample Size'),
       aes(x=`Sample Size`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(Comparison)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background.y = element_blank(),
          strip.text.y = element_blank())


p3 <- ggplot(mgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x))) %>%
       mutate(Comparison='Degree')),
       aes(x=Degree, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(Comparison)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(mgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('AUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x))) %>%
       mutate(Comparison='Number of Features')),
       aes(x=`Number of Features`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(Comparison)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.25,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background.y = element_blank(),
          strip.text.y = element_blank(), axis.ticks.y = element_blank())

library(cowplot)
plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('coxmgm_aupr_curves_surv.png', width=6, height=4, dpi=400)


ggplot(mgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
       group_by(Lambda, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at(c('skeletonRecall', 'moralPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=moralPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=moralPrecision_mean - 1.96 * moralPrecision_se,
           ymax=moralPrecision_mean + 1.96 * moralPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=4) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Adjacency Precision') +
    xlab('Adjacency Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    theme_bw()



ggplot(mgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')),
       aes(`Censoring Rate`, AUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=1), width=0.5) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('AUPRC') +
    xlab("% Censored") +
    theme_bw()

ggsave('coxmgm_AUPRC_by_samplesize.png', width=7, height=6, dpi=400)


ggplot(mgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500, 1000)),
       aes(`Censoring Rate`, AUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=1), width=0.5) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('AUPRC') +
    xlab("% Censored") +
    theme_bw()

ggsave('coxmgm_AUPRC_reduced.png', width=5, height=6, dpi=400)


ggplot(mgmBySampSizeMetrics %>%
       group_by(Lambda, `Graph Type`, `Edge Type`, `Censoring Rate`, `Sample Size`) %>%
       summarize_at(c('skeletonRecall', 'moralPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=moralPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=moralPrecision_mean - 1.96 * moralPrecision_se,
           ymax=moralPrecision_mean + 1.96 * moralPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Precision') +
    xlab('Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    theme_bw()

ggsave('coxmgm_pr_curves_by_samplesize_all_edgetypes.png', width=15, height=15, dpi=400)

ggplot(mgmBySampSizeSummaryMetrics,
       aes(`Censoring Rate`, AUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position='dodge') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('AUPRC') +
    theme_bw()

ggsave('coxmgm_AUPRC_by_samplesize_all_edgetypes.png', width=15, height=15, dpi=400)


## mgmBySampSizeSelectMetrics$F1 <- 2 * (mgmBySampSizeSelectMetrics$moralPrecision * mgmBySampSizeSelectMetrics$skeletonRecall) / (mgmBySampSizeSelectMetrics$moralPrecision + mgmBySampSizeSelectMetrics$skeletonRecall)

ggplot(mgmBySampSizeSelectMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       mutate_at('CoxMGM Select', factor, levels=c('Oracle', 'BIC', 'StARS', 'StEPS')),
       aes(`Sample Size`, F1, color=`Graph Type`,
           shape=`Censoring Rate`, lty=`Censoring Rate`)) +
    ## geom_bar(stat='summary', fun.data='mean_se', position='dodge', width=0.9) +
    geom_line(stat='summary') +
    geom_pointrange(stat='summary', fun.data='mean_se',
                    fun.args=list(mult=1.96), fatten=2) +
    ## geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
    ##               position=position_dodge(width=0.9), width=0.8) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`CoxMGM Select`)) +
    ylab('F1') +
    scale_x_continuous(trans='log10') + 
    ## scale_shape_manual(values=15:18) +
    labs(color='Graph Type', shape="% Censored", lty="% Censored") +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('coxmgm_F1_by_samplesize_reduced.png', width=6, height=3, dpi=400)


ggplot(mgmByNumFeatsSelectMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       mutate_at('CoxMGM Select', factor, levels=c('Oracle', 'BIC', 'StARS', 'StEPS')),
       aes(`Number of Features`, F1, color=`Graph Type`,
           shape=`Censoring Rate`, lty=`Censoring Rate`)) +
    ## geom_bar(stat='summary', fun.data='mean_se', position='dodge', width=0.9) +
    geom_line(stat='summary') +
    geom_pointrange(stat='summary', fun.data='mean_se',
                    fun.args=list(mult=1.96), fatten=2) +
    ## geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
    ##               position=position_dodge(width=0.9), width=0.8) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`CoxMGM Select`)) +
    ylab('F1') +
    scale_x_continuous(trans='log10') + 
    ## scale_shape_manual(values=15:18) +
    labs(color='Graph Type', shape="% Censored", lty="% Censored") +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('coxmgm_F1_by_numFeats_reduced.png', width=6, height=3, dpi=400)


ggplot(mgmByDegreeSelectMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       mutate_at('CoxMGM Select', factor, levels=c('Oracle', 'BIC', 'StARS', 'StEPS')),
       aes(Degree, F1, color=`Graph Type`,
           shape=`Censoring Rate`, lty=`Censoring Rate`)) +
    ## geom_bar(stat='summary', fun.data='mean_se', position='dodge', width=0.9) +
    geom_line(stat='summary') +
    geom_pointrange(stat='summary', fun.data='mean_se',
                    fun.args=list(mult=1.96), fatten=2) +
    ## geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
    ##               position=position_dodge(width=0.9), width=0.8) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`CoxMGM Select`)) +
    ylab('F1') +
    ## scale_x_continuous(trans='log10') + 
    ## scale_shape_manual(values=15:18) +
    labs(color='Graph Type', shape="% Censored", lty="% Censored") +
    theme_bw(base_size=7)

ggsave('coxmgm_F1_by_degree_reduced.png', width=6, height=3, dpi=400)



ggplot(mgmBySampSizeSelectMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
       mutate_at('CoxMGM Select', factor, levels=c('Oracle', 'BIC', 'StARS', 'StEPS')),
       aes(`CoxMGM Select`, F1, color=`Graph Type`, shape=`Censoring Rate`)) +
    ## geom_bar(stat='summary', fun.data='mean_se', position='dodge', width=0.9) +
    ## geom_line() +
    geom_pointrange(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                    position=position_dodge(width=0.9), fatten=3) +
    ## geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
    ##               position=position_dodge(width=0.9), width=0.8) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`), scales='free_y') +
    ylab('F1') +
    ## scale_shape_manual(values=15:18) +
    labs(color='Graph', shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('coxmgm_F1_by_samplesize_reduced.png', width=10, height=6, dpi=400)



ggplot(mgmBySampSizeSelectMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
       mutate_at('CoxMGM Select', factor, levels=c('Oracle', 'BIC', 'StARS', 'StEPS')),
       aes(`Censoring Rate`, F1, fill=`Graph Type`, shape=`CoxMGM Select`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_pointrange(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=1)) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('F1') +
    labs(shape="Select") +
    theme_bw()

ggsave('coxmgm_F1_by_samplesize.png', width=15, height=10, dpi=400)


ggplot(mgmBySampSizeSelectMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500)) %>%
       mutate_at('CoxMGM Select', factor, levels=c('Oracle', 'BIC', 'StARS', 'StEPS')),
       aes(`Censoring Rate`, F1, fill=`Graph Type`, pattern=`CoxMGM Select`)) +
    geom_bar_pattern(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position='dodge') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('F1') +
    theme_bw()

ggsave('coxmgm_F1_reduced.png', width=15, height=10, dpi=400)

write.csv(mgmBySampSizeMetrics, 'sim_results/coxmgm_by_samplesize.csv')
write.csv(mgmBySampSizeSummaryMetrics, 'sim_results/coxmgm_by_samplesize_summary.csv')
write.csv(mgmBySampSizeSelectMetrics, 'sim_results/coxmgm_by_samplesize_select.csv')

mgmBySampSizeMetrics <- read.csv('sim_results/coxmgm_by_samplesize.csv',
                                 row.names=1, check.names=F)
mgmBySampSizeSummaryMetrics <- read.csv('sim_results/coxmgm_by_samplesize_summary.csv',
                                        row.names=1, check.names=F)
mgmBySampSizeSelectMetrics <- read.csv('sim_results/coxmgm_by_samplesize_select.csv',
                                       row.names=1, check.names=F)


#### CausalCoxMGM precision recall curves (adjacency & orientation) across sample size

p <- 110
d <- 4
alphas <- c(0.001, 0.005, 0.01, 0.05, 0.1)

causalmgmBySampSizeMetrics <- data.frame()
causalmgmBySampSizeSummaryMetrics <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (N in sampSizes) {
            for (idx in 0:19) {
                ## ig.path <- readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))

                for (select in c('StEPS')) {
                    causal.graphs <- list()
                    for (alpha in alphas) {
                        causal.graphs[[which(alphas==alpha)]] <- tryCatch(expr={
                            loadGraph(
                                paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                                       '_', idx, '_', cr, '_', N, '.a',
                                       gsub('0[.]', '', as.character(alpha)),
                                       '.txt'))
                        }, error= function(cond) { NULL })
                    }
                    if (length(causal.graphs) > 0) {
                        singleCausalMgmMetrics <- do.call(rbind.data.frame, mclapply(causal.graphs, allMetricsByEdgeType, true=true, mc.cores=5))
                        singleCausalMgmMetrics$`Graph Type` <- gt
                        singleCausalMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                        singleCausalMgmMetrics$`Sample Size` <- as.numeric(N)
                        singleCausalMgmMetrics$Degree <- d
                        singleCausalMgmMetrics$`Number of Features` <- p
                        singleCausalMgmMetrics$`CoxMGM Select` <- select
                        singleCausalMgmMetrics$`Index` <- idx
                        causalmgmBySampSizeMetrics[nrow(causalmgmBySampSizeMetrics)+1:nrow(singleCausalMgmMetrics),colnames(singleCausalMgmMetrics)] <- singleCausalMgmMetrics
                        
                        singleCausalMgmSummaryMetrics <- data.frame()
                        for (edgetype in edgetypes[-8]) {
                            singleCausalMgmSummaryMetrics <- rbind(
                                singleCausalMgmSummaryMetrics,
                                data.frame(`Edge Type`=edgetype,
                                           skeletonAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'skeletonPrecision']),
                                           orientationAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'orientationRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'orientationPrecision']),
                                           causalAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'causalRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'causalPrecision']),
                                           check.names=F))
                        }
                        singleCausalMgmSummaryMetrics$`Graph Type` <- gt
                        singleCausalMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                        singleCausalMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                        singleCausalMgmSummaryMetrics$Degree <- d
                        singleCausalMgmSummaryMetrics$`Number of Features` <- p
                        singleCausalMgmSummaryMetrics$`CoxMGM Select` <- select
                        singleCausalMgmSummaryMetrics$`Index` <- idx
                        causalmgmBySampSizeSummaryMetrics[nrow(causalmgmBySampSizeSummaryMetrics)+1:7,colnames(singleCausalMgmSummaryMetrics)] <- singleCausalMgmSummaryMetrics
                    }
                }
            }
        }
    }
}



#### CausalCoxMGM precision recall curves (adjacency & orientation) across degree

p <- 110
N <- '00500'
alphas <- c(0.001, 0.005, 0.01, 0.05, 0.1)

causalmgmByDegreeMetrics <- data.frame()
causalmgmByDegreeSummaryMetrics <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (d in c(2,4,6)) {
            for (idx in 0:19) {
                ## ig.path <- readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))

                for (select in c('StEPS')) {
                    causal.graphs <- list()
                    for (alpha in alphas) {
                        causal.graphs[[which(alphas==alpha)]] <- tryCatch(expr={
                            loadGraph(
                                paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                                       '_', idx, '_', cr, '_', N, '.a',
                                       gsub('0[.]', '', as.character(alpha)),
                                       '.txt'))
                        }, error= function(cond) { NULL })
                    }
                    if (length(causal.graphs) > 0) {
                        singleCausalMgmMetrics <- do.call(rbind.data.frame, mclapply(causal.graphs, allMetricsByEdgeType, true=true, mc.cores=5))
                        singleCausalMgmMetrics$`Graph Type` <- gt
                        singleCausalMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                        singleCausalMgmMetrics$`Sample Size` <- as.numeric(N)
                        singleCausalMgmMetrics$Degree <- d
                        singleCausalMgmMetrics$`Number of Features` <- p
                        singleCausalMgmMetrics$`CoxMGM Select` <- select
                        singleCausalMgmMetrics$`Index` <- idx
                        causalmgmByDegreeMetrics[nrow(causalmgmByDegreeMetrics)+1:nrow(singleCausalMgmMetrics),colnames(singleCausalMgmMetrics)] <- singleCausalMgmMetrics
                        
                        singleCausalMgmSummaryMetrics <- data.frame()
                        for (edgetype in edgetypes[-8]) {
                            singleCausalMgmSummaryMetrics <- rbind(
                                singleCausalMgmSummaryMetrics,
                                data.frame(`Edge Type`=edgetype,
                                           skeletonAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'skeletonPrecision']),
                                           orientationAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'orientationRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'orientationPrecision']),
                                           causalAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'causalRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'causalPrecision']),
                                           check.names=F))
                        }
                        singleCausalMgmSummaryMetrics$`Graph Type` <- gt
                        singleCausalMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                        singleCausalMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                        singleCausalMgmSummaryMetrics$Degree <- d
                        singleCausalMgmSummaryMetrics$`Number of Features` <- p
                        singleCausalMgmSummaryMetrics$`CoxMGM Select` <- select
                        singleCausalMgmSummaryMetrics$`Index` <- idx
                        causalmgmByDegreeSummaryMetrics[nrow(causalmgmByDegreeSummaryMetrics)+1:7,colnames(singleCausalMgmSummaryMetrics)] <- singleCausalMgmSummaryMetrics
                    }
                }
            }
        }
    }
}


#### CausalCoxMGM precision recall curves (adjacency & orientation) across number of features

d <- 4
N <- '00500'
alphas <- c(0.001, 0.005, 0.01, 0.05, 0.1)

causalmgmByNumFeatsMetrics <- data.frame()
causalmgmByNumFeatsSummaryMetrics <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (p in numFeats) {
            for (idx in 0:19) {
                if (idx >= 10 && p==550) {
                    next
                }
                ## ig.path <- readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))

                for (select in c('StEPS')) {
                    causal.graphs <- list()
                    for (alpha in alphas) {
                        causal.graphs[[which(alphas==alpha)]] <- tryCatch(expr={
                            loadGraph(
                                paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                                       '_', idx, '_', cr, '_', N, '.a',
                                       gsub('0[.]', '', as.character(alpha)),
                                       '.txt'))
                        }, error= function(cond) { NULL })
                    }
                    if (length(causal.graphs) > 0) {
                        singleCausalMgmMetrics <- do.call(rbind.data.frame, mclapply(causal.graphs, allMetricsByEdgeType, true=true, mc.cores=5))
                        singleCausalMgmMetrics$`Graph Type` <- gt
                        singleCausalMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                        singleCausalMgmMetrics$`Sample Size` <- as.numeric(N)
                        singleCausalMgmMetrics$Degree <- d
                        singleCausalMgmMetrics$`Number of Features` <- p
                        singleCausalMgmMetrics$`CoxMGM Select` <- select
                        singleCausalMgmMetrics$`Index` <- idx
                        causalmgmByNumFeatsMetrics[nrow(causalmgmByNumFeatsMetrics)+1:nrow(singleCausalMgmMetrics),colnames(singleCausalMgmMetrics)] <- singleCausalMgmMetrics
                        
                        singleCausalMgmSummaryMetrics <- data.frame()
                        for (edgetype in edgetypes[-8]) {
                            singleCausalMgmSummaryMetrics <- rbind(
                                singleCausalMgmSummaryMetrics,
                                data.frame(`Edge Type`=edgetype,
                                           skeletonAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'skeletonPrecision']),
                                           orientationAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'orientationRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'orientationPrecision']),
                                           causalAUPR=calcAuPR(
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'causalRecall'],
                                               singleCausalMgmMetrics[singleCausalMgmMetrics$`Edge Type`==edgetype,'causalPrecision']),
                                           check.names=F))
                        }
                        singleCausalMgmSummaryMetrics$`Graph Type` <- gt
                        singleCausalMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                        singleCausalMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                        singleCausalMgmSummaryMetrics$Degree <- d
                        singleCausalMgmSummaryMetrics$`Number of Features` <- p
                        singleCausalMgmSummaryMetrics$`CoxMGM Select` <- select
                        singleCausalMgmSummaryMetrics$`Index` <- idx
                        causalmgmByNumFeatsSummaryMetrics[nrow(causalmgmByNumFeatsSummaryMetrics)+1:7,colnames(singleCausalMgmSummaryMetrics)] <- singleCausalMgmSummaryMetrics
                    }
                }
            }
        }
    }
}



ggplot(causalmgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `CoxMGM Select`) %>%
       summarize_at(c('skeletonRecall', 'skeletonPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=skeletonPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=skeletonPrecision_mean - 1.96 * skeletonPrecision_se,
           ymax=skeletonPrecision_mean + 1.96 * skeletonPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Adjacency Precision') +
    xlab('Adjacency Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave('causalcoxmgm_adj_pr_curves_bySampSize.png', width=6, height=3, dpi=400)


ggplot(causalmgmByNumFeatsMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD'), `Sample Size`==500) %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Number of Features`, `CoxMGM Select`) %>%
       summarize_at(c('skeletonRecall', 'skeletonPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=skeletonPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=skeletonPrecision_mean - 1.96 * skeletonPrecision_se,
           ymax=skeletonPrecision_mean + 1.96 * skeletonPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Adjacency Precision') +
    xlab('Adjacency Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Number of Features`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave('causalcoxmgm_adj_pr_curves_byNumFeats.png', width=6, height=3, dpi=400)

ggplot(causalmgmByDegreeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD'), `Sample Size`<=1000) %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Degree`, `CoxMGM Select`) %>%
       summarize_at(c('skeletonRecall', 'skeletonPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=skeletonPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=skeletonPrecision_mean - 1.96 * skeletonPrecision_se,
           ymax=skeletonPrecision_mean + 1.96 * skeletonPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Adjacency Precision') +
    xlab('Adjacency Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Degree`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave('causalcoxmgm_adj_pr_curves_byDegree.png', width=6, height=3, dpi=400)


ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `CoxMGM Select`, `Sample Size`) %>%
       summarize_at(c('skeletonAUPR', 'orientationAUPR'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=skeletonAUPR_mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=skeletonAUPR_mean - 1.96 * skeletonAUPR_se,
           ymax=skeletonAUPR_mean + 1.96 * skeletonAUPR_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


ggplot(causalmgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `CoxMGM Select`, `Sample Size`, Degree) %>%
       summarize_at(c('skeletonAUPR', 'orientationAUPR'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=skeletonAUPR_mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=skeletonAUPR_mean - 1.96 * skeletonAUPR_se,
           ymax=skeletonAUPR_mean + 1.96 * skeletonAUPR_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(causalmgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              !(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `CoxMGM Select`, `Sample Size`, `Number of Features`) %>%
       summarize_at(c('skeletonAUPR', 'orientationAUPR'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=skeletonAUPR_mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=skeletonAUPR_mean - 1.96 * skeletonAUPR_se,
           ymax=skeletonAUPR_mean + 1.96 * skeletonAUPR_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))




ggplot(causalmgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD'),
              `Sample Size` %in% c(100, 250, 500)) %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `CoxMGM Select`) %>%
       summarize_at(c('skeletonRecall', 'skeletonPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=skeletonPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=skeletonPrecision_mean - 1.96 * skeletonPrecision_se,
           ymax=skeletonPrecision_mean + 1.96 * skeletonPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`CoxMGM Select`)) +
    geom_pointrange(fatten=4) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Adjacency Precision') +
    xlab('Adjacency Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    theme_bw()

ggsave('causalcoxmgm_adj_pr_curves_reduced.png', width=9, height=10, dpi=400)


p1 <- ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(causalmgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(causalmgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.25,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_adjacency_aupr_curves_reduced.png', width=6, height=4, dpi=400)

p1 <- ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(causalmgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(causalmgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.25,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_adjacency_aupr_curves_all.png', width=6, height=3, dpi=400)


p1 <- ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(causalmgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.25,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(causalmgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('skeletonAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ## ymin=mean - 1.96 * se,
           ## ymax=mean + 1.96 * se,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Adjacency AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.25,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_adjacency_aupr_curves_surv.png', width=6, height=4, dpi=400)




p1 <- ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Orientation AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.225,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(causalmgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Causal AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.225,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(causalmgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Causal AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.225,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_causal_aupr_curves_reduced.png', width=6, height=4, dpi=400)


p1 <- ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Orientation AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.225,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(causalmgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Causal AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.225,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(causalmgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Causal AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.225,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_causal_aupr_curves_all.png', width=6, height=3, dpi=400)


p1 <- ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Orientation AUPRC') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.225,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(causalmgmByDegreeSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Causal AUPRC') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.225,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(causalmgmByNumFeatsSummaryMetrics %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       summarize_at('causalAUPR',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Causal AUPRC') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0.225,1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_causal_aupr_curves_surv.png', width=6, height=4, dpi=400)



normSHD <- function(df) {
    numC <- df$`Number of Features` * 5/11
    numD <- df$`Number of Features` * 5/11
    numS <- df$`Number of Features` * 1/11

    df[,'Normalized SHD'] <- NA

    df[df$`Edge Type`=='ALL','Normalized SHD'] <- 2 * df[df$`Edge Type`=='ALL','SHD'] / (df[df$`Edge Type`=='ALL','Number of Features'] * (df[df$`Edge Type`=='ALL','Number of Features']))

    df[df$`Edge Type`%in%c('SC','SD'),'Normalized SHD'] <- df[df$`Edge Type`%in%c('SC','SD'),'SHD'] / (numC[df$`Edge Type`%in%c('SC','SD')]*numS[df$`Edge Type`%in%c('SC','SD')])

    df[df$`Edge Type`=='SURV','Normalized SHD'] <- df[df$`Edge Type`=='SURV','SHD'] / ((numC[df$`Edge Type`=='SURV']+numD[df$`Edge Type`=='SURV'])*numS[df$`Edge Type`=='SURV'])

    df[df$`Edge Type`%in%c('CC','DD'),'Normalized SHD'] <- 2 * df[df$`Edge Type`%in%c('CC','DD'),'SHD'] / (numC[df$`Edge Type`%in%c('CC','DD')]*(numC[df$`Edge Type`%in%c('CC','DD')]-1))

    df[df$`Edge Type`=='CD','Normalized SHD'] <- df[df$`Edge Type`=='CD','SHD'] / (numC[df$`Edge Type`=='CD']*numD[df$`Edge Type`=='CD'])
    
    return(df)
}


normSHD <- function(df) {
    numC <- df$`Number of Features` * 5/11
    numD <- df$`Number of Features` * 5/11
    numS <- df$`Number of Features` * 1/11

    df[,'Normalized SHD'] <- NA

    df[df$`Edge Type`=='ALL','Normalized SHD'] <- 2 * df[df$`Edge Type`=='ALL','SHD'] / (df[df$`Edge Type`=='ALL','Number of Features'] * df[df$`Edge Type`=='ALL','Degree'])

    df[df$`Edge Type`%in%c('SC','SD'),'Normalized SHD'] <- 2 * df[df$`Edge Type`%in%c('SC','SD'),'SHD'] / (numS[df$`Edge Type`%in%c('SC','SD')] * df[df$`Edge Type`%in%c('SC','SD'),'Degree'])

    df[df$`Edge Type`=='SURV','Normalized SHD'] <- df[df$`Edge Type`=='SURV','SHD'] / (numS[df$`Edge Type`=='SURV'] * df[df$`Edge Type`=='SURV','Degree'])

    df[df$`Edge Type`%in%c('CC','DD','CD'),'Normalized SHD'] <- 2 * df[df$`Edge Type`%in%c('CC','DD', 'CD'),'SHD'] / (numC[df$`Edge Type`%in%c('CC','DD','CD')] * df[df$`Edge Type`%in%c('CC','DD','CD'),'Degree'])
    
    return(df)
}

p1 <- ggplot(normSHD(causalmgmBySampSizeMetrics) %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('Normalized SHD') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0,2.1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(normSHD(causalmgmByDegreeMetrics) %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('SHD') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0,2.1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(normSHD(causalmgmByNumFeatsMetrics) %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('SHD') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0,2.1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_normshd_curves_reduced.png', width=6, height=4, dpi=600)


p1 <- ggplot(normSHD(causalmgmBySampSizeMetrics) %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('Normalized SHD') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0,2.1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(normSHD(causalmgmByDegreeMetrics) %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('SHD') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0,2.1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(normSHD(causalmgmByNumFeatsMetrics) %>%
       filter(`Edge Type` %in% c('ALL')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('SHD') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0,2.1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_normshd_curves_all.png', width=6, height=3, dpi=600)



p1 <- ggplot(normSHD(causalmgmBySampSizeMetrics) %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ## xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           ## xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('Normalized SHD') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0,2.1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(normSHD(causalmgmByDegreeMetrics) %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(`Edge Type`!='NoSURV') %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, Degree) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('SHD') +
    xlab('Degree') +
    facet_grid(rows=vars(`Edge Type`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0,2.1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

p2 <- ggplot(normSHD(causalmgmByNumFeatsMetrics) %>%
       filter(`Edge Type` %in% c('SC', 'SD')) %>%
       ## filter(!(`Number of Features`==550 & `Graph Type`=='SF')) %>%
       group_by(`Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `Number of Features`) %>%
       ## mutate(`Normalized SHD`=SHD / (Degree * `Number of Features` / 2)) %>%
       summarize_at('Normalized SHD',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    ## geom_errorbar() +
    geom_line() +
    ylab('SHD') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Edge Type`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    ylim(c(0,2.1)) +
    theme_bw(base_size=7) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('causalcoxmgm_normshd_curves_surv.png', width=6, height=4, dpi=600)



ggplot(causalmgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
              ## `Sample Size` %in% c(100, 250, 500, 1000)) %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `CoxMGM Select`) %>%
       summarize_at(c('orientationRecall', 'orientationPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=orientationRecall_mean, y=orientationPrecision_mean,
           xmin=orientationRecall_mean - 1.96 * orientationRecall_se,
           xmax=orientationRecall_mean + 1.96 * orientationRecall_se,
           ymin=orientationPrecision_mean - 1.96 * orientationPrecision_se,
           ymax=orientationPrecision_mean + 1.96 * orientationPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Orientation Precision') +
    xlab('Orientation Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('causalcoxmgm_orient_pr_curves_reduced.png', width=7, height=6, dpi=400)

ggplot(causalmgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')) %>%
              ## `Sample Size` %in% c(100, 250, 500, 1000)) %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `CoxMGM Select`) %>%
       summarize_at(c('causalRecall', 'causalPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=causalRecall_mean, y=causalPrecision_mean,
           xmin=causalRecall_mean - 1.96 * causalRecall_se,
           xmax=causalRecall_mean + 1.96 * causalRecall_se,
           ymin=causalPrecision_mean - 1.96 * causalPrecision_se,
           ymax=causalPrecision_mean + 1.96 * causalPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Causal Precision') +
    xlab('Causal Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('causalcoxmgm_causal_pr_curves_reduced.png', width=7, height=6, dpi=400)


ggplot(causalmgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500),
              `CoxMGM Select`=='StEPS') %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `CoxMGM Select`) %>%
       summarize_at(c('skeletonRecall', 'skeletonPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=skeletonPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=skeletonPrecision_mean - 1.96 * skeletonPrecision_se,
           ymax=skeletonPrecision_mean + 1.96 * skeletonPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`CoxMGM Select`)) +
    geom_pointrange(fatten=4) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Adjacency Precision') +
    xlab('Adjacency Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    theme_bw()

ggsave('causalcoxmgm_adj_pr_curves_reduced_StEPS.png', width=9, height=10, dpi=400)

ggplot(causalmgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500),
              `CoxMGM Select`=='StEPS') %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `CoxMGM Select`) %>%
       summarize_at(c('orientationRecall', 'orientationPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=orientationRecall_mean, y=orientationPrecision_mean,
           xmin=orientationRecall_mean - 1.96 * orientationRecall_se,
           xmax=orientationRecall_mean + 1.96 * orientationRecall_se,
           ymin=orientationPrecision_mean - 1.96 * orientationPrecision_se,
           ymax=orientationPrecision_mean + 1.96 * orientationPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`CoxMGM Select`)) +
    geom_pointrange(fatten=4) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Orientation Precision') +
    xlab('Orientation Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    theme_bw()

ggsave('causalcoxmgm_orient_pr_curves_reduced_StEPS.png', width=9, height=10, dpi=400)

ggplot(causalmgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500),
              `CoxMGM Select`=='StEPS') %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `CoxMGM Select`) %>%
       summarize_at(c('directionalRecall', 'directionalPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=directionalRecall_mean, y=directionalPrecision_mean,
           xmin=directionalRecall_mean - 1.96 * directionalRecall_se,
           xmax=directionalRecall_mean + 1.96 * directionalRecall_se,
           ymin=directionalPrecision_mean - 1.96 * directionalPrecision_se,
           ymax=directionalPrecision_mean + 1.96 * directionalPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`CoxMGM Select`)) +
    geom_pointrange(fatten=4) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Directional Precision') +
    xlab('Directional Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    theme_bw()

ggsave('causalcoxmgm_directional_pr_curves_reduced_StEPS.png', width=9, height=10, dpi=400)



ggplot(mgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV')),
       aes(`Censoring Rate`, AUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position='dodge') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('AUPRC') +
    theme_bw()

ggsave('coxmgm_AUPRC_by_samplesize.png', width=15, height=10, dpi=400)

ggplot(mgmBySampSizeMetrics %>%
       group_by(Lambda, `Graph Type`, `Edge Type`, `Censoring Rate`, `Sample Size`) %>%
       summarize_at(c('skeletonRecall', 'moralPrecision'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=skeletonRecall_mean, y=moralPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_se,
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_se,
           ymin=moralPrecision_mean - 1.96 * moralPrecision_se,
           ymax=moralPrecision_mean + 1.96 * moralPrecision_se,
           color=`Graph Type`, lty=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_errorbarh(height=0) +
    geom_line() +
    ylab('Precision') +
    xlab('Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    theme_bw()

ggsave('coxmgm_pr_curves_by_samplesize_all_edgetypes.png', width=15, height=10, dpi=400)

library(ggpattern)

ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500, 1000)),
       aes(`Censoring Rate`, skeletonAUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge',
                     pattern_fill='white', pattern_key_scale_factor = 0.5,
                     pattern_density=0.3) +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=1), width=0.5) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('Adjacency AUPRC') +
    xlab("% Censored") + 
    theme_bw()

ggsave('causalcoxmgm_adj_AUPRC_reduced.png', width=5, height=6, dpi=400)


ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500, 1000)),
       aes(`Censoring Rate`, orientationAUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=1), width=0.5) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('Orientation AUPRC') +
    xlab("% Censored") + 
    theme_bw()

ggsave('causalcoxmgm_orient_AUPRC_reduced.png', width=5, height=6, dpi=400)

ggplot(causalmgmBySampSizeSummaryMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500, 1000)),
       aes(`Censoring Rate`, causalAUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=1), width=0.5) +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`)) +
    ylab('Causal AUPRC') +
    xlab("% Censored") + 
    theme_bw()

ggsave('causalcoxmgm_causal_AUPRC_reduced.png', width=5, height=6, dpi=400)


ggplot(causalmgmBySampSizeMetrics %>%
       filter(`Edge Type` %in% c('ALL', 'SC', 'SD', 'SURV'),
              `Sample Size` %in% c(100, 250, 500, 1000)) %>%
       group_by(Alpha, `Graph Type`, `Edge Type`, `Censoring Rate`,
                `Sample Size`, `CoxMGM Select`) %>%
       summarize_at(c('SHD'),
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=log10(Alpha), y=mean,
           ymin=mean - 1.96 * se,
           ymax=mean + 1.96 * se,
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('SHD') +
    ## xlab('Orientation Recall') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Sample Size`), scales='free_y') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('causalcoxmgm_SHD_reduced.png', width=7, height=6, dpi=400)




#### CausalCoxMGM vs LASSO Feature Selection across sample size

p <- 110
d <- 4
alphas <- c(0.001, 0.005, 0.01, 0.05, 0.1)

featSelectBySampSizeMetrics <- data.frame()

for (gt in rev(graphTypes)) {
    for (N in sampSizes) {
        for (cr in censorRates) {
            for (idx in 0:19) {
                ## ig.path <- readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))

                ## lassoRes <- readRDS(paste0('out/coxLASSO_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))

                ## singleFeatSelectMetrics <- lassoFeatSelectEval(true, lassoRes)
                ## singleFeatSelectMetrics$`Graph Type` <- gt
                ## singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                ## singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                ## singleFeatSelectMetrics$Degree <- d
                ## singleFeatSelectMetrics$`Number of Features` <- p
                ## singleFeatSelectMetrics$`CoxMGM Select` <- NA
                ## singleFeatSelectMetrics$`Index` <- idx
                ## featSelectBySampSizeMetrics <- rbind(featSelectBySampSizeMetrics, singleFeatSelectMetrics)
                

                ## for (select in c('StEPS')) {
                ## causal.graphs <- list()
                ## for (alpha in alphas) {
                ##     causal.graphs[[which(alphas==alpha)]] <- loadGraph(
                ##         paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                ##                '_', idx, '_', cr, '_', N, '.a',
                ##                gsub('0[.]', '', as.character(alpha)),
                ##                '.txt'))
                ## }
                select <- "StEPS"
                causal.graphs <- list()
                for (alpha in alphas) {
                    causal.graphs[[which(alphas==alpha)]] <- tryCatch(expr={
                        loadGraph(
                            paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                                   '_', idx, '_', cr, '_', N, '.a',
                                   gsub('0[.]', '', as.character(alpha)),
                                   '.txt'))
                    }, error= function(cond) { NULL })
                }

                if (length(causal.graphs)>0) {
                    
                    singleFeatSelectMetrics <- do.call(rbind.data.frame, mclapply(causal.graphs, causalFeatSelectEval, true=true, mc.cores=5))
                    singleFeatSelectMetrics$`Graph Type` <- gt
                    singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                    singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                    singleFeatSelectMetrics$Degree <- d
                    singleFeatSelectMetrics$`Number of Features` <- p
                    singleFeatSelectMetrics$`CoxMGM Select` <- select
                    singleFeatSelectMetrics$`Index` <- idx
                    featSelectBySampSizeMetrics <- rbind(featSelectBySampSizeMetrics, singleFeatSelectMetrics)

                    lassoRes <- readRDS(paste0('out/coxLASSO_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))

                    singleFeatSelectMetrics <- lassoFeatSelectEval(true, lassoRes)
                    singleFeatSelectMetrics$`Graph Type` <- gt
                    singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                    singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                    singleFeatSelectMetrics$Degree <- d
                    singleFeatSelectMetrics$`Number of Features` <- p
                    singleFeatSelectMetrics$`CoxMGM Select` <- NA
                    singleFeatSelectMetrics$`Index` <- idx
                    featSelectBySampSizeMetrics <- rbind(featSelectBySampSizeMetrics, singleFeatSelectMetrics)
                }
            }
        }
    }
}


featSelectBySampSizeMetrics$`Feature Select` <- factor(featSelectBySampSizeMetrics$`Feature Select`, levels=rev(c('Min', '1SE', 'MB', 'NB')), labels=rev(c('LASSO Min', 'LASSO 1SE', 'Causal MB', 'Causal NB')))


#### CausalCoxMGM vs LASSO Feature Selection across number of features

N <- '00500'
d <- 4
alphas <- c(0.001, 0.005, 0.01, 0.05, 0.1)

featSelectByNumFeatsMetrics <- data.frame()

for (gt in rev(graphTypes)) {
    for (p in numFeats) {
        for (cr in censorRates) {
            for (idx in 0:19) {
                if (p==550 && idx >= 10) next
                ## ig.path <- readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))

                ## lassoRes <- readRDS(paste0('out/coxLASSO_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))

                ## singleFeatSelectMetrics <- lassoFeatSelectEval(true, lassoRes)
                ## singleFeatSelectMetrics$`Graph Type` <- gt
                ## singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                ## singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                ## singleFeatSelectMetrics$Degree <- d
                ## singleFeatSelectMetrics$`Number of Features` <- p
                ## singleFeatSelectMetrics$`CoxMGM Select` <- NA
                ## singleFeatSelectMetrics$`Index` <- idx
                ## featSelectByNumFeatsMetrics <- rbind(featSelectByNumFeatsMetrics, singleFeatSelectMetrics)
                

                ## for (select in c('StEPS')) {
                ## causal.graphs <- list()
                ## for (alpha in alphas) {
                ##     causal.graphs[[which(alphas==alpha)]] <- loadGraph(
                ##         paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                ##                '_', idx, '_', cr, '_', N, '.a',
                ##                gsub('0[.]', '', as.character(alpha)),
                ##                '.txt'))
                ## }
                select <- "StEPS"
                causal.graphs <- list()
                for (alpha in alphas) {
                    causal.graphs[[which(alphas==alpha)]] <- tryCatch(expr={
                        loadGraph(
                            paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                                   '_', idx, '_', cr, '_', N, '.a',
                                   gsub('0[.]', '', as.character(alpha)),
                                   '.txt'))
                    }, error= function(cond) { NULL })
                }

                if (length(causal.graphs)>0) {
                    
                    singleFeatSelectMetrics <- do.call(rbind.data.frame, mclapply(causal.graphs, causalFeatSelectEval, true=true, mc.cores=5))
                    singleFeatSelectMetrics$`Graph Type` <- gt
                    singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                    singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                    singleFeatSelectMetrics$Degree <- d
                    singleFeatSelectMetrics$`Number of Features` <- p
                    singleFeatSelectMetrics$`CoxMGM Select` <- select
                    singleFeatSelectMetrics$`Index` <- idx
                    featSelectByNumFeatsMetrics <- rbind(featSelectByNumFeatsMetrics, singleFeatSelectMetrics)

                    lassoRes <- readRDS(paste0('out/coxLASSO_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))

                    singleFeatSelectMetrics <- lassoFeatSelectEval(true, lassoRes)
                    singleFeatSelectMetrics$`Graph Type` <- gt
                    singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                    singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                    singleFeatSelectMetrics$Degree <- d
                    singleFeatSelectMetrics$`Number of Features` <- p
                    singleFeatSelectMetrics$`CoxMGM Select` <- NA
                    singleFeatSelectMetrics$`Index` <- idx
                    featSelectByNumFeatsMetrics <- rbind(featSelectByNumFeatsMetrics, singleFeatSelectMetrics)
                }
            }
        }
    }
}


featSelectByNumFeatsMetrics$`Feature Select` <- factor(featSelectByNumFeatsMetrics$`Feature Select`, levels=rev(c('Min', '1SE', 'MB', 'NB')), labels=rev(c('LASSO Min', 'LASSO 1SE', 'Causal MB', 'Causal NB')))

#### CausalCoxMGM vs LASSO Feature Selection across number of features

N <- '00500'
p <- 110
alphas <- c(0.001, 0.005, 0.01, 0.05, 0.1)

featSelectByDegreeMetrics <- data.frame()

for (gt in rev(graphTypes)) {
    for (d in degrees) {
        for (cr in censorRates) {
            for (idx in 0:19) {
                if (p==550 && idx >= 10) next
                ## ig.path <- readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))

                ## lassoRes <- readRDS(paste0('out/coxLASSO_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))

                ## singleFeatSelectMetrics <- lassoFeatSelectEval(true, lassoRes)
                ## singleFeatSelectMetrics$`Graph Type` <- gt
                ## singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                ## singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                ## singleFeatSelectMetrics$Degree <- d
                ## singleFeatSelectMetrics$`Number of Features` <- p
                ## singleFeatSelectMetrics$`CoxMGM Select` <- NA
                ## singleFeatSelectMetrics$`Index` <- idx
                ## featSelectByDegreeMetrics <- rbind(featSelectByDegreeMetrics, singleFeatSelectMetrics)
                

                ## for (select in c('StEPS')) {
                ## causal.graphs <- list()
                ## for (alpha in alphas) {
                ##     causal.graphs[[which(alphas==alpha)]] <- loadGraph(
                ##         paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                ##                '_', idx, '_', cr, '_', N, '.a',
                ##                gsub('0[.]', '', as.character(alpha)),
                ##                '.txt'))
                ## }
                select <- "StEPS"
                causal.graphs <- list()
                for (alpha in alphas) {
                    causal.graphs[[which(alphas==alpha)]] <- tryCatch(expr={
                        loadGraph(
                            paste0('out/mpccoxmgm', select, '_', gt, '_deg', d, '_p', p,
                                   '_', idx, '_', cr, '_', N, '.a',
                                   gsub('0[.]', '', as.character(alpha)),
                                   '.txt'))
                    }, error= function(cond) { NULL })
                }

                if (length(causal.graphs)>0) {
                    
                    singleFeatSelectMetrics <- do.call(rbind.data.frame, mclapply(causal.graphs, causalFeatSelectEval, true=true, mc.cores=5))
                    singleFeatSelectMetrics$`Graph Type` <- gt
                    singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                    singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                    singleFeatSelectMetrics$Degree <- d
                    singleFeatSelectMetrics$`Number of Features` <- p
                    singleFeatSelectMetrics$`CoxMGM Select` <- select
                    singleFeatSelectMetrics$`Index` <- idx
                    featSelectByDegreeMetrics <- rbind(featSelectByDegreeMetrics, singleFeatSelectMetrics)

                    lassoRes <- readRDS(paste0('out/coxLASSO_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))

                    singleFeatSelectMetrics <- lassoFeatSelectEval(true, lassoRes)
                    singleFeatSelectMetrics$`Graph Type` <- gt
                    singleFeatSelectMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                    singleFeatSelectMetrics$`Sample Size` <- as.numeric(N)
                    singleFeatSelectMetrics$Degree <- d
                    singleFeatSelectMetrics$`Number of Features` <- p
                    singleFeatSelectMetrics$`CoxMGM Select` <- NA
                    singleFeatSelectMetrics$`Index` <- idx
                    featSelectByDegreeMetrics <- rbind(featSelectByDegreeMetrics, singleFeatSelectMetrics)
                }
            }
        }
    }
}


featSelectByDegreeMetrics$`Feature Select` <- factor(featSelectByDegreeMetrics$`Feature Select`, levels=rev(c('Min', '1SE', 'MB', 'NB')), labels=rev(c('LASSO Min', 'LASSO 1SE', 'Causal MB', 'Causal NB')))


ggplot(featSelectBySampSizeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)),
       aes(`Censoring Rate`, F1, fill=`Feature Select`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=0.9), width=0.5) +
    facet_grid(rows=vars(`Graph Type`), cols=vars(`Sample Size`)) +
    scale_fill_brewer(palette='Paired') +
    labs(fill="Method") +
    theme_bw(base_size=7)

ggsave('featSelect_bySampSize.png', width=6, height=4, dpi=400)


ggplot(featSelectByNumFeatsMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)),
       aes(`Censoring Rate`, F1, fill=`Feature Select`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=0.9), width=0.5) +
    facet_grid(rows=vars(`Graph Type`), cols=vars(`Number of Features`)) +
    scale_fill_brewer(palette='Paired') +
    labs(fill="Method") +
    theme_bw(base_size=7) +
    theme(legend.position="none")

ggsave('featSelect_byNumFeats.png', width=3, height=4, dpi=400)


ggplot(featSelectByDegreeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)),
       aes(`Censoring Rate`, F1, fill=`Feature Select`)) +
    geom_bar(stat='summary', fun.data='mean_se', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96),
                  position=position_dodge(width=0.9), width=0.5) +
    facet_grid(rows=vars(`Graph Type`), cols=vars(Degree)) +
    scale_fill_brewer(palette='Paired') +
    labs(fill="Method") +
    theme_bw(base_size=7) +
    theme(legend.position="none")

ggsave('featSelect_byDegree.png', width=3, height=4, dpi=400)


ggplot(featSelectBySampSizeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)),
       aes(x=`Sample Size`, y=F1, shape=`Censoring Rate`, lty=`Censoring Rate`, color=`Graph Type`)) +
    geom_line(stat='summary', fun.data='mean_se') +
    geom_pointrange(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96)) +
    facet_grid(rows=vars(`Feature Select`)) +
    scale_x_continuous(trans='log10') +
    labs(fill="Method") +
    theme_bw(base_size=7)



p1 <- ggplot(featSelectBySampSizeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, `Sample Size`) %>%
       summarize_at('F1',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('F1') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Feature Select`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.15, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(featSelectByDegreeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, Degree) %>%
       summarize_at('F1',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('F1') +
    xlab('Degree') +
    facet_grid(rows=vars(`Feature Select`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.15, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


p2 <- ggplot(featSelectByNumFeatsMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, `Number of Features`) %>%
       summarize_at('F1',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('F1') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Feature Select`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.15, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())


plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('featSelect_F1_curves.png', width=6, height=4, dpi=400)



p1 <- ggplot(featSelectBySampSizeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, `Sample Size`) %>%
       summarize_at('Precision',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Precision') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Feature Select`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.1, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(featSelectByDegreeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, Degree) %>%
       summarize_at('Precision',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Precision') +
    xlab('Degree') +
    facet_grid(rows=vars(`Feature Select`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.1, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


p2 <- ggplot(featSelectByNumFeatsMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, `Number of Features`) %>%
       summarize_at('Precision',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Precision') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Feature Select`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.1, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())


plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('featSelect_Precision_curves.png', width=6, height=4, dpi=400)



p1 <- ggplot(featSelectBySampSizeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, `Sample Size`) %>%
       summarize_at('Recall',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Sample Size`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Recall') +
    xlab('Sample Size') +
    facet_grid(rows=vars(`Feature Select`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.1, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank())


p3 <- ggplot(featSelectByDegreeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, Degree) %>%
       summarize_at('Recall',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=Degree, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Recall') +
    xlab('Degree') +
    facet_grid(rows=vars(`Feature Select`)) +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.1, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


p2 <- ggplot(featSelectByNumFeatsMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`)) %>%
       group_by(Alpha, `CoxMGM Select`, `Feature Select`,
                `Censoring Rate`, `Graph Type`, `Number of Features`) %>%
       summarize_at('Recall',
                    c(mean=mean, se=function(x) sd(x)/sqrt(length(x)))),
       aes(x=`Number of Features`, y=mean,
           ymin=pmax(mean - 1.96 * se, 0),
           ymax=pmin(mean + 1.96 * se, 1),
           color=`Graph Type`, lty=`Censoring Rate`, shape=`Censoring Rate`)) +
    geom_pointrange(fatten=2) +
    geom_line() +
    ylab('Recall') +
    xlab('Number of Features') +
    facet_grid(rows=vars(`Feature Select`)) +
    scale_x_continuous(trans='log10') +
    labs(lty="% Censored", shape="% Censored") +
    theme_bw(base_size=7) +
    ylim(c(0.1, 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          axis.title.y=element_blank(), axis.text.y = element_blank(),
          legend.position='none', strip.background = element_blank(),
          strip.text = element_blank(), axis.ticks.y = element_blank())


plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.1,1,1.4))
ggsave('featSelect_Recall_curves.png', width=6, height=4, dpi=400)



ggplot(featSelectBySampSizeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`),
              `Feature Select` != 'MB'),
       aes(`Censoring Rate`, Precision, fill=`Feature Select`)) +
    geom_bar(stat='summary', fun.data='mean_sdl', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96), position='dodge') +
    facet_grid(rows=vars(`Graph Type`), cols=vars(`Sample Size`)) +
    scale_fill_brewer(palette='Paired') +
    theme_bw()

ggplot(featSelectBySampSizeMetrics %>%
       filter(Alpha==0.05 | is.na(Alpha),
              `CoxMGM Select`=='StEPS' | is.na(`CoxMGM Select`),
              `Feature Select` != 'MB'),
       aes(`Censoring Rate`, Recall, fill=`Feature Select`)) +
    geom_bar(stat='summary', fun.data='mean_sdl', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_se', fun.args=list(mult=1.96), position='dodge') +
    facet_grid(rows=vars(`Graph Type`), cols=vars(`Sample Size`)) +
    scale_fill_brewer(palette='Paired') +
    theme_classic()


#### CoxMGM precision recall curves across graph degree

p <- 110
## d <- 4
N <- '00500'

mgmByDegreeMetrics <- data.frame()
mgmByDegreeSummaryMetrics <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (d in degrees) {
            for (idx in 0:9) {
                ig.path <- readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))
                singleMgmMetrics <- do.call(rbind.data.frame, mclapply(ig.path$graphs, prMetricsByEdgeType, true=true, mc.cores=16))
                singleMgmMetrics$`Graph Type` <- gt
                singleMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmMetrics$`Sample Size` <- as.numeric(N)
                singleMgmMetrics$Degree <- d
                singleMgmMetrics$`Number of Features` <- p
                singleMgmMetrics$`Index` <- idx
                mgmByDegreeMetrics <- rbind(mgmByDegreeMetrics, singleMgmMetrics)
                
                singleMgmSummaryMetrics <- data.frame()
                for (edgetype in edgetypes) {
                    singleMgmSummaryMetrics <- rbind(
                        singleMgmSummaryMetrics,
                        data.frame(`Edge Type`=edgetype,
                                   AUPR=calcAuPR(
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'moralPrecision']),
                                   check.names=F))
                }
                singleMgmSummaryMetrics$`Graph Type` <- gt
                singleMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSummaryMetrics$Degree <- d
                singleMgmSummaryMetrics$`Number of Features` <- p
                singleMgmSummaryMetrics$`Index` <- idx
                mgmByDegreeSummaryMetrics <- rbind(mgmByDegreeSummaryMetrics,
                                                     singleMgmSummaryMetrics)
            }
        }
    }
}


ggplot(mgmByDegreeMetrics %>%
       group_by(Lambda, `Graph Type`, `Edge Type`, `Censoring Rate`, Degree) %>%
       summarize_at(c('skeletonRecall', 'moralPrecision'), c(mean=mean, sd=sd)),
       aes(x=skeletonRecall_mean, y=moralPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_sd / sqrt(10),
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_sd / sqrt(10),
           ymin=moralPrecision_mean - 1.96 * moralPrecision_sd / sqrt(10),
           ymax=moralPrecision_mean + 1.96 * moralPrecision_sd / sqrt(10),
           color=`Graph Type`, lty=`Censoring Rate`)) +
    geom_errorbar() +
    geom_errorbarh() +
    geom_line() +
    facet_grid(rows=vars(`Edge Type`), cols=vars(Degree)) +
    theme_classic()

ggplot(mgmByDegreeSummaryMetrics,
       aes(`Censoring Rate`, AUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_sdl', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_sdl', fun.args=list(mult=1.96), position='dodge') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(Degree)) +
    theme_classic()


#### CoxMGM precision recall curves across number of features

## p <- 110
d <- 4
N <- '00500'

mgmByNumFeatsMetrics <- data.frame()
mgmByNumFeatsSummaryMetrics <- data.frame()

for (gt in graphTypes) {
    for (cr in censorRates) {
        for (p in numFeats) {
            for (idx in 0:9) {
                ig.path <- readRDS(paste0('out/coxmgmPath_', gt, '_deg', d, '_p', p, '_', idx, '_', cr, '_', N, '.rds'))
                true <- loadGraph(paste0('true_out/sout_', gt, '_deg', d, '_p', p, '_', idx, '.txt'))
                singleMgmMetrics <- do.call(rbind.data.frame, mclapply(ig.path$graphs, prMetricsByEdgeType, true=true, mc.cores=16))
                singleMgmMetrics$`Graph Type` <- gt
                singleMgmMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmMetrics$`Sample Size` <- as.numeric(N)
                singleMgmMetrics$Degree <- d
                singleMgmMetrics$`Number of Features` <- p
                singleMgmMetrics$`Index` <- idx
                mgmByNumFeatsMetrics <- rbind(mgmByNumFeatsMetrics, singleMgmMetrics)
                
                singleMgmSummaryMetrics <- data.frame()
                for (edgetype in edgetypes) {
                    singleMgmSummaryMetrics <- rbind(
                        singleMgmSummaryMetrics,
                        data.frame(`Edge Type`=edgetype,
                                   AUPR=calcAuPR(
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'skeletonRecall'],
                                       singleMgmMetrics[singleMgmMetrics$`Edge Type`==edgetype,'moralPrecision']),
                                   check.names=F))
                }
                singleMgmSummaryMetrics$`Graph Type` <- gt
                singleMgmSummaryMetrics$`Censoring Rate` <- paste0(100 * as.numeric(cr), '%')
                singleMgmSummaryMetrics$`Sample Size` <- as.numeric(N)
                singleMgmSummaryMetrics$Degree <- d
                singleMgmSummaryMetrics$`Number of Features` <- p
                singleMgmSummaryMetrics$`Index` <- idx
                mgmByNumFeatsSummaryMetrics <- rbind(mgmByNumFeatsSummaryMetrics,
                                                     singleMgmSummaryMetrics)
            }
        }
    }
}


ggplot(mgmByNumFeatsMetrics %>%
       group_by(Lambda, `Graph Type`, `Edge Type`, `Censoring Rate`, `Number of Features`) %>%
       summarize_at(c('skeletonRecall', 'moralPrecision'), c(mean=mean, sd=sd)),
       aes(x=skeletonRecall_mean, y=moralPrecision_mean,
           xmin=skeletonRecall_mean - 1.96 * skeletonRecall_sd / sqrt(10),
           xmax=skeletonRecall_mean + 1.96 * skeletonRecall_sd / sqrt(10),
           ymin=moralPrecision_mean - 1.96 * moralPrecision_sd / sqrt(10),
           ymax=moralPrecision_mean + 1.96 * moralPrecision_sd / sqrt(10),
           color=`Graph Type`, lty=`Censoring Rate`)) +
    geom_errorbar() +
    geom_errorbarh() +
    geom_line() +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Number of Features`)) +
    theme_classic()

ggplot(mgmByNumFeatsSummaryMetrics,
       aes(`Censoring Rate`, AUPR, fill=`Graph Type`)) +
    geom_bar(stat='summary', fun.data='mean_sdl', position='dodge') +
    geom_errorbar(stat='summary', fun.data='mean_sdl', fun.args=list(mult=1.96), position='dodge') +
    facet_grid(rows=vars(`Edge Type`), cols=vars(`Number of Features`)) +
    theme_classic()

