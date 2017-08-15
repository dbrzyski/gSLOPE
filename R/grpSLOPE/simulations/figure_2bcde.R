#---------------------------------------------------------------------------
# This code reproduces the simulation studies of Figures 2b-e in
# D. Brzyski, A.Gossmann, W. Su, M. Bogdan (2016), "Group SLOPE -
# - adaptive selection of groups of predictors"
# Preprint available at https://arxiv.org/abs/1610.04960
#---------------------------------------------------------------------------

library(grpSLOPE)
library(dplyr)

# Adjust the number of cores to your particular system
library(doParallel)
registerDoParallel(cores = as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")))


#--- Set up global parameters for the simulation

# auxilliary function to get (group-wise) FDP and power from one grpSLOPE solution object
get_FDP_and_power <- function(result, true.relevant){
  truepos <- length(intersect(result$selected, true.relevant))
  falsepos <- length(result$selected) - truepos
  FDP <- falsepos / max(1, length(result$selected))
  pow <- truepos / length(true.relevant)
  return(c("FDP" = FDP, "pow" = pow))
}

# set X to identity matrix, b/c gSLOPE with orthogonal design is equivalent to a problem with identity mat.
p <- 5000
X <- diag(rep(1,p))

# set the grouping structure
group <- c(rep(1:200, each=3),
           rep(201:400, each=4),
           rep(401:600, each=5),
           rep(601:800, each=6),
           rep(801:1000, each=7))
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)
n.group <- length(group.id)

# determine signal strength, such as used in Figure 2 in Brzyski et. al. (2016)
Bfun <- function(l) {
  sqrt(4*log(n.group) / (1 - n.group^(-2/l)) - l)
}
signal.strength <- sum(Bfun(group.length)) / sum(sqrt(group.length))

# considered numbers of truly relevant groups
n.relevant <- floor(seq(1, 250, length=11))

# how many times the simulation is repeated
n.replications <- 300


#--- Simulation main loop

set.seed(20170812)

# run the simulations n.replications times at each level of n.relevant,
# and with each combination of lambda="max", lambda="mean", fdr=0.1, fdr=0.05

parallel.results <- vector(mode="list")
for (k in 1:length(n.relevant)) {
  cat(paste("sparsity level", k, "started"))

  parallel.results[[k]] <- foreach(i=1:n.replications, .combine = rbind) %dopar% {
    cat(".")

    # generate coeffient vector, pick relevant groups at random
    b <- rep(0, p)
    n.signif <- n.relevant[k]
    ind.relevant <- sample(1:n.group, n.signif)
    for (j in ind.relevant) { b[group.id[[j]]] <- signal.strength }

    # generate the response vector
    y <- X %*% b + rnorm(p, sd=1)

    # get Group SLOPE solutions with different lambda and fdr values
    lambda.max.1 <- grpSLOPE(X=X, y=y, group=group, fdr=0.1,
                             lambda="max", sigma=1, verbose=FALSE,
                             orthogonalize=FALSE, normalize=FALSE)
    lambda.max.05 <- grpSLOPE(X=X, y=y, group=group, fdr=0.05,
                              lambda="max", sigma=1, verbose=FALSE,
                              orthogonalize=FALSE, normalize=FALSE)
    lambda.mean.1 <- grpSLOPE(X=X, y=y, group=group, fdr=0.1,
                              lambda="mean", sigma=1, verbose=FALSE,
                              orthogonalize=FALSE, normalize=FALSE)
    lambda.mean.05 <- grpSLOPE(X=X, y=y, group=group, fdr=0.05,
                               lambda="mean", sigma=1, verbose=FALSE,
                               orthogonalize=FALSE, normalize=FALSE)

    # get the FDPs and powers of the grpSLOPE solutions
    true.relevant <- names(group.id)[ind.relevant]
    FDPs.and.powers <- rep(NA, 10)
    names(FDPs.and.powers) <- c("n.relevant", "replication",
                                "lambda.max.1.FDP", "lambda.max.1.power",
                                "lambda.max.05.FDP", "lambda.max.05.power",
                                "lambda.mean.1.FDP", "lambda.mean.1.power",
                                "lambda.mean.05.FDP", "lambda.mean.05.power")
    FDPs.and.powers["n.relevant"] <- length(true.relevant)
    FDPs.and.powers["replication"] <- i
    FDPs.and.powers[c("lambda.max.1.FDP", "lambda.max.1.power")] <- get_FDP_and_power(lambda.max.1, true.relevant)
    FDPs.and.powers[c("lambda.max.05.FDP", "lambda.max.05.power")] <- get_FDP_and_power(lambda.max.05, true.relevant)
    FDPs.and.powers[c("lambda.mean.1.FDP", "lambda.mean.1.power")] <- get_FDP_and_power(lambda.mean.1, true.relevant)
    FDPs.and.powers[c("lambda.mean.05.FDP", "lambda.mean.05.power")] <- get_FDP_and_power(lambda.mean.05, true.relevant)

    # return the results
    FDPs.and.powers
  }

  cat("done\n")
}

# collect results in a data frame
results <- data.frame() %>% tbl_df
for(k in 1:length(n.relevant)) {
  results <- bind_rows(results, tbl_df(parallel.results[[k]]))
}

# save the results
save(results, file = "../RData/figure_2bcde_simulation_results.RData")
