#---------------------------------------------------------------------------
# This code reproduces the simulation studies of Figures 3a-c in
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

# set X to identity matrix, b/c gSLOPE with orthogonal design is equivalent to a problem with identity mat.
p <- 5000
X <- diag(rep(1,p))
fdr <- 0.05

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

    # NOTE: the convenience function grpSLOPE does not allow to change the weights
    # manually. So we will generate the lambda sequence, the weight vectors, and
    # apply the proximal gradient solver for Group SLOPE by hand as three separate steps.

    # generate lambda
    lambda.mean <- lambdaGroupSLOPE(fdr=fdr, group=group,
                                    wt=sqrt(group.length),
                                    method="mean")

    # considered vectors of weight per coefficient
    sqrt_wt <- rep(NA, p)
    for (j in 1:n.group) {
      sqrt_wt[group.id[[j]]] <- sqrt(group.length[j])
    }
    len_wt <- sqrt_wt^2
    one_wt <- rep(1, p)

    # get Group SLOPE solutions with different choices of weights
    sqrt_wt_fit <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=group,
                                                    wt=sqrt_wt, lambda=lambda.mean,
                                                    verbose=FALSE)
    len_wt_fit <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=group,
                                                   wt=len_wt, lambda=lambda.mean,
                                                   verbose=FALSE)
    one_wt_fit <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=group,
                                                   wt=one_wt, lambda=lambda.mean,
                                                   verbose=FALSE)

    # store the sizes of the selected groups

    avail.length <- unique(group.length)
    selected.group.length <- rep(0, 3*length(avail.length))
    names(selected.group.length) <- c(paste0("sqrt_wt_", avail.length),
                                      paste0("len_wt_", avail.length),
                                      paste0("one_wt_", avail.length))

    for (j in 1:length(group.id)) {
      length_j <- length(group.id[[j]])
      if (sum(abs(sqrt_wt_fit$x[group.id[[j]]])) > 0) {
        selected.group.length[paste0("sqrt_wt_", length_j)] <-
          selected.group.length[paste0("sqrt_wt_", length_j)] + 1
      }
      if (sum(abs(len_wt_fit$x[group.id[[j]]])) > 0) {
        selected.group.length[paste0("len_wt_", length_j)] <-
          selected.group.length[paste0("len_wt_", length_j)] + 1
      }
      if (sum(abs(one_wt_fit$x[group.id[[j]]])) > 0) {
        selected.group.length[paste0("one_wt_", length_j)] <-
          selected.group.length[paste0("one_wt_", length_j)] + 1
      }
    }

    selected.group.length.fraction <-
      selected.group.length / max(1, sum(selected.group.length))

    # return the results
    selected.group.length.fraction
  }

  cat("done\n")
}

# collect results in a data frame
results <- data.frame() %>% tbl_df
for(k in 1:length(n.relevant)) {
  results <- bind_rows(results, tbl_df(parallel.results[[k]]))
}

# save the results
save(results, file = "../RData/figure_3_simulation_results.RData")
