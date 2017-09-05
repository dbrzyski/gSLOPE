#--------------------------------------------------------------------------
# This code is based on the simulation study of Figure 4 in D. Brzyski, 
# W. Su, M. Bogdan (2015), "Group SLOPE - adaptive selection of 
# groups of predictors" (http://arxiv.org/abs/1511.09078)
#--------------------------------------------------------------------------

library(grpSLOPE)
library(dplyr)
library(tidyr)
library(grpreg)

# Adjust the number of cores to your particular system
library(doParallel)
registerDoParallel(cores = as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")))

#--- Set up global parameters for the simulation

set.seed(20170828)

# auxilliary function to get (group-wise) FDP and power from one grpSLOPE solution object
get_FDP_and_power <- function(result, true.relevant){
  truepos <- length(intersect(result$selected, true.relevant))
  falsepos <- length(result$selected) - truepos
  FDP <- falsepos / max(1, length(result$selected))
  pow <- truepos / length(true.relevant)
  return(c("FDP" = FDP, "pow" = pow))
}

# set the grouping structure
n.group <- 1000
group.length <- rbinom(n.group, 1000, 0.008)
group.length[group.length == 0] <- 1 # avoid having groups of size 0
group <- c()
for (i in 1:n.group) {
  group <- c(group, rep(i, group.length[i]))
}
group.id <- getGroupID(group)

# design matrix dimensions
n <- 5000
p <- length(group)

# determine signal strength, such as used in Figure 1 in Brzyski et. al. (2015)
Bfun <- function(l) {
  sqrt(4*log(n.group) / (1 - n.group^(-2/l)) - l)
}
signal.strength <- sum(Bfun(group.length)) / n.group

# considered numbers of truly relevant groups
n.relevant <- floor(seq(3, 60, length=7))

# how many times the simulation is repeated
n.replications <- 200

#--- Simulation main loop

# run the simulations n.replications times at each level of n.relevant,
# and with each combination of lambda="max", lambda="mean", fdr=0.1, fdr=0.05
parallel.results <- vector(mode="list")
for (k in 1:length(n.relevant)) {
  cat(paste("sparsity level", k, "started"))

  parallel.results[[k]] <- foreach(i=1:n.replications, .combine = rbind) %dopar% {
    cat(".")

    # generate and normalize the model matrix
    X <- matrix(rnorm(n*p), n, p)
    X <- scale(X, center=TRUE, scale=FALSE)
    X <- apply(X, 2, function(x) x/sqrt(sum(x^2)) )

    # generate coeffient vector, pick relevant groups at random
    b <- rep(0, p)
    n.signif <- n.relevant[k]
    ind.relevant <- sample(1:n.group, n.signif)
    for (j in ind.relevant) {
      signals <- runif(group.length[j]) + 0.1
      X1 <- as.matrix(X[ , group.id[[j]]]) %*% signals
      b[group.id[[j]]] <- (signal.strength / norm(as.matrix(X1), "f")) * signals
    }

    # generate the response vector
    y <- X %*% b + rnorm(n)

    # get Group SLOPE solutions with different lambda and fdr values
    # (this has the same gFDR controlling properties with orthogonalize=FALSE too)
    lambda.1 <- grpSLOPE(X = X, y = y, group = group, fdr = 0.1)
    lambda.05 <- grpSLOPE(X = X, y = y, group = group, fdr = 0.05)

    # get CV Group LASSO solution
    cvgrplasso <- cv.grpreg(X, y, group, penalty = "grLasso", family = "gaussian")

    # get the FDPs and powers of the grpSLOPE solutions
    true.relevant <- names(group.id)[ind.relevant]
    FDPs.and.powers <- rep(NA, 8)
    names(FDPs.and.powers) <- c("n.relevant", "replication",
                                "lambda.1.FDP", "lambda.1.power",
                                "lambda.05.FDP", "lambda.05.power",
                                "grpLASSO.FDP", "grpLASSO.power")
    FDPs.and.powers["n.relevant"] <- length(true.relevant)
    FDPs.and.powers["replication"] <- i
    FDPs.and.powers[c("lambda.1.FDP", "lambda.1.power")] <-
      get_FDP_and_power(lambda.1, true.relevant)
    FDPs.and.powers[c("lambda.05.FDP", "lambda.05.power")] <-
      get_FDP_and_power(lambda.05, true.relevant)

    # get FDP and power for the CV Group LASSO solution
    cvgrplasso.coef <- coef(cvgrplasso)[-1]
    total.discoveries.cvgrplasso <- 0
    true.discoveries.cvgrplasso <- 0
    false.discoveries.cvgrplasso <- 0
    for (i in 1:n.group) {
      if(sqrt(sum(cvgrplasso.coef[group.id[[i]]]^2)) > 1e-6){
        total.discoveries.cvgrplasso <- total.discoveries.cvgrplasso + 1
        if(names(group.id)[i] %in% true.relevant) {
          true.discoveries.cvgrplasso <- true.discoveries.cvgrplasso + 1
        } else {
          false.discoveries.cvgrplasso <- false.discoveries.cvgrplasso + 1
        }
      }
    }
    FDPs.and.powers[c("grpLASSO.FDP", "grpLASSO.power")] <-
      c(false.discoveries.cvgrplasso/max(total.discoveries.cvgrplasso, 1),
        true.discoveries.cvgrplasso/length(true.relevant))

    FDPs.and.powers
  }

  cat("done\n")
}

#--- Collect and summarize simulation results

# collect results in a data frame
results <- data_frame()
for(k in 1:length(n.relevant)) {
  results <- bind_rows(results, tbl_df(parallel.results[[k]]))
}

save(results, parallel.results, n.replications, n.relevant, p, group.id,
     file = "../RData/figure_4_simulation_results.RData")
