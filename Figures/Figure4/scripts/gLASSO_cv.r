###############################################################
###############################################################

library(gglasso)
m = 1000;
n = 5000;


path ="path to \\Figures\\Figure4\\scripts\\"
save_path ="path to the directory when results should be saved"

setwd(path)
Lgths  <- read.table("Lgths.txt");
Lgths  <- as.numeric(Lgths);
Lgths  <- as.matrix(Lgths);
I  <- read.table("I.txt");
I  <- as.numeric(I);
I  <- as.matrix(I);

p    <- length(I);
K    <- c(3, 10, 20, 30, 40, 50, 60);
iter <- 200;
eps  <-1e-6

# center with 'apply()'
center_apply <- function(x) {
  apply(x, 2, function(y) y - mean(y));
}
# normalize with 'apply()'
normalize_apply <- function(x) {
  apply(x, 2, function(y) y/sqrt(sum(y*y)) );
}

FDR        = rep(c(0), times= length(K));
POWER      = rep(c(0), times= length(K));
stdFDR     = rep(c(0), times= length(K));
stdPOWER   = rep(c(0), times= length(K));
LNGTHs     = rep(c(0), times= length(K));
LNGTHs_t   = rep(c(0), times= length(K));
LAMB       = rep(c(0), times= length(K));

SigStrg = sqrt( 4*log(m)/(1-m^(-2/Lgths)) - Lgths );
SigStrg = mean(SigStrg);

###############################################################
###############################################################

for (jj in 1:length(K)) {

k            = K[jj];
power    = rep(c(0), times= iter);
fdr      = rep(c(0), times= iter);
lamb     = rep(c(0), times= iter);


for (ii in 1:iter) {
  print(ii)
  X <- rnorm(n*p);
  X <- matrix(X,n,p);
  X <- center_apply(X);
  X <- normalize_apply(X);
  
  for (dd in 1:m) {
    Idd     <- which(I==dd)
    A       <- X[,Idd]
    Q       <- qr.Q(qr(A))
    X[,Idd] <- Q
  }
  
  ## true beta generating
  beta   <- rep(c(0), times= p);
  beta   <- as.matrix(beta,p,1);
  group_supp <- sample(1:m, k, replace = FALSE)
  z <- rnorm(n);
  z <- matrix(z,n,1);
  
  for (gg in 1:k) {
    ggroup       <- which(I==group_supp[gg])
    ggsize       <- length(ggroup)
    ggsignals    <- runif(ggsize)+0.1
    beta[ggroup] <- ggsignals*SigStrg/sqrt(sum(ggsignals^2));
  }
  
  ## observations generating
  y = X%*%beta+z;
  
  ## solution of gLASSO with cross validation
  cv=cv.gglasso(X, y, I, lambda = NULL, pred.loss = "L2", nfolds = 5, intercept=FALSE);
  gLASSO = coef(cv$gglasso.fit, s = cv$lambda.min);
  gLASSO = gLASSO[2:(p+1)];
  gLASSO_I = rep(c(0), times= m);
  gLASSO_I   <- as.matrix(gLASSO_I,m,1);
  for (ggg in c(1:m)) {
    ggroup       <- which(I==ggg)
    gLASSO_I[ggg] <-sum(abs(gLASSO[ggroup]))
  }
  all_discoveries   = length(which(gLASSO_I>eps));
  true_discoveries  = length(intersect(which(gLASSO_I>eps), group_supp));
  false_discoveries = all_discoveries - true_discoveries;
  
  ## statistics
  fdr[ii]   = false_discoveries/max(1,all_discoveries);
  power[ii] = true_discoveries/k;
  if (ii == 1){
    lngths <- Lgths[which(gLASSO_I>eps)]
  } else {
    lngths <- c(lngths, Lgths[which(gLASSO_I>eps)])
  }
  if (ii == 1){
    lngths_t <- Lgths[intersect(which(gLASSO_I>eps), group_supp)]
  } else {
    lngths_t <- c(lngths_t,Lgths[intersect(which(gLASSO_I>eps), group_supp)])
  }
  
  lamb[ii] <- cv$lambda.min

}

FDR[jj]      = mean(fdr);
POWER[jj]    = mean(power);
stdFDR[jj]   = sd(fdr);
stdPOWER[jj] = sd(power);
LNGTHs[jj]   = mean(lngths);
LNGTHs_t[jj] = mean(lngths_t);
LAMB[jj]     = mean(lamb);

## saving
setwd(save_path)
write(FDR, "gFDR_gLASSO_CV.txt", sep="\t")
write(POWER, "POWER_gLASSO_CV.txt", sep="\t")
write(stdFDR, "stdgFDR_gLASSO_CV.txt", sep="\t")
write(stdPOWER, "stdPOWER_gLASSO_CV.txt", sep="\t")
write(LNGTHs, "LNGTHS_gLASSO_CV.txt", sep="\t")
write(LNGTHs_t, "LNGTHS_gLASSO_CV_t.txt", sep="\t")
write(LAMB, "LAMB_CV.txt", sep="\t")
}



