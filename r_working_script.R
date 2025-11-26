# This is the script to collect all the working codes:
#=====================================================

# In this script, 5 different methods for reference-based deconvolution
#   will be used, as well as 2 statistical metrics to compare the results
#   with the given true propotions

# The methods used are:
#   1) non-negative least squares
#   2) robust partial correlation
#   3) support vector regression
#   4) quadratic programming
#   5) constrained projection

# The metrics being used are:
#   a) R squared
#   b) Root mean squared error (RMSE)


              #=======================================
              #=======================================
              ## PART 1: REFERENCE BASED DECONVOLUTION
              #=======================================
              #=======================================


#=====================================================
# Install and download packages
#=====================================================

install.packages("BiocManager")
BiocManager::install("deconvR")
install.packages("MASS")
BiocManager::install("EpiDISH")

library(deconvR)
library(nnls)
library(MASS)
library(EpiDISH)



#=====================================================
# Read data
#=====================================================

bulk0 <- read.csv("Pseudobulk0.csv")
bulk1 <- read.csv("Pseudobulk1.csv")
bulk2 <- read.csv("Psuedobulk2.csv")
bulk3 <- read.csv("Pseudobulk3.csv")
bulk4 <- read.csv("Pseudobulk4.csv")
bulk5 <- read.csv("Pseudobulk5.csv")
bulk6 <- read.csv("Pseudobulk6.csv")

refdata <- read.csv("refmatrix.csv")

true_props2 <- read.csv("trueprops.csv")
true_props3 <- read.csv("trueprops3.csv")
true_props4 <- read.csv("trueprops4.csv")
true_props5 <- read.csv("trueprops5.csv")
true_props6 <- read.csv("trueprops6.csv")

#true_props0 <- as.data.frame(lapply(true_props0, function(x) as.numeric(as.character(x))))
#true_props0 <- t(true_props0)
#true_props3 <- as.data.frame(lapply(true_props3, function(x) as.numeric(as.character(x))))
#true_props3 <- t(true_props3)
#true_props4 <- as.data.frame(lapply(true_props4, function(x) as.numeric(as.character(x))))
#true_props4 <- t(true_props4)
#true_props5 <- as.data.frame(lapply(true_props5, function(x) as.numeric(as.character(x))))
#true_props5 <- t(true_props5)
#true_props6 <- as.data.frame(lapply(true_props6, function(x) as.numeric(as.character(x))))
#true_props6 <- t(true_props6)


#=====================================================
# Preprocessing data
#=====================================================



#=====================================================
#Statistical metrics
#=====================================================

# The first metric used is R-squared. This metric calculates the variance in the 
#   estimated proprtions that can be explained by the true proportions.
#   It identifies similar patterns. The drawback here is that it does not account for 
#   absolute correctness but rather relative accuracy.

# The second metric used is RMSE. This metric calculates the absolute error between
#   the estimated proportions and the true proportions. The squared error prevents 
#   cancellation of positive and negative errors and makes it easier to interpret
#   as the RMSE is between 0 and 1. RMSE of 0.1 means that the predictions are 10% off.

                #=====================
                # Metric a) R squared
                #=====================

calculate_r_squared <- function(true_props, results) {
  
  #remove the NA values
  valid_idx <- !is.na(true_props) & !is.na(results)
  true_props <- true_props[valid_idx]
  results <- results[valid_idx]
  
  if(length(true_props) == 0) return(NA)
  
  #calculate the sum of squared residuals
  ss_res <- sum((true_props - results)^2)
  
  #calculate the total sum of squares
  ss_tot <- sum((true_props - mean(true_props))^2)
  
  #calculate R^2
  if(ss_tot == 0) {
    r_squared <- NA
  } else {
    r_squared <- 1 - (ss_res / ss_tot)
  }
  
  return(r_squared)
}


                #================
                # Metric b) RMSE
                #================

calculate_rmse <- function(true_props, results) {
  
  true_props <- as.matrix(true_props)
  results <- as.matrix(results)
  
  squared_errors <- (true_props - results)^2
  mse <- mean(squared_errors, na.rm = TRUE)
  rmse <- sqrt(mse)
  return(rmse)
}



#=====================================================
#Deconvolution methods
#=====================================================

                #================
                # Method 1) NNLS
                #================

bulk_list <- list(bulk2, bulk3, bulk4, bulk5, bulk6)
true_list <- list(true_props2, true_props3, true_props4, true_props5, true_props6)
nnls_results <- list()

for (i in seq_along(bulk_list)) {
  nnls_results[[i]] <- deconvolute(
  reference = refdata,
  bulk = bulk_list[[i]],
  model = "nnls"
)
}

nnls_results

#nnls_b0 <- nnls_results[1]
#nnls_b3 <- nnls_results[4]
#nnls_b4 <- nnls_results[5]
#nnls_b5 <- nnls_results[6]
#nnls_b6 <- nnls_results[7]

#r_squared <- list()

#r_squared[1] <- calculate_r_squared(true_props0, nnls_b0[[1]]$proportions)
#r_squared[2] <- calculate_r_squared(true_props3, nnls_b3[[1]]$proportions)
#r_squared[3] <- calculate_r_squared(true_props4, nnls_b4[[1]]$proportions)
#r_squared[4] <- calculate_r_squared(true_props5, nnls_b5[[1]]$proportions)
#r_squared[5] <- calculate_r_squared(true_props6, nnls_b6[[1]]$proportions)

#r_squared

calculate_r_squared(true_props2, nnls_results[[1]]$proportions)
calculate_rmse(true_props2, nnls_results[[1]]$proportions)

calculate_r_squared(true_props3, nnls_results[[2]]$proportions)
calculate_rmse(true_props3, nnls_results[[2]]$proportions)

calculate_r_squared(true_props4, nnls_results[[3]]$proportions)
calculate_rmse(true_props4, nnls_results[[3]]$proportions)

calculate_r_squared(true_props5, nnls_results[[4]]$proportions)
calculate_rmse(true_props5, nnls_results[[4]]$proportions)

calculate_r_squared(true_props6, nnls_results[[5]]$proportions)
calculate_rmse(true_props6, nnls_results[[5]]$proportions)

true_props2 <- as.matrix(true_props2)
est_props2 <- nnls_results[[1]]$proportions

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for(i in 1:ncol(true_props2)){
  cell_type <- colnames(true_props2)[i]
  
  plot(true_props2[, i], est_props2[, i],
       xlab = "True Proportions",
       ylab = "Estimated Proportions",
       main = cell_type,
       pch = 19,
       col = "steelblue",
       xlim = c(0, 1),
       ylim = c(0, 1))
  
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  
  cor_val <- cor(true_props2[, i], est_props2[, i])
  text(0.1, 0.9, paste("r=", round(cor_val, 3)), cex = 0.9)
}



                #===============
                # Method 2) rpc
                #===============

refdata <- as.matrix(refdata)
bulk3 <- as.matrix(bulk3)

refdata <- refdata[, -1]
refdata <- as.matrix(refdata)

bulk3 <- bulk3[, -1]
bulk3 <- as.matrix(bulk3)

rpc_b3 <- epidish(beta.m = bulk3, ref.m = refdata, method = "RPC")


qr_rank <- qr(refdata)$rank
ncol_ref <- ncol(refdata)
cat("Rank:", qr_rank, "  Columns:", ncol_ref, "\n")

round(cor(refdata), 3)

# These are extreme high correlations, meaning that most of the cell-type profiles
#   are almost linear combinations of each other

                #===============
                # Method 3) svr
                #===============

svr_results <- list()

for (i in seq_along(bulk_list)) {
  svr_results[[i]] <- deconvolute(
  reference = refdata,
  bulk = bulk_list[[i]],
  model = "svr"
)
}

svr_results

svr_b2 <- svr_results[3]

calculate_r_squared(true_props2, svr_results[[1]]$proportions)
calculate_rmse(true_props2, svr_results[[1]]$proportions)

calculate_r_squared(true_props3, svr_results[[2]]$proportions)
calculate_rmse(true_props3, svr_results[[2]]$proportions)

calculate_r_squared(true_props4, svr_results[[3]]$proportions)
calculate_rmse(true_props4, svr_results[[3]]$proportions)

calculate_r_squared(true_props5, svr_results[[4]]$proportions)
calculate_rmse(true_props5, svr_results[[4]]$proportions)

calculate_r_squared(true_props6, svr_results[[5]]$proportions)
calculate_rmse(true_props6, svr_results[[5]]$proportions)


true_props2 <- as.matrix(true_props2)
est_props2 <- svr_results[[1]]$proportions

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for(i in 1:ncol(true_props2)){
  cell_type <- colnames(true_props2)[i]
  
  plot(true_props2[, i], est_props2[, i],
       xlab = "True Proportions",
       ylab = "Estimated Proportions",
       main = cell_type,
       pch = 19,
       col = "steelblue",
       xlim = c(0, 1),
       ylim = c(0, 1))
  
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  
  cor_val <- cor(true_props2[, i], est_props2[, i])
  text(0.1, 0.9, paste("r=", round(cor_val, 3)), cex = 0.9)
}


                #==============
                # Method 4) QP
                #==============

qp_results <- list()

for (i in seq_along(bulk_list)) {
  qp_results[[i]] <- deconvolute(
    reference = refdata,
    bulk = bulk_list[[i]],
    model = "qp"
  )
}

qp_results

calculate_r_squared(true_props2, qp_results[[1]]$proportions)
calculate_rmse(true_props2, qp_results[[1]]$proportions)

calculate_r_squared(true_props3, qp_results[[2]]$proportions)
calculate_rmse(true_props3, qp_results[[2]]$proportions)

calculate_r_squared(true_props4, qp_results[[3]]$proportions)
calculate_rmse(true_props4, qp_results[[3]]$proportions)

calculate_r_squared(true_props5, qp_results[[4]]$proportions)
calculate_rmse(true_props5, qp_results[[4]]$proportions)

calculate_r_squared(true_props6, qp_results[[5]]$proportions)
calculate_rmse(true_props6, qp_results[[5]]$proportions)

true_props2 <- as.matrix(true_props2)
est_props2 <- qp_results[[1]]$proportions

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for(i in 1:ncol(true_props2)){
  cell_type <- colnames(true_props2)[i]
  
  plot(true_props2[, i], est_props2[, i],
       xlab = "True Proportions",
       ylab = "Estimated Proportions",
       main = cell_type,
       pch = 19,
       col = "steelblue",
       xlim = c(0, 1),
       ylim = c(0, 1))
  
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  
  cor_val <- cor(true_props2[, i], est_props2[, i])
  text(0.1, 0.9, paste("r=", round(cor_val, 3)), cex = 0.9)
}


true_props3 <- as.matrix(true_props3)
est_props3 <- svr_results[[2]]$proportions

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for(i in 1:ncol(true_props3)){
  cell_type <- colnames(true_props3)[i]
  
  plot(true_props3[, i], est_props3[, i],
       xlab = "True Proportions",
       ylab = "Estimated Proportions",
       main = cell_type,
       pch = 19,
       col = "steelblue",
       xlim = c(0, 1),
       ylim = c(0, 1))
  
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  
  cor_val <- cor(true_props3[, i], est_props3[, i])
  text(0.1, 0.9, paste("r=", round(cor_val, 3)), cex = 0.9)
}


                #==============
                # Method 5) CP
                #==============

cp_bulk3 <- epidish(beta.m = bulk3, ref.m = refdata, method = 'CP')

# The quadratic term D = R^tR is not positive definite, it is too collinear 
#   to invert. This is the same issue as the quadratic programming:
#   the reference matrix has extremely correlated columns





              #=======================================
              #=======================================
              ## PART 1: REFERENCE FREE DECONVOLUTION
              #=======================================
              #=======================================


#=====================================================
# Install and download packages
#=====================================================

#install.packages("BiocManager")
#BiocManager::install(version = "3.22")
#BiocManager::install(c("BiocGenerics", "S4Vectors", "IRanges", "SummarizedExperiment"))
#BiocManager::install("MeDeCom")
#library(MeDeCom)

install.packages("devtools")
devtools::install_github("CompEpigen/MeDeCom")
library(MeDeCom)



#=====================================================
#Deconvolution methods
#=====================================================

bulk5 <- read.csv("Pseudobulk5.csv", row.names = 1)
bulk5 <- as.matrix(bulk5)
bulk5 <- apply(bulk5, 2, as.numeric)

Ks <- 6:8                  
lambdas <- c(0, 0.01, 0.1) 
itermax <- 100
nfolds <- 2
ninit <- 5

medecom_results <- MeDeCom::runMeDeCom(
  bulk5,
  Ks = Ks,
  lambdas = lambdas,
  Ninit = ninit,
  iter_max = itermax,
  Nfold = nfolds,
  NFolds = nfolds,
  cores = 1   
)

# The pseudo cell type profiles
lmc <- getLMCs(medecom_results)

# Estimated cell type proportions
props <- getProportions(medecom_results)

# Plot the proportions
plotProportions(medecom_results)


