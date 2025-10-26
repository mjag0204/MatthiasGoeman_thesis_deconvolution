install.packages("BiocManager")
BiocManager::install("deconvR")
install.packages("MASS")
BiocManager::install("EpiDISH")

library(deconvR)
library(nnls)
library(MASS)
library(EpiDISH)


##-------------
## READ DATA
##-------------

bulk0 <- read.csv("Pseudobulk0.csv", row.names = 1, check.names = FALSE)
bulk1 <- read.csv("Pseudobulk1.csv", row.names = 1, check.names = FALSE)
bulk2 <- read.csv("Psuedobulk2.csv", row.names = 1, check.names = FALSE)
bulk3 <- read.csv("Pseudobulk3.csv", row.names = 1, check.names = FALSE)
bulk4 <- read.csv("Pseudobulk4.csv", row.names = 1, check.names = FALSE)
bulk5 <- read.csv("Pseudobulk5.csv", row.names = 1, check.names = FALSE)
bulk6 <- read.csv("Pseudobulk6.csv", row.names = 1, check.names = FALSE)

refdata <- read.csv("refmatrix.csv", row.names = 1, check.names = FALSE)


##--------------------------------------
## convert the data to correct matrices
##--------------------------------------

  #convert to matrix or data.frame
  #-------------------------------

    refdata <- as.matrix(refdata)
    bulk0 <- as.matrix(bulk0)
    bulk1 <- as.matrix(bulk1)
    bulk2 <- as.matrix(bulk2)
    bulk3 <- as.matrix(bulk3)
    bulk4 <- as.matrix(bulk4)
    bulk5 <- as.matrix(bulk5)
    bulk6 <- as.matrix(bulk6)
    
    #refdata <- as.data.frame(refdata)
    #bulk0 <- as.data.frame(bulk0)

  #add common features to get the shared CpG IDs
  #---------------------------------------------

# check if the bulk data and the reference data have the same
# amount of CpG IDs as row names
# common features gives the rowID's that appear in both datasets

    #common_features <- intersect(rownames(bulk0), rownames(refdata))
    #bulk0 <- bulk0[common_features, ]
    #refdata <- refdata[common_features, ]


    #refdata$Feature <- rownames(refdata)
    #bulk0$Feature <- rownames(bulk0)


    #refdata <- refdata[, c(ncol(refdata), 1:(ncol(refdata)-1))]
    #bulk0 <- bulk0[, c(ncol(bulk0), 1:(ncol(bulk0)-1))]


##----------------------------------
## run the NNLS-based deconvolution
##----------------------------------

results <- deconvolute(
  reference = refdata,
  bulk = bulk0,
  model = "nnls"
)


##-------------------------------------
## run manual NNLS-based deconvolution
##-------------------------------------

manual_deconvolve <- function(reference, bulk) {
  
  n_samples <- ncol(bulk)
  n_celltypes <- ncol(reference)
  
  # make the results matrix for later
  proportions <- matrix(0, nrow = n_celltypes, ncol = n_samples)
  rownames(proportions) <- colnames(reference)
  colnames(proportions) <- colnames(bulk)
  
  #solve the NNLS for each sample
  for (i in 1:n_samples) {
    
    bulk_sample <- bulk[, i]
    
    # a true/false method to remove the NA or infinite numbers
    # only the 'good' values are now stored under values
    values <- !is.na(bulk_sample) & !is.infinite(bulk_sample) & 
              apply(reference, 1, function(x) all(!is.na(x) & !is.infinite(x)))
    
    nnls_results <- nnls(reference[values, ], bulk_sample[values])
    
    proportions[, i] <- nnls_results$x
    
    proportions[, i] <- proportions[, i] / sum(proportions[, i])
  }
  
  return(proportions)
}


##----------------------------------
## run the deconvolution function
##----------------------------------

results0 <- manual_deconvolve(refdata, bulk0)
print("Cell type proportions bulk 0")
print(results)

results1 <- manual_deconvolve(refdata, bulk1)
print("Cell type proportions bulk 1")
print(results1)

results2 <- manual_deconvolve(refdata, bulk2)
print("Cell type proportions bulk 0")
print(results2)

results3 <- manual_deconvolve(refdata, bulk3)
print("Cell type proportions bulk 3")
print(results3)

results4 <- manual_deconvolve(refdata, bulk4)
print("Cell type proportions bulk 4")
print(results4)

results5 <- manual_deconvolve(refdata, bulk5)
print("Cell type proportions bulk 5")
print(results5)

results6 <- manual_deconvolve(refdata, bulk6)
print("Cell type proportions bulk 6")
print(results6)


##----------------------------------------------------------
## Run the manual coded Robust Partial Correlation Function
##----------------------------------------------------------

rpc <- function(bulk_sample, reference_matrix) {
  
  n_celltypes <- ncol(reference_matrix)
  
  cor_vector <- rep(0, n_celltypes)
  for i in 1:n_celltypes {
    cor_vector[i] <- cor(bulk_sample, reference_matrix[, i],
                         method = "spearman", use = "complete.obs")
  }
  
  ref_cor_matrix <- cor(reference_matrix, method = "spearman", 
                        use = "complete.obs")
  tryCatch({
    inv_cor <- ginv(ref_cor_matrix)
    props <- inv_cor %*% cor_vector
    
    props[props < 0] <- 0
  })
  
}