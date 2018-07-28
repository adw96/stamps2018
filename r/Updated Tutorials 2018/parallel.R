###############################################
#
# parallelisation
# By Lan Huong Nguyen
#
###############################################

# In Option B we calculated the relative abundance of the taxa
# using for loops
# Now try to perform the same but in parallel

# If you haven't done it before, install the 'foreach' and 'doParallel' packages
# install.packages("foreach")
# install.packages("doParallel")

# Load the packages 
library(doParallel)

# Set the number of cores
no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)

# Use foreach function to compute the relative abundances
relative_abundances <- 
  foreach(current_index = 1:(dim(abundances)[2]), 
          .combine = cbind)  %dopar% {
            number_reads <- sum(abundances[, current_index]) 
            abundances[, current_index]/number_reads
          }


# Check if results match the standard calculations
all(relative_abundances == scale(abundances, center = F, scale = colSums(abundances)))

# Note that:'.combine' argument tells foreach function how to combine results
# from multiple threads. '.combine' can be specified as a function or 
# a character string e.g. c, cbind, rbind, list, '+', '*'

# Also observe that variables from the local environment are by default available
# for each thread without specification. This is not true of the variables
# in the parent environment, e.g. the following gives an error:

test <- function () {
  foreach(current_index = 1:(dim(abundances)[2]), 
          .combine = cbind)  %dopar% {
            number_reads <- sum(abundances[, current_index]) 
            abundances[, current_index]/number_reads
          }
}
relative_abundances <- test()

# To export specific variables to the cluster use the 'export' argument
test <- function () {
  foreach(current_index = 1:(dim(abundances)[2]), 
          .combine = cbind,
          .export = "abundances")  %dopar% {
            number_reads <- sum(abundances[, current_index]) 
            abundances[, current_index]/number_reads
          }
}
relative_abundances <- test()

# Check if results match the standard calculations
all(relative_abundances == scale(abundances, center = F, scale = colSums(abundances)))

# When you are done working in paralel close the cluster so that resources 
# are returned to the operating system
stopCluster(cl)


# In R  lapply() is a faster and cleaner way to perform for loop tasks. 
# You can speed things up even further by using parallel equivalent of 
# lapply() and  sapply(), parLapply() and parSapply() respectively. These
# functions are available in the package, 'parallel'.

# install.packages("parallel")

library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

# Export variable 'abundances' to the cluster
clusterExport(cl, "abundances")

# If you are using some specific packages inside parLapply, you need to load 
# them through  clusterEvalQ(), e.g. clusterEvalQ(cl, library(phyloseq))

# Use parLapply() function
relative_abundances <- 
  parLapply(cl, 
            1:(dim(abundances)[2]),
            function(current_index) {
              print(class(current_index))
              number_reads <- sum(abundances[, current_index]) 
              abundances[, current_index]/number_reads
            })

# Close the cluster
stopCluster(cl)

# Currently, 'relative_abundances' is a list. To converti it to a matrix do
# the following:
relative_abundances <- do.call("cbind", relative_abundances)

# Check if results match the standard calculations
all(relative_abundances == scale(abundances, center = F, scale = colSums(abundances)))

# For more details on parallelisms check the following link:
# http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
