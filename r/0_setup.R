###############################################
# 
# Setup
#
###############################################

# First, set your working directory to wherever you want to 

utils::download.file("LINK", "R_tutorial.zip")
utils::unzip("R_tutorial.zip")

# We will move into the folder to make sure that everyone is working within the same directory
setwd("./STAMPS2018")

# Before moving on, check to make sure you have 7 ".R" files and 2 ".txt" files
list.files()

# If you have that, then load in the data we'll be using
covariates <- read.csv("FWS_covariates.txt", sep = "\t")
abundances <- read.csv("FWS_OTUs.txt", sep = "\t", row.names = 1, header = T)

# Let's install 
install.packages("dplyr")
install.packages("magrittr")
install.packages("devtools")
devtools::install_github("adw96/breakaway")
install.packages("ggplot2")
install.packages("phyloseq")
install.packages("parallel")
install.packages("foreach")
install.packages("doParallel")


