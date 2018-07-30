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