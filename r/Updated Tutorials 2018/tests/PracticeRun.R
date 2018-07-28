setwd("/Users/paulinetrinh/Documents/GitHub/stamps2018/r/Updated Tutorials 2018")
covariates <- read.table("/Users/paulinetrinh/stamps2018/r/FWS_covariates.txt", header=TRUE, sep="\t", as.is=TRUE)
abundances <- read.table("/Users/paulinetrinh/stamps2018/r/FWS_OTUs.txt", header=TRUE, row.names = 1, sep="\t", as.is=TRUE)

abundances[, 1] + 5
exp(abundances[, 1])


# what about for more complex functions?

###############################
# LOOPS
###############################

# In Option A we calculated the relative abundance of the taxa
# using scale()
# This was a good way to do it, are we going to try to do the
# same thing with a loop

relative_abundances <- abundances # create a matrix of the same size as abundances
for (current_index in 1:(dim(abundances)[2])) { # cycling through the columns...
  number_reads <- sum(abundances[, current_index]) # find the column total
  relative_abundance <- abundances[, current_index]/number_reads # use it to find relative abundance
  relative_abundances[, current_index] <- relative_abundance # save the current column into our new matrix
}

# we can track that this matches our original calculation with
all(relative_abundances == scale(abundances, center = F, scale = colSums(abundances)))

# what does all() do?
?all

# Exercise: use a for loop to find which genus was observed most frequently
# when combining all samples

#create a matrix of 1 column to hold our total number of sequences
mostcommongenus <- abundances[,1, drop=FALSE] #create a matrix of only 1 columns

for (i in 1:(dim(abundances)[1])) {
  genus_reads <- sum(abundances[i,]) #find the row total counts 
  mostcommongenus[i,] <- genus_reads #insert each row count into our 
}

rownames(mostcommongenus)[which.max(mostcommongenus[,1])]



############################
# apply
############################
# Unfortunately, for loops, while intuitive and broadly useful, 
# are relatively slow in R

# apply() is a much better way of doing this!

fast_total_reads <- apply(abundances, 2, sum)
# what does the second argument 2 do? Compare to second argument of 1
# Find out with ?apply 
# 1 indicates rows, 2 indicates columns so sum the columns is what's happening with fast_total_reads

# sweep() is related to apply()
sweep(abundances, 2, fast_total_reads, "/")
# "go down the columns of abundances, applying "/" (division) to each column,
# with the second argument of the division being the element of
# fast_total_reads corresponding to the column number of abundances
all(sweep(abundances, 2, fast_total_reads, "/") == relative_abundances)


# careful! dividing a matrix by a vector divides across the rows, 
# not the columns!
all(abundances/fast_total_reads == relative_abundances)

# mapply(), tapply(), lapply() are all variations on apply for different 
# situations

# In writing this tutorial, I just learnt about prop.table,
# which turns a matrix into its relative abundances across a fixed margin.
# Cool!
all(prop.table(as.matrix(abundances), margin = 2) == relative_abundances)


##############################################
# 
# Lists
#
###############################################

# Lists are great! Lists allow you to combine many different types of data
# of many different sizes into groups

# We're going to make a list corresponding to all of the information
# about the first sample

sample1 <- list()
sample1$name <- covariates$SampleName[1]
sample1$season <- covariates$Season[1]
sample1$place <- covariates$Location[1]
sample1
sample1$counts <- abundances[,1]

# to see all the information in a list, go
sample1
# and to see a particular element of the list, 
# place the name of that element after the sign $
sample1$name
# you can also call the nth element of the list 
sample1[1]

# you can also make a list of lists!
all_samples <- list()
all_samples$sample1 <- sample1
all_samples$sample2 <- list()
all_samples$sample2$name <- covariates$SampleName[2]
all_samples[[1]]

# Exercise: Use a loop to loop through all of the samples, 
# creating a list of the information corresponding to that sample.
# (name, location, month, relative abundance table, abundance table)
for(i in 1:3){
  assign(paste("Sample", i, sep = ""), i)    
}
ls()

name <- paste("Sample", i, sep="") 



all_samples[[2]][[1]] # will call the first list

alist <- list(2)
alist <- list(3)
all_samples <- list()
for (i in 1:(dim(covariates)[1])) { # cycle go through each row 
  for (j in 1:(dim(covariates)[2])){ # within each row go through each column
    anotherlist[j] <- list(covariates[i,j])
  }
  all_samples[i] <- list(anotherlist)
} 


l1<- as.list(c(1,2,3,4,5))

# Hint:
all_samples[[1]]
all_samples[[2]]
# Double brackets can also be used to refer to the elements of a list, 
# instead of naming them individually
names(relative_abundances)


# Bonus: do this with a apply loop (maybe after finishing the function
# writing section)




###############################################
# 
# Function writing
#
###############################################

# writing functions is easy in R!
my_first_function <- function() {
  cat("Hello World!")
}
my_first_function()
# functions can have no arguments, or several arguments:
my_second_function <- function(times) {
  counter <- 1
  while (counter <= times) { # this is a different type of loop
    cat("Hello World!\n")
    counter <- counter + 1
  }
}
my_second_function(3)

# If you are dealing with one OTU table, turning it into
# relative abundance table is pretty easy. But what if you have 
# to do it for many OTU tables?

# Exercise: Write a function that takes in an OTU table
# and returns the relative abundance table

#My Answer
otu_function <- function(table) {
  apply(abundances, 2, function(x) x/sum(x))
}
relativeabundance_fxn <- otu_function(abundances) 
  
# Note that you can construct functions inline:
all(apply(abundances, 2, function(x) x/sum(x)) == relative_abundances)

# Save your work!

# Congratulations! You are at a level with your R understanding
# that will get you through most of STAMPS! If you feel so inclined, 
# have a look through Option C, but at this point, you probably deserve
# a cold beverage and a break!



###############################################
#
# source and system
#
###############################################

# Exercise: Write a script to separate out the different 
# taxonomy data in shell
# and run it through R

# Hint: Use source() to run an R script in this session
# Also works for URLs!
?source

# Hint: Use system() to run a system command. 
# In the alpha diversity lab we will see an example
# of this (running CatchAll via the command line)
?system
