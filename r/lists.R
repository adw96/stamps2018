###############################################
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

# you can also make a list of lists!
all_samples <- list()
all_samples$sample1 <- sample1
all_samples$sample2 <- list()
all_samples$sample2$name <- covariates$SampleName[2]
all_samples

# Exercise: Use a loop to loop through all of the samples, 
# creating a list of the information corresponding to that sample.
# (name, location, month, relative abundance table, abundance table)

# Hint:
all_samples[[1]]
all_samples[[2]]
# Double brackets can also be used to refer to the elements of a list, 
# instead of naming them individually
names(all_samples)


# Bonus: do this with a apply loop (maybe after finishing the function
# writing section)