##########################################
#
#  Introduction to R 
#  Option B: For users with some experience with R
#  A tutorial for STAMPS 2017
# 
#########################################
#
#  This tutorial is written by Amy Willis, July 2017
#  Contributions also made by [please feel free to contribute!]
#
#########################################
#
# If you downloaded this from the website, 
# be sure to rename as "STAMPS_Intro2R_OptionB.R", 
# R and R Studio will not interpret is as an R script with ".txt" suffix

#########################################
#
# Disclaimer: The level of hand-holding in this tutorial will be
# less than in Option A. There will be some expectation on you to
# attempt to solve your issues. By all means request help from the TAs,
# but in order to obtain mastery of any language it is important to practice 
# correcting yourself, and now is a great opportunity to practice :)
# 
#########################################


###############################################
# 
# Get some context 
#
###############################################
# Set your working directory (refer to Option A if you're not sure how)
setwd("~/Documents/MyProjectDirectory") 

# If you already worked through Option A, feel free to skip until the next section

# In Option A, we loaded in some covariate and OTU data and briefly
# made some plots and looked at some summaries
covariates <- read.table("FWS_covariates.txt", header=TRUE, sep="\t", as.is=TRUE)
abundances <- read.table("FWS_OTUs.txt", header=TRUE, row.names = 1, sep="\t", as.is=TRUE)

# take a moment to have a look at the data here

###############################################
# 
# ggplot
#
###############################################

# Last night you learned how to install a package from CRAN
# We're now going to use the package ggplot2

# If this doesn't work for you, look at the Installing R website
library(ggplot2)

# check that the following produces a plot
ggplot(data.frame("x"=c(1,2,3), "y"=c(4,4.6, 5)))+geom_point(aes(x,y))
# and don't worry about what this means for now

# base graphics in R, such as plot() and boxplot(), are not very pretty
# enter: ggplot! (everyone calls the package ggplot2 "ggplot") 

# ggplot makes simple things difficult, difficult things simple, and
# everything beautiful

# ggplot has a different structure to base graphics
# here is how you make a ggplot
abundances <- data.frame(abundances)
ggplot(abundances) 

# nothing happened! why?
# we need to add components, such as a scatterplot or a histogram, 
# to build up a plot
ggplot(abundances, aes(x = JPA_Jan)) +
  geom_histogram()

# this says that "abundances" is the data frame that we want to plot, 
# JPA_Jan is the only variable that we care about
# and that we want a histogram

# obviously there are lots of taxa that were not observed in this sample, hence the zeros
# let's restrict our attention to taxa that we did see
# Subsetting 2 different variables is done like this
ggplot(subset(abundances, JPA_Jan > 0), aes(JPA_Jan)) + 
  geom_histogram()

# the first argument is the data frame of interest
# aes() specifies which column we want
# and the geom_histogram() says that we want a histogram

# The range of this histogram is huge, so we want to focus on the distribution of
# OTUs that were observed once or more but less than 60 times, we would do
ggplot(subset(abundances, JPA_Jan > 0 & JPA_Jan < 60), aes(JPA_Jan)) + 
  geom_histogram()

# a scatterplot shows how correlated the abundances are between january and february
ggplot(subset(abundances, JPA_Jan > 0 & JPA_Feb > 0), 
       aes(JPA_Jan, JPA_Feb)) + 
  geom_point()
# again, ggplot() specifies what we want plotted
# geom_point() specifies that we want a scatterplot

ggplot(subset(abundances, JPA_Jan > 0 & JPA_Feb > 0), 
       aes(JPA_Jan, JPA_Feb)) + 
  geom_point() + 
  theme_classic() # changes the styling
  
ggplot(abundances, aes(JPA_Jan, JPA_Feb)) + 
  geom_point() + 
  labs(title="My first scatterplot") + # adds a title
  xlab("January OTU counts") + # adds labels 
  ylab("February OTU counts")
  
# one tedious thing about ggplot is that all of the information that
# you want to plot must be in the same data frame
# So incorporating Season (from the data frame "covariates") 
# to colour the points our scatterplot is not trivial

# We need to create a new data frame that contains all this information
# together... 

# Exercise: Use the internet to learn how to change the colour
# of the points. Create a new data frame that contains a
# logical [TRUE/FALSE] variable indicating if the January abundance
# is higher than 4000. Use this variable to differentially colour the single
# observation that fits this category

# Don't look at the solution until you've given it a good attempt 
# (this is equally good practice for troubleshooting your problems)
# but the solution is available in STAMPS_Intro2R_OptionB_solutions.R

###############################################
# 
# loops
#
###############################################

# R is a vectorised language, so for a function that takes 1 argument, 
# you can apply across a vector
abundances[, 1] + 5
exp(abundances[, 1])

# what about for more complex functions?

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

###############################################
# 
# The apply family
#
###############################################

# Unfortunately, for loops, while intuitive and broadly useful, 
# are relatively slow in R

# apply() is a much better way of doing this!

fast_total_reads <- apply(abundances, 2, sum)
# what does the second argument 2 do? Compare to second argument of 1
# Find out with ?apply

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

# Note that you can construct functions inline:
all(apply(abundances, 2, function(x) x/sum(x)) == relative_abundances)

# Save your work!

# Congratulations! You are at a level with your R understanding
# that will get you through most of STAMPS! If you feel so inclined, 
# have a look through Option C, but at this point, you probably deserve
# a cold beverage and a break!