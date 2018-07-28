###############################################
# 
# additional R-code from previous R-tutorials B & C
#
###############################################

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


##########################################
#
#  Introduction to R
#  Option C: For users with a lot of experience with R
#  A tutorial for STAMPS 2017
#
#########################################
#
#  This tutorial was written in July 2017
#  A collaborative effort by 
#  Julia Fukuyama, Lan Huong Nguyen, 
#  Kris Sankaran, and Amy Willis
# 
#  [please feel free to contribute!]
#
#########################################
#
# Disclaimer: TA support is primarily for Options A and B.
# Please make every effort to resolve any difficulties that you encounter
# on your own. That said, we hope that you learn something new from this
# tutorial!
#
#########################################
# Set your working directory and load data
setwd("~/Documents/MyProjectDirectory")
covariates <- read.table("FWS_covariates.txt", header=TRUE, sep="\t", as.is=TRUE)
abundances <- read.table("FWS_OTUs.txt", header=TRUE, row.names = 1, sep="\t", as.is=TRUE)

# take a moment to have a look at the data here

###############################################
#
# more ggplot, reshape, cast, melt
# By Kris Sankaran
#
###############################################
library(ggplot2)
library(reshape2) # melt and cast
library(dplyr) # filter, group_by, and mutate

# let's try to plot the abundance of the dominant 2
# taxa against each other, and colour the points
# based on Season

# We need to reshape the data so that what will go
# in the x and y for the scatterplot are their own
# columns. We need to *melt* the data. This takes it
# from "wide" (one sample per column) to tall
# (different samples stacked on top of one another).

# First let's retrieve the taxonomic information,
# using strsplit to split the rownames according to
# the ";" character
taxa <- strsplit(rownames(abundances), split = ";")
taxa <- do.call(rbind, taxa) # bind the list into a matrix
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# We can use the %>% operation so that we don't need
# to keep track of temporary variables. Melt converts
# the data from wide to tall
melted_abundances <- data.frame(taxa, abundances) %>%
  melt(
    id.vars = colnames(taxa),
    variable.name = "SampleName",
    value.name = "abundance"
  )

# We can incorporate the season information using
# left_join, and get the taxa abundances using
# group_by_at and summarise
melted_abundances <- melted_abundances %>%
  group_by_at(colnames(taxa)) %>%
  mutate(total = sum(abundance)) %>%
  left_join(covariates) %>%
  ungroup()

# To get the most abundant taxa, we filter down
# to the top two totals. Then, we cast the data.frame
# so that there is one column for what will become
# each axis in the plot
top_abundances <- melted_abundances %>%
  filter(total >= sort(rowSums(abundances), decreasing = TRUE)[2]) %>%
  dcast(SampleName + Season ~ Genus, value.var = "abundance")

# Now we can make the plot! Notice that most of
# the work was in defining data with the appropriate
# input format
ggplot(top_abundances) +
  geom_point(aes(x = Burkholderia, y = Ralstonia, col = Season))