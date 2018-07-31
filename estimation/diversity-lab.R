## diversity-lab.R
## A script to introduce you to alpha diversity estimation and comparison

## Lab author: Amy Willis

## breakaway authors: Amy Willis, Kathryn Barger, John Bunge, 
##                        Bryan Martin, 2012+

## Start by checking what directory you are in with
getwd()
## Set it to your local copy: setwd([your directory name here])


library(tidyverse)

## Let's use data distributed from phyloseq
library(phyloseq)
data("GlobalPatterns")
GlobalPatterns

## Have a look through GlobalPatterns samples
GlobalPatterns %>% sample_data


# Let's just look at the water samples of GlobalPatterns
# To speed things up, let's just aggregate taxa to the 
# order level. So we're going to estimate *order* level 
# diversity 
# This might take a while, but it will speed things up later
water <- GlobalPatterns %>%
  subset_samples(SampleType %in% c("Freshwater", 
                                   "Freshwater (creek)", 
                                   "Ocean", 
                                   "Sediment (estuary)")) %>%
  tax_glom("Order")
water



# phyloseq has some inbuilt tools for exploring alpha
# diversity, but they're not great. They understate richness,
# understate uncertainty, and don't allow hypothesis testing

# Enter: breakaway :D
# breakaway was specifically designed for these tasks

devtools::install_github("adw96/breakaway")
library(breakaway)
# If this doesn't work, run install.packages(devtools) and try again
# If this doesn't work, read the error message carefully 
# and try to debug it yourself (this is what you have to do at home -- so
# it's good to practice here)
# If this doesn't work, call a TA over


## Let's look at the observed richness of the water samples
observed_c <- sample_richness(water)
summary(observed_c)
plot(observed_c, water, color = "SampleType")

# Hmmmm, but what if we observed these samples at different depth?
# Depth may be confounded with observed richness. 
# Let's practice ggplot and phyloseq to look at this

data.frame("observed_richness" = (observed_c %>% summary)$estimate,
           "depth" = phyloseq::sample_sums(water), # Easter egg! Phyloseq's function to get depth
           "type" = water %>% sample_data %>% get_variable("SampleType")) %>%
  ggplot(aes(x = depth, y = observed_richness, color = type)) +
  geom_point()

# So our lowest depth samples (Sediment) had lower richness
# compared to Freshwater (creek)!
# Not surprising!

# Recall that the best way to adjust for this is by estimating the 
# number of missing species using a species richness estimate.

ba <- breakaway(water)
ba 
plot(ba, water, color = "SampleType")

## Cool! How do these estimates work?

## Let's look at just one sample
tr <- water %>% subset_samples(X.SampleID == "TRRsed1")
tr

# what's the structure of this dataset?
fc <- tr %>% otu_table %>% make_frequency_count_table
# this is the frequency count table for this dataset
fc %>% head(10)
# So there are 11 singletons (Orders observed once) here

# Quiz: how many times was the most common Order observed?
# (Hint: what does tail() do?)
# (A: 8863 times)

# Let's fit the breakaway to this sample
ba_tr <- breakaway(fc)
ba_tr
# this is an alpha diversity estimate -- a special class for
# alpha diversity estimates

# breakaway picks lots of models and chooses the best
# Which model did it pick?
ba_tr$model

# Kemp -- that's a non-mixed Poisson model. Cool!

# Kemp models work by fitting a legit probabilistic model to
# a transformation of the data. We can plot the transformation and the fit:
ba_tr %>% plot
# (Not super important -- but if you want to know more, check out
# the paper: Willis & Bunge (2015), Biometrics)

# breakaway is easy to use! You can run it on one frequency table, but
# in general you'll run it on a phyloseq object, like here:
ba <- breakaway(water)
ba 

# you can plot using our default plotting function...
plot(ba, water, color = "SampleType")

# ... or you can take the estimates and turn them into 
# data frames so you can plot them yourself
summary(ba) 
# these are in the same order as the phyloseq samples:
water %>% otu_table %>% sample_names

# Also, since the breakaway package implements lots of
# species richness estimates, you could choose a different one
# e.g., if you wanted something more stable
water %>%
  chao_bunge %>%
  plot(water, color = "SampleType")

# However, note that species richness is a challenging problem, 
# and that error bars will generally be large

# grrrrr don't use this one
water %>%
  chao1 %>%
  plot(water, color = "SampleType")
# (but note that these are the real error bars on this estimate)

# In many cases the absolute number of species
# isn't as interesting as comparing ecosystems. 
# Let's test the hypothesis that different types of water systems
# have the same microbial diversity
betta()

#### Let's move on to estimating other diversity indices

# Species richness counts all species equally
# However, if a species is rare you may think that
# it doesn't play a role in the community.
#  Another alpha diversity index, called the Shannon
# index, works similarly to species richness but it
# down weights the importance of  rare taxa

# Since rare taxa may be dubious, the Shannon index
# is very popular in microbial ecology

# For the reasons discussed in lecture, it's important to
# estimate the Shannon diversity using a valid
# estimator of the Shannon diversity.
# It's even more important to come up with
# the standard error

# Shockingly, until recently, there were no tools
# to estimate Shannon diversity
# in the presence of an ecological/microbial network!

# DivNet is a new tool that allows you to estimate this

# Check out our preprint: Willis & Martin (2018+), bioRxiv

# Let's load the package DivNet
devtools::install_github("adw96/DivNet")
library(DivNet)

# DivNet is very flexible, but by default
# DivNet will estimate the
# microbial network and use it only to
# adjust the standard errors on diversity estimates
dv_water <- divnet(water)
dv_water
plot(dv_water$shannon)

# Let's compare this to the naive approach of
# just "plugging in" the observed proportions
# to the Shannon diversity  formula

plot(water %>% breakaway::shannon)

# You will notice that the estimates are the same
# as previously, but the error bars differ significantly
# i.e., *there are error bars*

# The error bars were the same as previously
# because we haven't told DivNet anything about the experimental
# design. Here we observed
# 4 different water systems, so we will add this
# as a covariate
dv_water <- water %>%
  divnet(x = "SampleType")

# 


# To check that we are in the right directory, we load in the data:
otu_table <- read.table("mask_data.txt",header=T)
covariates <- read.table("mask_covariates.txt",header=T)
## If an error was thrown, that means that the data was not in the
## working directory. Ensure that the data is in filepath specified by
## getwd(), or setwd([where you have saved your data])

## This is data from 63 different lakes

## Have a look at the first few rows/columns your data and your covariates
otu_table[1:4,1:4]
covariates[1:4,]
## data is very big! To see how big, go
dim(otu_table)

# let's look at 1 of the variables in more detail
head(covariates$Year)
## Huh? Levels? What type of object is Year, anyway?
class(covariates$Year)
## Year is what's called a "Factor," i.e. a categorical variable. But
## what if we want to consider it as a continuous variable? We need to
## force it into a numeric. We do this as follows
Year_num <- as.numeric(substr(as.character(covariates$Year),4,5))
Year_num
## Feel like you're falling behind? Don't stress! If you want to learn the
## details you can go back later and figure out exactly what the functions
## as.character, substr, and as.numeric all do, and why I chose to combine
## them in this way.

## Take a moment to congratulate yourself on loading a package,
## opening data, and exploring the dataset! That's huge! Woohoo!


##### CREATE FREQUENCY TABLES

# How could we begin to summarise these 63 lakes?There are
# 1549 different taxa to look at!

# Let's start by looking at the species richness of the samples

## Load the package we need, breakaway. 
# Because we want the version of the package on github, 
# we need to use another package to be able to get it from github
install.packages("devtools")
require(devtools)
install_github("adw96/breakaway", force = T) # download it
require(breakaway) # load it


## As mentioned in the lecture, the format of the data
# needed for species richness estimation is frequency count tables
# Let's get them now
frequency_count_list <- build_frequency_count_tables(otu_table)

## Have a look at the frequency tables ("frequency count data"). I've put them
## in a list so we need to use double square brackets [[63]] to refer to the 63rd one
frequency_count_list[[63]]
## Interpretation: In this sample, there were 57 different species observed
## only once (singletons), 25 different species observed only twice, ...,
## 1 species observed 171 times

##### ESTIMATE SPECIES RICHNESS
## Let's run breakaway on the first frequency count table
breakaway(frequency_count_list[[1]])
## You should get some output to screen, including your estimate & s.e., and a
## plot of the fits to the ratios. Note that it is not a fit to the frequencies,
## it is a fit to the ratios of frequencies. You would never need to include
## this type of plot in one of your papers. It is solely for you to check
## for model misspecification. What's model misspecification? If the black
## diamonds don't remotely follow the pattern of the white circles, that's
## model misspecification.

## The reference for breakaway is
## Willis, A. and Bunge, J. (2015). Estimating Diversity via Frequency Ratios. Biometrics 71 1042-1049.
## Feel free to email me with questions!

## Sometimes, breakaway's usual procedure doesn't work, that is, it gives
## a negative estimate, which is of course silly. In that case, breakaway
## returns a different model's result. It's called the WLRM. There isn't a
## picture. Here is an example of a case where breakaway returns the WLRM.
breakaway(frequency_count_list[[60]])
## breakaway can defer to the WLRM for several reasons. Perhaps there are
## too many singletons. Perhaps there isn't a long enough tail. Perhaps
## there is false diversity. In this case, there was probably not enough data.

## Let's see if this failure was sensitive to the singleton count by running
## breakaway_nof1. This requires no singleton count (implicit is that the
## singleton count was erroneous) and predicts it from the other frequencies.
## Here is an example.
breakaway_nof1(frequency_count_list[[60]][-1,])

## The reference for this method:
## Willis, A. (2016). Species richness estimation with high diversity but spurious singletons.

## breakaway_nof1 is an exploratory tool for assessing sensitivity of
## breakaway to the singleton count. You should not use it for diversity
## estimation -- only diversity exploration.




## Let's move on to looking at the Objective Bayes procedures of
## Kathryn Barger. It's Oprah's favourite species richness estimator!

## There are 4 different types of objective bayes estimates, due to
## 4 different models. If you have time play with all of them!
## For now we are just going to look at the negative binomial.

#objective_bayes_poisson(frequency_count_list[[60]])$results
#objective_bayes_geometric(frequency_count_list[[60]])$results
#objective_bayes_mixedgeo(frequency_count_list[[60]])$results

#install.packages("MASS"); require(MASS) ## you need MASS for this to work
# if R had to restart, go `require(breakaway)` again
bayesian_results <- objective_bayes_negbin(frequency_count_list[[1]], answers = T)
## (Don't worry about those warnings. We're working on them -- sorry!)

## That's a lot of information! Bayesians are very good at generating
## a lot of information. This is because rather than a single estimate
## you see the distribution of estimates (remember that the Bayesian
## paradigm believes the parameter to be random => it has a distribution)


## Let's talk about some of the information:
bayesian_results$results
# $results mode.N/mean.N/median.N : the mode/mean / median estimate
# $results L/UCI.N:  A 95% percent interval estimate for the richness
# The picture at the bottom: The distribution of estimates

## The reference for this method is
## Barger & Bunge. (2010). Objective Bayesian estimation for the number
##    of species. Bayesian Analysis. 5(4), 765-785.

## EVENNESS

## The above discussion focused exclusively on richness. Let's look at
## the variability of evenness estimates

## Let's calculate the plug in estimate of Shannon diversity
shannon(frequency_count_list[[1]])

## 4.17, huh? Let's look at how variable it is
set.seed(2) # the following functions are random, so let's set the seed (allows reproducibility)
resample_estimate(otu_table[,1], shannon)
resample_estimate(otu_table[,1], shannon)
resample_estimate(otu_table[,1], shannon)

## Hmmm, doesn't look too variable! Let's look at a lot of them
par(mfrow=c(1,1))
hist(replicate(200, resample_estimate(otu_table[,1], shannon)))
## (replicate says: "do this 200 times")
## Yikes, that's some negative skew! That suggests that we have
## risks in randomly observing really low Shannon diversity estimates
## This is an important thing to keep in mind when analysing data
## eg. Did we just have bad luck with the sample from the patient
##     taking the drug? Or is the drug causing the effect?

## BTW, you can use the above to look sample size variability as well
## If you have very different numbers of reads across samples, this can$
## introduce additional variability. Let's look at the distribution
## of reads in the current dataset.
ns <- unlist(lapply(frequency_count_list, function(x) sum(x[,2])))
hist(ns)
## Well, but a lot of variability! Lets account for it in looking at
## our distribution of Shannon diversity estimates
set.seed(8)
hist(replicate(500, resample_estimate(otu_table[,1], shannon, my_sample_size = ns)))
## That's a lot more variability than we saw originally!

## Going forward with your analyses, don't forget to account for
## variability in your estimates. Bootstrap standard errors are better
## than nothing!
sd(replicate(500, resample_estimate(otu_table[,1], shannon, my_sample_size = ns)))

## An ongoing mission of the breakaway package is to develop better estimates for
## diversity indices than the current "plug-in" estimates
# The function alpha_better() is currently under development --
# it's a way to upscale alpha diversity for unobserved species
# Have a play with it!

alpha_better(frequency_count_list[[1]], 0) # Species richness
alpha_better(frequency_count_list[[1]], 2) # Inverse Simpson
# compare to
hill(frequency_count_list[[1]], 0)
hill(frequency_count_list[[1]], 2)

## Let's see it over a range
x <- seq(from = 0.001, to = 5, length.out = 100)
plot(x, hill(frequency_count_list[[1]], x), type = "l")
points(x, alpha_better(frequency_count_list[[1]], x), col="red", type="l")
## These are Hill numbers -- a function over alpha diversity:
## Species richness is on the LHS (low q values), evenness on the RHS
## You can see how plug-in alpha div (black) is low, because it 
## doesn't account for unobserved species. The red (upscaled) 
## estimates correct for this -- thus why they're higher!

# Come talk to me if you want to know more -- details and docs 
# coming soon... But you should still pester me :) Basically,
# the estimates are good but I'm still figuring out the standard
# errors


# Please consider following the breakaway package development on github
# to get updates about our progress!

## COMPARISONS

## We spent a lot of time looking at each sample's alpha diversity values

## Let's do some comparisons across samples... accounting for variability
## of course!

## We're going to do this procedure for shannon evenness, and estimate
## it using the plug-in estimate. This is not meant to imply that you should 
## be interested in Shannon! Just that some people are.

## Feel free to substitute in whatever alpha diversity/evenness/richness
## procedure you're interested in!

## To do this we will start by iterating on Shannon on every one of our samples.
## It may take a minute or two to run. Feel like a stretch? Go ahead!

## This section may take a while. Do you know if your neighbour likes bagels?
## Would they order a croissant over a bagel? Now is a great time to find out!
estimates_shannon <- matrix(NA,nrow=dim(otu_table)[2],ncol=4)
rownames(estimates_shannon) <- colnames(otu_table)
colnames(estimates_shannon) <- c("shannon_est","shannon_seest","shannon_lcb","shannon_ucb")
for (i in 1:dim(otu_table)[2]) {
  samples <- replicate(500, resample_estimate(otu_table[,i], shannon, my_sample_size = ns))
  estimates_shannon[i,1] <- mean(samples)
  estimates_shannon[i,2] <- sd(samples)
  estimates_shannon[i,3:4] <- quantile(samples, c(0.025, 0.975))
}
## Here gives us our estimates and standard errors so we can have a look
estimates_shannon[,1:2]
## We can even plot intervals. The x-axis is just against the
## enumeration of the lakes. It is not meaningful.
## The following plotting command may be different to ones you have seen
## before. It's a really quick way of visualizing the error in your samples
## and quickly spotting outliers. NOTE: Outliers have SMALL LINES, which means
## high precision, and are generally far away from points near them.
betta_pic(estimates_shannon[,1], estimates_shannon[,2])

## No obvious outliers here.

## Lets look at the effect of summer samples
col_by_seasons <- ifelse(covariates$Season=="Autum","black",ifelse(covariates$Season=="Spring","pink","red"))
betta_pic(estimates_shannon[,1], estimates_shannon[,2], mycol = col_by_seasons)
legend("bottom",c("Spring","Summer","Autumn"),col=c("pink","red","black"),cex=0.8,lwd=3,box.col=F)
## Don't forget that because we are plotting diversity *estimates*, we
## need to plot *lines* (i.e. confidence intervals) not points (point
## estimates). That's very important.


## We're now going to create our "design" matrix, a.k.a. "X matrix."
## To do this we're going to cheat a little. We're going to use an
## existing R method, lm, to save us the hassle.
## We are going to investigate the effect of temperature and site
covar_matrix <- model.matrix(lm(rnorm(dim(covariates)[1])~covariates$Site*(covariates$Season=="Summ")))
head(covar_matrix)
colnames(covar_matrix) <- c("Int", "SiteP", "Summer", "SitePSummer")
## Details (skip if uninterested): lm(y~x+z) is a function that fits a
## regression line to y using the variables x and z. It has a very nice
## interface that saves you doing the hard work (aka "math") to create your
## design matrix. The function we are going to use has no such nice interface.
## For that reason, we steal the design matrix out of lm()'s implementation.
## lm(y~x*z) looks at interactions. In this case, we are seeing if the
## evenness trend as a function of temperature changes depending on the site

## Easter egg: if you can put your own data (from your own research) in the
## *exact* same form as the above example, you don't need to understand what
## happens below. You can just copy it :)
## Hint: The biggest problem you may have in implementing the above is
## the ordering of the data and the covariates. The following:
colnames(otu_table)==rownames(covariates)
## being true in every spot is a necessary but not sufficient condition for
## the following to work properly.

#### MODEL SPECIES RICHNESS, TEST SIGNIFICANCE, INTERPRET RESULTS
## Let's go ahead and try to fit our model
results <- betta(estimates_shannon[,1],estimates_shannon[,2],covar_matrix)
## Let's take a look at the results
results$table
## Interpretation is idential to regression.

## Here are some things to quickly note:
## Firstly, significance of Summer means that the average
## Shannon diversity of the summer samples is significantly
## less than non-summer AT SITE L. It's less by 0.42 on
## average.

## The non-signif of the SiteP variables means there
## is no significant difference between Sites P and L,
## in both summer and not... AFTER ACCOUNTING FOR 
## ESTIMATION ERROR!

## Notice that a regression on the Shannon estimates
## notes more evidence against the "site has no affect" null:
summary(lm(estimates_shannon[,1] ~ covariates$Site*(covariates$Season=="Summ")))$coef[,c(1,4)]
## Not such a big deal in this case, but may be with more
## marginal results/more variables.

## You should always use betta() when modelling and doing
## inference on functions of your OTU table!

## sigsq_u (in full results output) being significant means there is still *heterogeneity*
## in the lakes. "Not all lakes have the same Shannon diversity
## even after accounting for season & site."
##  (microscale heterogeneity, pH, chemical factors, depth...?)


## If you are unfamiliar with regression analysis, we would recommend
## taking an introductory course in regression analysis. You will
## not regret it (after some period of time)!

### Congratulations! You have finished almost all of the tutorial!

## An important note: I'm not pushing Shannon.
## You can repeat all of the above analysis with your own
## favourite index. 


## Exercise: define a function to calculate the inverse
## of the plug-in Simpson index  estimate

## To estimate the bootstrap mean and standard error
## of the Simpson index, we do the following
simpson_resamples <- replicate(100, resample_estimate(otu_table[,1], simpson))
hist(simpson_resamples)
mean(simpson_resamples)
sd(simpson_resamples)

## this is just for the first column
## you would do this for every sample (in a loop)
## then use these as inputs to betta

## The below shows you how to run CatchAll and read the output into R
## CatchAll is very fussy, because it was built for an outdated operating system
##

##### Optional: CATCHALL

# CatchAll is an old program for estimating species richness.
# It's very good, but it is not open source.
# I (Amy) am in the process of writing it into R, because it is currently a real pain
# to use. Keep checking back at github.com/adw96/CatchAll for progress

# First, download the CatchAll executable from
# http://northeastern.edu/catchall/downloads/CatchAll.4.0/CatchAllGUI.exe
# and place it in your current directory

# Next, download
# runcatchall.sh
# from
# https://github.com/adw96/stamps
# and place it in your current directory

## CatchAll is computionally intensive. For this reason, we will run it
## on only 4 of the samples.
## We are going to create a directory to write our frequency tables to,
## write our frequency tables, run CatchAll on them, read the CatchAll
## output back into R, then analyze the same data again using the
## CatchAll richness estimates

## We are going to be clever and run command line commands via
## R. So the below command is equivalent to writing "mkdir catchalltables"
## in the command line.
system("mkdir catchalltables")

## CatchAll is computionally intensive. For this reason, we will run it
## on only 4 of the samples.
## We are going to create a directory to write our frequency tables to,
## write our frequency tables, run CatchAll on them, read the CatchAll
## output back into R, then analyze the same data again using the
## CatchAll richness estimates
for (i in 1:4) {
  tmp <- frequency_count_list[[i]] ## catchall is really fussy! Need numbers not strings/factors
  tmp2 <- cbind(as.numeric(levels(tmp[,1]))[-1],tmp[,2])
  write.table(tmp2,paste(getwd(), "/catchalltables/",
                         names(frequency_count_list)[i],".csv",sep=""),sep=",",row.names = FALSE,col.names = FALSE)
}
system("cp runcatchall.sh catchalltables")
system("cp CatchAllCmdL.exe catchalltables")
## This next step runs CatchAll on every file; it may take a few minutes
system("cd catchalltables; ./runcatchall.sh")

### IF THIS DIDN'T WORK *and* YOU WANT TO LEARN TO USE CATCHALL**
## 1. If you're on a Mac/Unix: You probably don't have mono (an emulator)
##    Go to http://www.mono-project.com/ and install it. 
## 2. If you're on Windows: Hmmm. You're going to need a different version. 
##    Download CatchAllcmdW.exe from http://www.northeastern.edu/catchall/downloads.html
##    Then run CatchAllcmdW.exe [inputfilename] [outputpath]
## 3. If neither of the above help: Sorry! Call me/a TA over. Or, come back to
##    it later and open up butterfly_BestModelsAnalysis.csv (an example output)
##    to see

## ** If CatchAll isn't so important, feel free to check out. You've learnt
## a lot already!

## We just ran CatchAll without interfacing with the command line.
## Woohoo! Now we'll read them in.

## We tell R to search in the directory "catchalltables" for any files ending
## in "*BestModelsAnalysis.csv" We then pull out the relevant columns
## into a matrix called estimates_catchall
txt.sources = list.files(path=paste(getwd(),"/catchalltables/",sep=""),
                         pattern="*BestModelsAnalysis.csv",recursive=TRUE)
myread <- function(x)  {
  y <- read.csv(paste("catchalltables/",x,sep=""),stringsAsFactors=FALSE)
  return(as.numeric(y[1,4:7]))
}
estimates_catchall <- matrix(sapply(txt.sources,myread),byrow=TRUE,ncol=4)
rownames(estimates_catchall) <- paste("Sample", 1:4, sep="")
colnames(estimates_catchall) <- c("catchall_est","catchall_seest","catchall_lcb","catchall_ucb")
estimates_catchall

## We can now run betta in the same way!
betta(estimates_catchall[,1], estimates_catchall[,2])$table
## Note that these estimates are homogeneous (all the same!)
## probably because there are only 4
## Here we haven't included any covariates (it's unlikely that we
## could find anything meaningful with 4 data points), but if you did
## this with your own data you would include covariates in the same
## way that we did for breakaway.

## THANKS FOR STAYING THROUGH TO THE END! I'd love to hear if you
## hated it/loved it/worst day ever -- please give me feedback to 
## better help you!! :) Amy
