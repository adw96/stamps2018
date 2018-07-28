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
# This was a good way to do it, but are we going to try to do the
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
