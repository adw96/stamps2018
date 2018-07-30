# STAMPS 2018 R Tutorials 

Welcome! This subdirectory contains R Tutorials for STAMPS 2018 @ MBL. 

1. Make sure you have downloaded R and RStudio.

2. Open RStudio.

3. Set your filepath to a location on your computer that you in which you are comfortable saving files using:
``` r
setwd("YOUR FILEPATH HERE")
```

4. Run the following code in the R console: 
``` r
utils::download.file("https://github.com/mblstamps/stamps2018/blob/master/R_tutorial/R_tutorial.zip?raw=true", "R_tutorial.zip")
utils::unzip("R_tutorial.zip")
```

5. Move into the working directory that was created:
``` r
setwd("STAMPS2018")
```

6. Before moving on, check to make sure you have 7 ".R" files and 2 ".txt" files:
``` r
list.files()
```

7. If you have that, then load in the data we'll be using for this tutorial:
``` r
covariates <- read.csv("FWS_covariates.txt", sep = "\t")
abundances <- read.csv("FWS_OTUs.txt", sep = "\t", row.names = 1, header = T)
```

8. Next, install all the packages we will be needing for this tutorial:
``` r 
install.packages("dplyr")
install.packages("magrittr")
install.packages("devtools")
install.packages("ggplot2")
install.packages("phyloseq")
install.packages("parallel")
install.packages("foreach")
install.packages("doParallel")
devtools::install_github("adw96/breakaway")
```

9. You are now set up! The tutorials are set up in the following order:
``` 
1_loops.R
2_apply.R
3_functions.R
4_lists.R
5_sourcesystem.R
6_parallel.R
7_markdown.R 
```

Happy coding!
