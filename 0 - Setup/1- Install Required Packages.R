 # Check for Packages

#clear environment
rm(list = ls())

# this.path
check_this_path<-require(this.path)
if (check_this_path==FALSE) {
  install.packages("this.path")
}

# tidyverse
check_tidyverse<-require(tidyverse)
if (check_tidyverse==FALSE) {
  install.packages("tidyverse")
}

# mvtnorm
check_mvtnorm<-require(mvtnorm)
if (check_mvtnorm==FALSE) {
  install.packages("mvtnorm")
}

# readxl 
check_readxl<-require(readxl)
if (check_readxl==FALSE) {
  install.packages("readxl")
}

# mice
check_mice<-require(mice)
if (check_mice==FALSE) {
  install.packages("mice")
}

check_doParallel<-require(doParallel)
if (check_doParallel==FALSE) {
  install.packages("doParallel")
}
 
# clear environment
rm(list = ls())









