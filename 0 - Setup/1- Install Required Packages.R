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


check_ggplot2<-require(ggplot2)
if (check_ggplot2==FALSE) {
  install.packages("check_ggplot2")
}

check_splines<-require(splines)
if (check_splines==FALSE) {
  install.packages("splines")
}

check_nlme<-require(nlme)
if (check_nlme==FALSE) {
  install.packages("nlme")
}

check_Matrix<-require(Matrix)
if (check_Matrix==FALSE) {
  install.packages("Matrix")
}

# clear environment
rm(list = ls())









