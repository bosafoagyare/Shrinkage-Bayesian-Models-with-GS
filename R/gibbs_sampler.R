#(1) %%%%%%%%%%%%%% LOAD LIBRARIES %%%%%%%%%%%%%%
#--------- Load algorithms script
source("sslab.R")      # contains spike--and-slab algos
source("blasso.R")     # contains Bayesian Lasso algos
source("utils.R")      # contains utility functions

#--------- Load, otherwise install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lars, tidyverse, foreach, coda, glmnet)


#--------- Check for OS and install/load appropriate parallel package
ifelse(Sys.info()["sysname"] == "Windows",
               pacman::p_load(doParallel), # for Windows
               pacman::p_load(doMC))       # for Mac/Linux
  

#--------- Check for OS and setup appropriate parallel system
cores = detectCores()
if(Sys.info()["sysname"] == "Windows"){
  cl <- makeCluster(cores-2)
  registerDoParallel(cl) # for Windows
}else{
  registerdoMC(cores-2)
}
       

#--------- Check for OS and close appropriate parallel system 
if(Sys.info()["sysname"] == "Windows") stopCluster(cl)














