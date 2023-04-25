#(1) %%%%%%%%%%%%%% LOAD LIBRARIES %%%%%%%%%%%%%%
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lars, tidyverse, foreach, coda)


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





#(2) %%%%%%%% UTILS %%%%%%%%%%%%%%%%%
#------- Data generators
beta.true <- function(){
  temp <- rep(0,(p*5))
  p.nz <- ceiling(5*p*s)
  if(p.nz>0){
    temp[1:p.nz] <- rt(p.nz,2)
  }
  return(temp)
}







