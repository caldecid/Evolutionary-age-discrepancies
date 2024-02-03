## METADATA ===============================================================
## Description: This script sources the working directory and the files 
## 
## R version: 4.0.2 for Windows
## Date: 2022-12-05
## Author: Carlos Calderon
##=======================================================================##

##libraries
library(ape)
library(phytools)
library(readr)
library(picante)
library(tidyverse)
library(phangorn)
library(TreeTools)
library(tidytree)
library(readr)
library(FossilSim)
library(TreeSim)
library(ggridges)
library(ggpubr)
library(multcompView)
library(cowplot)
library(scales)
library(future)
library(future.apply)
library(TreeSimGM)
library(HDInterval)
library(arules)
library(gridExtra)
library(grid)
library(spatstat)
library(modi)

# set working directory -----------
wd <- "/home/alunos/Documents/caldecid/Evolutionary-age-discrepancies"


### [] create variables for directory relative location -------

####### results

##scripts

scripts <- "results/scripts"

functions <- "Functions"


functions_dir <- list.files(paste0(scripts, "/", functions), full.names = TRUE)


for(i in seq_along(functions_dir)){
  source(functions_dir[i])
}


##data 

# results/data/raw
raw <- "results/data/raw"

#results/data/processed
pro <- "results/data/processed"



# results/data/meta
meta <- "results/data/metadata"


###### figures files

figures <- "/text/figures"


##### common folders

# quantifying error files
q_error <- "/quantifying_error"


##stochastic function files
st_sp <- "/stochastic_sp_age"

##conservation status files
c_status <- "/conservation_status"








