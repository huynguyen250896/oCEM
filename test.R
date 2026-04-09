source("R/optimizeCOM.R")
source("R/overlapCEM.R")
source("R/ICA-Zscore.R") # default value
library("tidyverse")



exp <- read.table("data_n_code/ecoli/ecoli_data.tsv", sep = "\t", header = FALSE)

num_pc <- optimizeCOM(data = exp, cores = 5, P = 50, verbose=F)

cem <-overlapCEM(data = exp, ncomp = num_pc)