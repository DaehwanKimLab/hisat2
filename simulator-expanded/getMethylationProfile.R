#!/bin/bash/Rscript

read.csv("temp_cpg_locs.txt",sep=",",header=FALSE) -> cpgLocations; 

methProbDF <- cbind(cpgLocations, runif(dim(cpgLocations)[1],0,1)); 
write.csv(methProbDF,file="temp_cpg_locs_probs_from_R.csv");
