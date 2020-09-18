# HISAT2 Read Simulator Expanded

## Author Information

The original simulation code was completed in 2015 by Dr. Daehwan Kim. 
Enhancements by Micah Thornton (Kim Lab) beginning 2020. 

## Overview 

The original HISAT2 Simulator for synthesizing reads from a given genetic sequence file will generate reads 
that come from one of several sequencing kinds.  The primary goal of this branch is to update the 
original simulator which does not currently allow for the simulation of BiSulfite reads for analysis 
for use with the original BSmooth statistical pipeline (Kasper Hansen et al 2012-2020) as well as the simulation
of reads of other kinds. 

## Original Simulator Information 


## Updates and Edits 

09-17-2020:  Added a class for methylation status of CpGs to be generated using a Bernoulli trials approach. 

09-18-2020:  In order to verify the proper functionality of the simulation script, one may use the included Makefile by running 'make Tests'; 
the resulting test information will be outputted to a log file called sim-test-res.log. This is a work in progress, for the intention of testing
various simulated sequencing reads. 

