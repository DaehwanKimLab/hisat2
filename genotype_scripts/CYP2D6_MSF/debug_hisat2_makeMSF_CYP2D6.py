#!/usr/bin/env python

import os, sys, subprocess, re
import inspect
from argparse import ArgumentParser, FileType

def extractSeq(faFile):
    seq = ""
    for line in faFile:
        if line.startswith(">"):
            continue

        seq += line.strip()

    return seq

def makeVarDict(fname):
    alleleVarDict = {}

    allLines = [line.strip() for line in fname]
    assert allLines[1].upper().startswith("CYP")
    alleleVarDict[allLines[1]] = ["None"] # first allele is reference allele
    
    for line in allLines[2:]:
        assert line.upper().startswith("CYP")
        alleleName = line.split("\t")[0]
        
        try:
            varList = line.split("\t")[1].split(',')
        except IndexError:
            continue
        
        try:
            assert not alleleName in alleleVarDict
            alleleVarDict[alleleName] = varList
        except:
            print >> sys.stdout, ("Warning, %s allele is already represented" % alleleName)
            alleleVarDict[alleleName] = alleleVarDict[alleleName] + varList

    return alleleVarDict

def main():
    cyp2d6_var_file = open("cyp2d6.web.output",'r')
    cyp2d6_var_dict = makeVarDict(cyp2d6_var_file)

##    for item in cyp2d6_var_dict.items():
##        print(item)

    cyp2d6_faFile = open("cyp2d6.fasta",'r')
    cyp2d6_seq = extractSeq(cyp2d6_faFile)

    for allele,varList in cyp2d6_var_dict.items():
        for var in varList:
            if not ">" in var:
                continue # @Daehwan - just look at snp for debugging (checking fasta sequence)

            pos = int(var[:-3])
            ntChange = var[-3:].replace('>','')

            assert re.compile('[ACGT][ACGT]', re.IGNORECASE).search(ntChange)

            # @Daehwan - checking to make sure nt changes are compatible with seq (shows database errors in nt positions)
            try:
		if pos > 0:
                	assert(cyp2d6_seq[pos + 1618] == ntChange[0]) # nt at pos in seq must match database
		else:
			assert(cyp2d6_seq[pos + 1619] == ntChange[0])
            except:
                print >> sys.stdout, "Warning: position %d in sequence contains %s, but expected %s from database" % (pos, cyp2d6_seq[pos + 1618], ntChange[0])

main()
