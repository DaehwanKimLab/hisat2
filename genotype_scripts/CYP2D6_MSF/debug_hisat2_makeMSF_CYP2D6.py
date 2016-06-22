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

    msfTable = []

##    for item in cyp2d6_var_dict.items():
##        print(item)

    cyp2d6_faFile = open("cyp2d6.fasta",'r')
    cyp2d6_seq = extractSeq(cyp2d6_faFile)

    for allele,varList in cyp2d6_var_dict.items():
        for var in varList:
            isSnp = False
            isDel = False
            isIns = False
        
            if ">" in var:
                isSnp = True
            elif "del" in var:
                isDel = True
            elif "ins" in var:
                isIns = True
            else:
                assert("None" in var)

            if isSnp:
                pos = int(var[:-3])
                ntChange = var[-3:].replace('>','')
                assert len(ntChange) == 2
                for nt in ntChange:
                    assert nt in "ACGT"

                # @Daehwan - checking to make sure nt changes are compatible with seq (shows database errors in nt positions)
                try:
                    if pos > 0:
                            assert(cyp2d6_seq[pos + 1618] == ntChange[0]) # nt at pos in seq must match database
                    else:
                            assert(cyp2d6_seq[pos + 1619] == ntChange[0])
                except:
                    print >> sys.stdout, "Warning: position %d in sequence contains %s, but expected %s from database" % (pos, cyp2d6_seq[pos + 1618] if pos > 0 else cyp2d6_seq[pos + 1619], ntChange[0])
                    print >> sys.stdout, "\tError occured on variation %s on allele %s" % (var, allele)
                    
            elif isDel:
                pos = var.split('del')[0].split('_')
                pos = [int(p) for p in pos]
                ntDel = var.split('del')[1]
                for nt in ntDel:
                    assert nt in "ACGT"

                if len(pos) > 1: # Multiple Deletions
                    assert len(pos) == 2
                    try:
                        assert pos[1] - pos[0] + 1 == len(ntDel)
                    except:
                        print >> sys.stdout, "Incorrect deletion data with %s on allele %s" % (var, allele)

                    try:
                        if pos[0] > 0:
                                assert(cyp2d6_seq[pos[0] + 1618 : pos[1] + 1618 + 1] == ntDel)
                        else:
                                assert(cyp2d6_seq[pos[0] + 1619 : pos[1] + 1619 + 1] == ntDel)
                    except:
                        print >> sys.stdout, "Warning, positions %d to %d in sequence contains %s, but expected %s from database" % \
                              (pos[0], pos[1], cyp2d6_seq[pos[0] + 1618 : pos[1] + 1618 + 1] if pos[0] > 0 else cyp2d6_seq[pos[0] + 1619 : pos[1] + 1619 + 1], ntDel)
                        print >> sys.stdout, "\tError occured on variation %s on allele %s" % (var, allele)

                else: # Single Deletion
                    assert len(ntDel) == 1
                    try:
                        if pos[0] > 0:
                                assert(cyp2d6_seq[pos[0] + 1618] == ntDel)
                        else:
                                assert(cyp2d6_seq[pos[0] + 1619] == ntDel)
                    except:
                        print >> sys.stdout, "Warning: position %d in sequence contains %s, but expected %s from database" % \
                              (pos[0], cyp2d6_seq[pos[0] + 1618] if pos[0] > 0 else cyp2d6_seq[pos[0] + 1619], ntDel)
                        print >> sys.stdout, "\tError occured on variation %s on allele %s" % (var, allele)
                        
            else:
                continue                 
                # Ignore insertions for now

main()
