#!/usr/bin/env python

import os, sys, subprocess, re
import inspect
from argparse import ArgumentParser, FileType

def create_map(seq):
    seq_map = {}
    count = 0
    for i in range(len(seq)):
        bp = seq[i]
        if bp == '.':
            continue
        assert bp in "ACGT"
        seq_map[count] = i
        count += 1
    return seq_map

def splitString(someStr,posList):
    posList.insert(0,-1)
    posList.append(len(someStr) - 1)
    splitStr = []
    for i in range(len(posList) - 1):
        left = posList[i] + 1
        right = posList[i+1] + 1
        splitStr.append(someStr[left:right])

    return splitStr

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
            alleleVarDict[alleleName] = set(varList)
        except:
            print >> sys.stdout, ("Warning, %s allele is already represented" % alleleName)
            alleleVarDict[alleleName] = alleleVarDict[alleleName] | set(varList)

    return alleleVarDict

def main():
    cyp_var_file = open("cyp2d6.web.output",'r')
    cyp_var_dict = makeVarDict(cyp_var_file)
             
    '''
    for item in _var_dict.items():
        print(item)
    '''
    
    cyp_faFile = open("cyp2d6.fasta",'r')
    cyp_seq = extractSeq(cyp_faFile)
    preBackbone_seq = ''
    

    msfTable = []
    alleleIDdict = {} # allele ID corresponds to index in msfTable

    # Building backbone structure (augment length with insertions)
    longestIns = {} # { key = position : value = length }
    for allele,varList in cyp_var_dict.items():
        for var in varList:
            if not "ins" in var:
                continue
            pos = var.split('ins')[0].split('_')
            pos = [int(p) for p in pos]
            ntIns = var.split('ins')[1]
            try:
                assert (len(pos) == 1) or (len(pos == 2) and pos[1] - pos[0] == 1)
            except:
                print >> sys.stdout, "Incorrect format for insertion: variation %s on allele %s" % (var, allele)

            # convert to position in string
            if pos[0] > 0:
                pos = pos[0] + 1618
            else:
                pos = pos[0] + 1619
                
            # Make dictionary of longest insertions
            if not pos in longestIns:
                longestIns[pos] = len(ntIns)
            else:
                if len(ntIns) > longestIns[pos]:
                    longestIns[pos] = len(ntIns)

    print(longestIns)
    posInsList = sorted(longestIns.keys())
    print(posInsList)
    
    splitSeq = splitString(cyp_seq,posInsList)
    posInsList = posInsList[1:-1]

    for i in range(len(posInsList)):
        # print(posInsList[i])
        splitSeq[i] += '.' * longestIns[posInsList[i]]

    for subseq in splitSeq:
        try:
            assert len(subseq) > 0 and not subseq.startswith('.')
            preBackbone_seq += subseq
        except:
            continue
    # pre-backbone built

    print(len(cyp_seq))
    print(len(preBackbone_seq))

############################################################################################################################################################            
'''
    for allele,varList in cyp_var_dict.items():
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
                            assert(cyp_seq[pos + 1618] == ntChange[0]) # nt at pos in seq must match database
                    else:
                            assert(cyp_seq[pos + 1619] == ntChange[0])
                except:
                    print >> sys.stdout, "Warning: position %d in sequence contains %s, but expected %s from database" % (pos, cyp_seq[pos + 1618] if pos > 0 else cyp_seq[pos + 1619], ntChange[0])
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
                                assert(cyp_seq[pos[0] + 1618 : pos[1] + 1618 + 1] == ntDel)
                        else:
                                assert(cyp_seq[pos[0] + 1619 : pos[1] + 1619 + 1] == ntDel)
                    except:
                        print >> sys.stdout, "Warning, positions %d to %d in sequence contains %s, but expected %s from database" % \
                              (pos[0], pos[1], cyp_seq[pos[0] + 1618 : pos[1] + 1618 + 1] if pos[0] > 0 else cyp_seq[pos[0] + 1619 : pos[1] + 1619 + 1], ntDel)
                        print >> sys.stdout, "\tError occured on variation %s on allele %s" % (var, allele)

                else: # Single Deletion
                    assert len(ntDel) == 1
                    try:
                        if pos[0] > 0:
                                assert(cyp_seq[pos[0] + 1618] == ntDel)
                        else:
                                assert(cyp_seq[pos[0] + 1619] == ntDel)
                    except:
                        print >> sys.stdout, "Warning: position %d in sequence contains %s, but expected %s from database" % \
                              (pos[0], cyp_seq[pos[0] + 1618] if pos[0] > 0 else cyp_seq[pos[0] + 1619], ntDel)
                        print >> sys.stdout, "\tError occured on variation %s on allele %s" % (var, allele)
                        
            else:
                continue

'''

main()
