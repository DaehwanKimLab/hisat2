#!/usr/bin/env python

import os, sys, subprocess, re
import inspect, operator
from argparse import ArgumentParser, FileType

def checkNTloc(fasta_fileName,var_fileName):
    seq = ""
    for line in open(fasta_fileName,'r'):
        if line[0] == '>':
            continue
        seq += line.strip()

    cyp_var_file = open(var_fileName,'r')
    cyp_var_dict = makeVarDict(cyp_var_file)
    cyp_var_file.close()

    print "len:", len(seq)
    varsPos = set()
    varsNeg = set()

    for varList in cyp_var_dict.values():
        for var in varList:
            if ">" in var: # is SNP
                posNt = int(var[:-3])
                ntChange = var[-3:].replace('>','')
                assert len(ntChange) == 2
                for nt in ntChange:
                    assert nt in "ACGT"

                if posNt > 0:
                    varsPos.add(str(posNt) + '->' + ntChange[0])
                else:
                    assert posNt < 0
                    varsNeg.add(str(posNt) + '->' + ntChange[0])
                    
            elif "del" in var: # is deletion
                posNt = var.split('del')[0].split('_')
                posNt = [int(p) for p in posNt]
                ntDel = var.split('del')[1]
                for nt in ntDel:
                    assert nt in "ACGT"

                if len(posNt) == 1: # single nt deletion
                    assert len(ntDel) == 1
                    if posNt[0] > 0:
                        varsPos.add(str(posNt[0]) + '->' + ntDel)
                    else:
                        assert posNt[0] < 0
                        varsNeg.add(str(posNt[0]) + '->' + ntDel)

                else: # mutliple nt deletion
                    assert len(posNt) == 2
                    assert posNt[1] - posNt[0] + 1 == len(ntDel)
                    ntDelList = list(ntDel)
                    for i in range(posNt[0],posNt[1] + 1):
                        if i > 0:
                            varsPos.add(str(i) + '->' + ntDelList.pop(0))
                        else:
                            assert i < 0
                            varsNeg.add(str(i) + '->' + ntDelList.pop(0))
                    assert len(ntDelList) == 0
                    
            else:
                assert ("ins" in var) or ("None" in var)
                continue
            
    '''varsPos = set(["310-G", "746-C", "843-T", "997-C", "1039-C", "1513-C", "1661-G", "1724-C", "1869-T", "1978-C", "2470-T", "2575-C", "2850-C", "3828-G", "4115-C", "4180-G"])'''
    
    scorePos = {} # { position offset : number of alignments } for positive positions
    for i in range(-len(seq), len(seq)):
        align_score = 0
        for var in varsPos:
            pos, base = var.split('->')
            pos = int(pos)
            
            try:
                seq[pos-i]
            except IndexError:
                continue
            
            if seq[pos-i] == base:
                align_score += 1

        scorePos[i] = align_score

    print "Positive postitions offset: %d" % \
          (len(seq) - max(scorePos.iteritems(), key=operator.itemgetter(1))[0])
    

    scoreNeg = {} # { position offset : number of alignments } for negative positions
    for i in range(-len(seq), len(seq)):
        align_score = 0
        for var in varsNeg:
            pos, base = var.split('->')
            pos = int(pos)
            
            try:
                seq[pos-i]
            except IndexError:
                continue
            
            if seq[pos-i] == base:
                align_score += 1

        scoreNeg[i] = align_score

    print "Negative postitions offset: %d" % \
          (len(seq) - max(scoreNeg.iteritems(), key=operator.itemgetter(1))[0])


'''
    for i in range(-len(seq), len(seq)):
        Pass = True
        for var in varsNeg:
            pos, base = var.split('->')
            pos = int(pos)
            if pos < i or pos - i >= len(seq):
                Pass = False
                break

            if seq[pos-i] != base:
                Pass = False
                break

        if Pass:
            print "Negative positions offset: %d" % -i
'''

def create_map(seq):
    seq_map = {}
    count = 0
    for i in range(len(seq)):
        bp = seq[i]
        if bp == '.':
            continue
        assert bp.upper() in "ACGT"
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

def makeSnp(oldSeq, pos, oldNt, newNt):
    assert oldSeq[pos] == oldNt
    newSeq = oldSeq[:pos] + newNt + oldSeq[pos+1:]
    assert len(newSeq) == len(oldSeq)
    return newSeq

def makeDel(oldSeq, left, right, toDel):
    assert right - left + 1 == len(toDel)
    assert oldSeq[left:right + 1] == toDel
    newSeq = oldSeq[:left] + '.'*len(toDel) + oldSeq[right + 1:]
    assert len(newSeq) == len(oldSeq)
    return newSeq
    
def makeIns(oldSeq,left,right,toIns):
    assert right - left - 1 >= len(toIns)
    for nt in oldSeq[left + 1:right]:
      assert nt == '.'
    remDots = right - left - 1 - len(toIns)
    newSeq = oldSeq[:left + 1] + toIns + '.'*remDots + oldSeq[right:]
    assert len(newSeq) == len(oldSeq)
    return newSeq
    

def main():
    cyp_var_file = open("cyp2d6.var",'r')
    cyp_var_dict = makeVarDict(cyp_var_file)
    cyp_var_file.close()
             
    '''
    for item in _var_dict.items():
        print(item)
    '''
    
    cyp_faFile = open("cyp2d6.fasta",'r')
    cyp_seq = extractSeq(cyp_faFile)
    cyp_faFile.close()
    preBackbone_seq = ''
    

    msfTable = {}

    # Building backbone structure (augment length with insertions)
    longestIns = {} # { key = position : value = length }
    for allele,varList in cyp_var_dict.items():
        for var in varList:
            if not "ins" in var:
                continue
            pos = var.split('ins')[0].split('_')
            pos = [int(p) for p in pos]
            ntIns = var.split('ins')[1]
            correctFormat = len(pos) == 2 and pos[1] - pos[0] == 1
            if not correctFormat:
                correctFormat = len(pos) == 1
            try:
                assert correctFormat
            except:
                print >> sys.stdout, "Incorrect format for insertion: variation %s on allele %s" % (var, allele)
                continue

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
        '''print(posInsList[i])'''
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
    print('\n\n')
    '''print >> open("preBackbone_seq.fasta",'w'), preBackbone_seq'''
    map_cyp = create_map(preBackbone_seq) # { Index of bp in original seq : Actual index in string }
    
    
############################################################################################################################################################

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
                isRef = True

            if isSnp:
                pos = int(var[:-3])
                dbPos = pos
                ntChange = var[-3:].replace('>','')
                assert len(ntChange) == 2
                for nt in ntChange:
                    assert nt in "ACGT"

                if pos > 0:
                    pos = pos + 1618
                else:
                    pos = pos + 1619

                if pos < 0:
                    print >> sys.stdout, "Warning: position %d out of bounds" % (dbPos)
                    print >> sys.stdout, "\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                    continue
                    
                try:
                    assert(preBackbone_seq[map_cyp[pos]] == ntChange[0]) # nt at pos in seq must match database
                except:
                    print >> sys.stdout, "Warning: position %d in sequence contains %s, but expected %s from database" % (dbPos, preBackbone_seq[map_cyp[pos]], ntChange[0])
                    print >> sys.stdout, "\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                    continue
                
                # Adding to msf table
                if not allele in msfTable:
                    msfTable[allele] = makeSnp(preBackbone_seq, map_cyp[pos], ntChange[0], ntChange[1])
                else:
                    msfTable[allele] = makeSnp(msfTable[allele], map_cyp[pos], ntChange[0], ntChange[1])
                    
            elif isDel:
                pos = var.split('del')[0].split('_')
                pos = [int(p) for p in pos]
                if len(pos) == 1: # Handle single deletion with format for multi deletion with one location (e.g. [1707] -> [1707,1707])  
                    pos.append(pos[0])
                assert len(pos) == 2
                dbPos = pos
                ntDel = var.split('del')[1]
                for nt in ntDel:
                    assert nt in "ACGT"

                for i in range(len(pos)):
                    if pos[i] > 0:
                        pos[i] = pos[i] + 1618
                    else:
                        pos[i] = pos[i] + 1619

                for i in range(len(pos)):
                    if pos[i] < 0:
                        print >> sys.stdout, "Warning: position %d out of bounds" % (dbPos[i])
                        print >> sys.stdout, "\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                        continue
                        
            
                try:
                    assert pos[1] - pos[0] + 1 == len(ntDel)
                except:
                    print >> sys.stdout, "Incorrect deletion data with %s on allele %s. Skipping variation." % (var, allele)
                    continue
                            
                try:
                    assert preBackbone_seq[ map_cyp[pos[0]] : map_cyp[pos[1]] + 1 ] == ntDel
                except:
                    print >> sys.stdout, "Warning, positions %d to %d in sequence contains %s, but expected %s from database" % \
                          (dbPos[0], dbPos[1], preBackbone_seq[ map_cyp[pos[0]] : map_cyp[pos[1]] + 1 ], ntDel)
                    print >> sys.stdout, "\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                    continue


                # Adding to msf table
                if not allele in msfTable:
                    msfTable[allele] = makeDel(preBackbone_seq, map_cyp[pos[0]], map_cyp[pos[1]], ntDel)
                else:
                    msfTable[allele] = makeDel(msfTable[allele], map_cyp[pos[0]], map_cyp[pos[1]], ntDel)

                        
            elif isIns:
                pos = var.split('ins')[0].split('_')
                pos = [int(p) for p in pos]
                if len(pos) == 1:
                    pos.append(pos[0] + 1)
                assert len(pos) == 2
                dbPos = pos
                assert pos[1] - pos[0] == 1
                ntIns = var.split('ins')[1]
                for nt in ntIns:
                    assert nt in "ACGT"

                for i in range(len(pos)):
                    if pos[i] > 0:
                        pos[i] = pos[i] + 1618
                    else:
                        pos[i] = pos[i] + 1619

                for i in range(len(pos)):
                    if pos[i] < 0:
                        print >> sys.stdout, "Warning: position %d out of bounds" % (dbPos[i])
                        print >> sys.stdout, "\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                        continue


                # Adding to msf table
                if not allele in msfTable:
                    msfTable[allele] = makeIns(preBackbone_seq, map_cyp[pos[0]], map_cyp[pos[1]], ntIns)
                else:
                    msfTable[allele] = makeIns(msfTable[allele], map_cyp[pos[0]], map_cyp[pos[1]], ntIns)


            else:
                assert isRef
                assert not allele in msfTable
                msfTable[allele] = preBackbone_seq


    msfFile = open('cyp2d6.msf','w')
    for allele,msf_seq in msfTable.items():
        print >> msfFile, "%s\t%s" % (allele, msf_seq)

    msfFile.close()

checkNTloc("cyp2d6.fasta","cyp2d6.var")
