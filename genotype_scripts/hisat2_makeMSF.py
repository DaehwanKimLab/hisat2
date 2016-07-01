#!/usr/bin/env python

import os, sys, subprocess, re
import inspect, operator
from argparse import ArgumentParser, FileType

def checkNTloc(fasta_fileName,var_fileName,gene_name):
    print "\nGene: %s" % gene_name
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
                    try:
                        assert posNt[1] - posNt[0] + 1 == len(ntDel)
                    except AssertionError:
                        print "Incorrect deletion format: %s" % (var)
                        sys.exit(1)
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
                seq[pos+i]
            except IndexError:
                continue
            
            if seq[pos+i] == base:
                align_score += 1

        scorePos[i] = align_score
    oSetPos = max(scorePos.iteritems(), key=operator.itemgetter(1))[0]
    print "Positive postitions offset: %d" % oSetPos
    print "Score: %d out of %d\n" % (scorePos[oSetPos], len(varsPos))
    

    print "Checking negative position offset: %d" % (oSetPos + 1)
    align_score = 0
    oSetNeg = oSetPos + 1
    for var in varsNeg:
        pos, base = var.split('->')
        pos = int(pos)
        
        try:
            seq[pos + oSetNeg]
        except IndexError:
            continue
        
        if seq[pos + oSetNeg] == base:
            align_score += 1
    print "Score: %d out of %d\n\n" % (align_score, len(varsNeg))

    '''
    scoreNeg = {} # { position offset : number of alignments } for negative positions
    for i in range(-len(seq), len(seq)):
        align_score = 0
        for var in varsNeg:
            pos, base = var.split('->')
            pos = int(pos)
            
            try:
                seq[pos+i]
            except IndexError:
                continue
            
            if seq[pos+i] == base:
                align_score += 1

        scoreNeg[i] = align_score
    oSetNeg = max(scoreNeg.iteritems(), key=operator.itemgetter(1))[0]
    print "Negative postitions offset: %d" % oSetNeg
    print "Score: %d out of %d\n\n" % (scoreNeg[oSetNeg], len(varsNeg))
    '''

    if len(varsNeg) == 0 and len(varsPos) != 0:
        return oSetPos, oSetNeg, float(scorePos[oSetPos])/float(len(varsPos)), 1.0
    elif len(varsNeg) != 0 and len(varsPos) == 0:
        return oSetPos, oSetNeg, 1.0, float(align_score)/float(len(varsNeg))
    elif len(varsNeg) == 0 and len(varsPos) == 0:
        return oSetPos, oSetNeg, 1.0, 1.0
    else:
        assert len(varsNeg) != 0 and len(varsPos) != 0
        return oSetPos, oSetNeg, float(scorePos[oSetPos])/float(len(varsPos)), float(align_score)/float(len(varsNeg))
        

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
    '''assert allLines[1].upper().startswith("CYP")
    alleleVarDict[allLines[1]] = ["None"] # first allele is reference allele'''
    
    for line in allLines[1:]:
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
    

def makeMSF(gene_name, oSetPos, oSetNeg):
    cyp_var_file = open("cyp_var_files/%s.var" % gene_name,'r')
    cyp_var_dict = makeVarDict(cyp_var_file)
    cyp_var_file.close()
             
    '''
    for item in _var_dict.items():
        print(item)
    '''
    
    cyp_faFile = open("cyp_fasta/%s.fasta" % gene_name,'r')
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
                print >> sys.stdout, "\tIncorrect format for insertion: variation %s on allele %s" % (var, allele)
                continue

            # convert to position in string
            if pos[0] > 0:
                pos = pos[0] + oSetPos
            else:
                pos = pos[0] + oSetNeg
                
            # Make dictionary of longest insertions
            if not pos in longestIns:
                longestIns[pos] = len(ntIns)
            else:
                if len(ntIns) > longestIns[pos]:
                    longestIns[pos] = len(ntIns)
    
    '''print(longestIns)'''
    posInsList = sorted(longestIns.keys())
    '''print(posInsList)'''
    
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

    '''print(len(cyp_seq))
    print(len(preBackbone_seq))
    print('\n\n')'''
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
                    pos = pos + oSetPos
                else:
                    pos = pos + oSetNeg

                if pos < 0 or pos > len(cyp_seq) - 1:
                    print >> sys.stdout, "\tWarning: position %d out of bounds" % (dbPos)
                    print >> sys.stdout, "\t\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                    continue
                    
                try:
                    assert(preBackbone_seq[map_cyp[pos]] == ntChange[0]) # nt at pos in seq must match database
                except:
                    print >> sys.stdout, "\tWarning: position %d in sequence contains %s, but expected %s from database" % (dbPos, preBackbone_seq[map_cyp[pos]], ntChange[0])
                    print >> sys.stdout, "\t\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
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
                        pos[i] = pos[i] + oSetPos
                    else:
                        pos[i] = pos[i] + oSetNeg

                skipDel = False
                for i in range(len(pos)):
                    if pos[i] < 0 or pos[i] > len(cyp_seq) - 1:
                        print >> sys.stdout, "\tWarning: position %d out of bounds" % (dbPos[i])
                        print >> sys.stdout, "\t\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                        skipDel = True

                if skipDel:
                    continue
                        
            
                try:
                    assert pos[1] - pos[0] + 1 == len(ntDel)
                except:
                    print >> sys.stdout, "\tIncorrect deletion data with %s on allele %s. Skipping variation." % (var, allele)
                    continue
                            
                try:
                    assert preBackbone_seq[ map_cyp[pos[0]] : map_cyp[pos[1]] + 1 ] == ntDel
                except:
                    print >> sys.stdout, "\tWarning, positions %d to %d in sequence contains %s, but expected %s from database" % \
                          (dbPos[0], dbPos[1], preBackbone_seq[ map_cyp[pos[0]] : map_cyp[pos[1]] + 1 ], ntDel)
                    print >> sys.stdout, "\t\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
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
                        pos[i] = pos[i] + oSetPos
                    else:
                        pos[i] = pos[i] + oSetNeg

                skipIns = False
                for i in range(len(pos)):
                    if pos[i] < 0 or pos[i] > len(cyp_seq) - 1:
                        print >> sys.stdout, "Warning: position %d out of bounds" % (dbPos[i])
                        print >> sys.stdout, "\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                        skipIns = True

                if skipIns:
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

    msfFile = open('cyp_msf/%s.msf' % gene_name,'w')
    for allele,msf_seq in msfTable.items():
        print >> msfFile, "%s\t%s" % (allele, msf_seq)

    msfFile.close()


def main():
    os.system('mkdir cyp_msf')

    cyp_gene_names = ['cyp1a1','cyp1a2','cyp1b1','cyp2a6',
                      'cyp2a13','cyp2b6','cyp2c8','cyp2c9',
                      'cyp2c19','cyp2d6','cyp2e1','cyp2f1',
                      'cyp2j2','cyp2r1','cyp2S1','cyp2w1',
                      'cyp3a4','cyp3a5','cyp3a7','cyp3a43',
                      'cyp4a11','cyp4a22','cyp4b1','cyp4f2',
                      'cyp5a1','cyp8a1','cyp19a1','cyp21a2',
                      'cyp26a1']

    oSetPos = 0
    oSetNeg = 0
    oSetScorePos = 0.0
    oSetScoreNeg = 0.0
        

    for gene_name in cyp_gene_names:
        oSetPos, oSetNeg, oSetScorePos, oSetScoreNeg = checkNTloc("cyp_fasta/%s.fasta" % gene_name,"cyp_var_files/%s.var" % gene_name,gene_name)
        if not (oSetScorePos > 0.95 and oSetScoreNeg > 0.95):
            print "Less than 95% bp match, skipping gene."
            continue
        
        makeMSF(gene_name, oSetPos, oSetNeg)

main()
