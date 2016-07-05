#!/usr/bin/env python

import os, sys, subprocess, re
import inspect, operator
from argparse import ArgumentParser, FileType

incorrect_msf_entries = []

def checkNTloc(fasta_fileName,var_fileName,gene_name):
    print("\n\nGene: %s" % gene_name)
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
        
def create_inv_map(seq):
    seq_map = {}
    count = 0
    for i in range(len(seq)):
        bp = seq[i]
        if bp == '.':
            continue
        assert bp.upper() in "ACGT"
        seq_map[i] = count
        count += 1
    return seq_map

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

def readMSF(msf_fname):
    msf_dict = {}
    for line in msf_fname:
        line = line.strip()
        allele_name = line.split('\t')[0]
        msf_seq = line.split('\t')[1]
        assert not allele_name in msf_dict
        msf_dict[allele_name] = msf_seq

    return msf_dict

def msf_removeIns(ref_seq, al_seq):
    assert len(ref_seq) == len(al_seq)
    ins_ind_list = []
    for i in range(len(ref_seq)):
        if ref_seq[i] == '.':
            ins_ind_list.append(i)

    ori_ref_seq = ref_seq.replace('.','')
    ori_al_seq = list(al_seq)

    for i in ins_ind_list:
        ori_al_seq[i] = '-'

    ori_al_seq = ''.join(ori_al_seq).replace('-','')

    assert len(ori_ref_seq) == len(ori_al_seq)
    return ori_ref_seq, ori_al_seq        


def msfToVarList(ref_seq, al_seq):
    var_list = []
    
    assert len(ref_seq) == len(al_seq)
    for bp in ref_seq: assert bp in "ACGT."
    for bp in al_seq: assert bp in "ACGT."
    inv_map = create_inv_map(ref_seq)
    
    ins_re = re.compile('[ACGT]\.+')
    ins_subStrPos = [(m.start(0), m.end(0)) for m in re.finditer(ins_re, ref_seq)] # list of duples of indicies of insertions in ref_seq
    ins_pos_length = [(tup[0], tup[1] - tup[0] - 1) for tup in ins_subStrPos]

    for tup in ins_pos_length:
        ins_pos, ins_length = tup[0], tup[1]
        ins_seq = al_seq[ins_pos + 1: ins_pos + ins_length  + 1]
        ins_seq = ins_seq.replace('.','')
        if len(ins_seq) == 0:
            continue
        ins_str_data = str(inv_map[tup[0]]) + '_' + str(inv_map[tup[0]] + 1) + 'ins' + ins_seq
        var_list.append(ins_str_data)

    # insertions finished
    
    ori_ref_seq, ori_al_seq = msf_removeIns(ref_seq, al_seq)

    for i in range(len(ori_ref_seq)):
        if ori_al_seq[i] == '.':
            continue 
        elif ori_al_seq[i] != ori_ref_seq[i]: # snp
            var_list.append(str(i) + ori_ref_seq[i] + '>' + ori_al_seq[i])

    del_subStrPos = [(m.start(0), m.end(0)) for m in re.finditer(ins_re, ori_al_seq)] # list of duples of indicies of deletions in ori_al_seq
    del_pos_length = [(tup[0], tup[1] - tup[0] - 1) for tup in del_subStrPos]

    for tup in del_pos_length:
        del_pos, del_length = tup[0], tup[1]
        del_seq = ori_ref_seq[del_pos + 1 : del_pos + del_length + 1]
        if del_length == 1:
            assert len(del_seq) == 1
            del_str_data = str(tup[0] + 1) + 'del' + del_seq
        else:
            del_str_data = str(tup[0] + 1) + '_' + str(tup[0] + tup[1]) + 'del' + del_seq
        var_list.append(del_str_data)

    # deletions finished

    return var_list

def checkMSFfile(gene_name, msf_fname, var_fname, fasta_filename):
    oSetPos, oSetNeg, oSet_pos_score, oSet_neg_score = checkNTloc(fasta_filename, var_fname, gene_name)
    # print("\n\nGene: %s" % gene_name)
    
    msf_file = open(msf_fname,'r')
    msf_dict = readMSF(msf_file) # { Allele name : MSF sequence }
    msf_file.close()

    var_file = open(var_fname,'r')
    var_dict = makeVarDict(var_file)
    var_file.close()
    
    fa_file = open(fasta_filename,'r')
    oriSeq = extractSeq(fa_file)
    fa_file.close()


    # Find reference allele
    ref_allele = ''
    for allele_name in var_dict.keys():
        if len(var_dict[allele_name]) == 1 and list(var_dict[allele_name])[0] == "None":
            assert ref_allele == ''
            ref_allele = allele_name
    assert not ref_allele == ''


    # Check if ref allele seq in msf matches fasta
    assert ref_allele in msf_dict
    '''print(len(oriSeq))
    print(len(msf_dict[ref_allele].replace('.','')))'''

    '''print >> open('oriSeq.txt','w'), oriSeq
    print >> open('msfSeq.txt','w'), msf_dict[ref_allele].replace('.','')'''

    try:
        assert msf_dict[ref_allele].replace('.','') == oriSeq
        print("Sequences match for reference allele %s" % ref_allele)
    except AssertionError:
        print("Warning: sequences do not match for reference allele %s" % ref_allele)
        sys.exit(1)


    # Check all alleles are included
    try:
        assert set(msf_dict.keys()).issubset(set(var_dict.keys()))
    except AssertionError:
        print("Extra alleles in MSF!\n")
        print(sorted(msf_dict.keys()))
        print("\n\n")
        print(sorted(var_dict.keys()))
        sys.exit(1)


    # Convert from database positions to sequence positions (using offset)
    for allele, var_list in var_dict.items():
        oSet_var_list = []
        for var in var_list:
            if '>' in var: # snp
                pos = int(var.split('>')[0][:-1])
                ntSnp = [var.split('>')[0][-1]]
                ntSnp.append(var.split('>')[1])
                assert len(ntSnp) == 2
                if pos > 0:
                    pos = pos + oSetPos
                else:
                    pos = pos + oSetNeg

                if pos < 0 or pos > len(oriSeq) - 1: # out of bounds
                    continue
                if oriSeq[pos] != ntSnp[0]: # mismatch
                    print('\tMismatch on variation %s' % var)
                    continue

                oSet_var = str(pos) + ntSnp[0] + '>' + ntSnp[1]
                oSet_var_list.append(oSet_var)

            elif 'del' in var: # deletion
                pos = var.split('del')[0].split('_')
                pos = [int(p) for p in pos]
                if len(pos) == 1: # Handle single deletion with format for multi deletion with one location (e.g. [1707] -> [1707,1707])  
                    pos.append(pos[0])
                assert len(pos) == 2
                ntDel = var.split('del')[1]
                for nt in ntDel:
                    assert nt in "ACGT"

                skipDel = False
                for i in range(len(pos)):
                    if pos[i] > 0:
                        pos[i] = pos[i] + oSetPos
                    else:
                        pos[i] = pos[i] + oSetNeg
                    if pos[i] < 0 or pos[i] > len(oriSeq) - 1: # out of bounds
                        continue
                if (oriSeq[ pos[0] : pos[1] + 1 ] != ntDel): # mismatch
                    print('\tMismatch on variation %s' % var)
                    continue

                if skipDel:
                    continue

                assert pos[1] - pos[0] + 1 == len(ntDel)

                oSet_var = 'del' + ntDel
                if pos[0] == pos[1]:
                    oSet_var = str(pos[0]) + oSet_var
                else:
                    oSet_var = str(pos[0]) + '_' + str(pos[1]) + oSet_var

                oSet_var_list.append(oSet_var)                        

            elif 'ins' in var: # insertion
                pos = var.split('ins')[0].split('_')
                pos = [int(p) for p in pos]
                if len(pos) == 1:
                    pos.append(pos[0] + 1)
                assert len(pos) == 2
                assert pos[1] - pos[0] == 1
                ntIns = var.split('ins')[1]
                for nt in ntIns:
                    assert nt in "ACGT"

                skipIns = False
                for i in range(len(pos)):
                    if pos[i] > 0:
                        pos[i] = pos[i] + oSetPos
                    else:
                        pos[i] = pos[i] + oSetNeg
                    if pos[i] < 0 or pos[i] > len(oriSeq) - 1: # out of bounds
                        skipIns = True

                if skipIns:
                    continue

                oSet_var = str(pos[0]) + '_' + str(pos[1]) + 'ins' + ntIns
                oSet_var_list.append(oSet_var)

            else:
                assert allele == ref_allele
                assert var == 'None'
                assert len(oSet_var_list) == 0
                oSet_var_list.append('None')

        var_dict[allele] = set(oSet_var_list)



    # Check variants created from MSF file against variants list
    for allele, msf_seq in msf_dict.items():
        if allele == ref_allele:
            continue
        msf_var_list = msfToVarList(msf_dict[ref_allele], msf_seq)
        print('\t' + str(var_dict[allele] == set(msf_var_list)) + '\t' + str(allele) + '\t' + str(msf_var_list))

        try:
            assert var_dict[allele] == set(msf_var_list)
        except AssertionError:
            incorrect_msf_entries.append(allele)
            print('\n')
            print('\t\tVar File:\t' + str(var_dict[allele]))
            print('\t\tMSF File:\t' + str(set(msf_var_list)))
            print('\t\tDifference:\t' + str(var_dict[allele] - set(msf_var_list)) + '\n')
            '''sys.exit(1)'''

def main():
  
    gene_names = ['cyp1a1','cyp1a2','cyp1b1','cyp2a6',
                  'cyp2a13','cyp2b6','cyp2c8','cyp2c9',
                  'cyp2c19','cyp2d6','cyp2e1','cyp2f1',
                  'cyp2j2','cyp2r1','cyp2S1','cyp2w1',
                  'cyp3a4','cyp3a5','cyp3a7','cyp3a43',
                  'cyp4a11','cyp4a22','cyp4b1','cyp4f2',
                  'cyp5a1','cyp8a1','cyp19a1','cyp21a2',
                  'cyp26a1']

    for gene_name in gene_names:
        checkMSFfile(gene_name, 'cyp_msf/%s.msf' % gene_name, 'cyp_var_files/%s.var' % gene_name, 'cyp_fasta/%s.fasta' % gene_name)

    print('\n\n%d incorrect msf entries on alleles %s' % (len(incorrect_msf_entries), str(incorrect_msf_entries)))


main()

