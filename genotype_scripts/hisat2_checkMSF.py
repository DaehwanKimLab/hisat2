#!/usr/bin/env python

import os, sys, subprocess, re
import inspect, operator
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

def checkMSFfile(gene_name, msf_fname, var_fname, fasta_filename):
    print("\n\nGene: %s" % gene_name)
    
    msf_file = open(msf_fname,'r')
    msf_dict = readMSF(msf_file)
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
    print(len(oriSeq))
    print(len(msf_dict[ref_allele].replace('.','')))

    '''print >> open('oriSeq.txt','w'), oriSeq
    print >> open('msfSeq.txt','w'), msf_dict[ref_allele].replace('.','')'''

    try:
        assert msf_dict[ref_allele].replace('.','') == oriSeq
        print("Sequences match for reference allele %s" % ref_allele)
    except AssertionError:
        print("Warning: sequences do not match for reference allele %s" % ref_allele)
    


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
    
##    cyp_var_file = open('cyp_var_files/%s.var' % gene_name,'r')
##    cyp_var_dict = makeVarDict(cyp_var_file)
##    cyp_var_file.close()
##
##    cyp_fasta = open('cyp_fasta/%s.fasta' % gene_name, 'r')
##    cyp_seq = extractSeq()
##    cyp_fasta.close()


    

main()
