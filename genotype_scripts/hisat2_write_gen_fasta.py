#!/usr/bin/env python

import os, sys, subprocess, re
import inspect
from argparse import ArgumentParser, FileType

def readMSF(msf_fname):
    msf_dict = {}
    for line in msf_fname:
        line = line.strip()
        allele_name = line.split('\t')[0]
        msf_seq = line.split('\t')[1]
        assert not allele_name in msf_dict
        msf_dict[allele_name] = msf_seq

    return msf_dict

def writeGenFasta(gene_name, msf_fname, line_length):
    msf_file = open(msf_fname,'r')
    msf_seq_dict = readMSF(msf_file)
    msf_file.close()

    gen_fasta_file = open('gen_fasta/%s_gen.fasta' % gene_name, 'w')
    
    for allele, seq in msf_seq_dict.items():
        seq = seq.replace('.','')
        print >> gen_fasta_file, ('>' + allele + ' ' + str(len(seq)) + ' bp')
        seq_lines = [seq[i:i+line_length] for i in range(0, len(seq), line_length)]
        print >> gen_fasta_file, ('\n'.join(seq_lines))

    gen_fasta_file.close()
    print('%s_gen.fasta completed' % gene_name)

def main():
    os.system('mkdir gen_fasta')
    gene_names = ['cyp1a1','cyp1a2','cyp1b1','cyp2a6',
                  'cyp2a13','cyp2b6','cyp2c8','cyp2c9',
                  'cyp2c19','cyp2d6','cyp2e1','cyp2f1',
                  'cyp2j2','cyp2r1','cyp2S1','cyp2w1',
                  'cyp3a4','cyp3a5','cyp3a7','cyp3a43',
                  'cyp4a11','cyp4a22','cyp4b1','cyp4f2',
                  'cyp5a1','cyp8a1','cyp19a1','cyp21a2',
                  'cyp26a1']

    for gene_name in gene_names:
        writeGenFasta(gene_name, 'cyp_msf/%s.msf' % gene_name, 60)

main()
