#!/usr/bin/env python

import sys, os, subprocess, re 
import inspect, random 
import math 
import glob
import urllib2
import xml.etree.ElementTree as etree
from argparse import ArgumentParser, FileType

"""
Download RBC data from website 
"""

def get_xml(url):
    file = urllib2.urlopen(url)
    data = file.read()
    file.close()

    tree = etree.fromstring(data)
    return tree

def get_website(url):
    tries = 0
    while tries < 4:
        tries += 1
        try:
            # print >> sys.stderr, 'Attempt %d of 3 to connect to URL' % tries
            resp = urllib2.urlopen(url)
        except IOError:
            # print >> sys.stderr, 'No Connection to URL'
            if tries == 3:
                raise ValueError()
        else:
            # print >> sys.stderr, 'Connection successful!'
            file = resp.read()
            return file.splitlines(True)

def get_geneRefSeq(locus = {}):
    refseq = {}
    for key in locus:
        refseq.update({ key : [] })

    webaddress = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene'
    try:
        website = get_website(webaddress)
    except ValueError:
        print >> sys.stderr, 'Cannot access refseq database at this time'
        raise ValueError()

    for line in website:
        line = line.strip().split('\t')
        if line[2] in refseq:
            refseq[line[2]].append(line[3])
      
    return refseq

def get_seqbyRef(access, gene_name = '', getall = False):
    webaddress = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=gb&sendto=on&id=' + access
    try:
        website = get_website(webaddress)
    except ValueError:
        # print >> sys.stderr, 'No data available for accession: %s' % access
        raise ValueError()

    seqline = ''
    exon_ranges = {}
    mRNAline = 1
    gene_hit = False
    gene_found = False
    start_seq = False
    for line in website:
        if line.startswith('//'): # end of fine
            break
        
        # get sequence
        line = line.strip()
        if start_seq:
            line = line.replace(' ','')
            seqline += re.sub(r'\d+', '', line.upper())
        if line.startswith('ORIGIN'):
            start_seq = True

        # Get exons
        if getall and gene_name not in exon_ranges:
            assert gene_name != ''
            if line.startswith('gene'):
                gene_range = line.split(' ')[-1].split('..')
                gene_hit = True
            elif line.startswith('/gene') and gene_hit:
                if gene_name in line.replace('"', '').split('='):
                    left, right = int(gene_range[0]) - 1, int(gene_range[1])
                    gene_found = True
            elif line.startswith('mRNA') and gene_found:
                raw_exons = re.findall('\(([^)]+)', line)
                if raw_exons[0][-1] == ',':
                    mRNAline += 1
                if mRNAline == 1:
                    exon_ranges.update({ gene_name : ''.join(raw_exons).split(',') })
                continue
            if mRNAline > 1:
                raw_exons.append(line.replace(')', ''))
                if raw_exons[-1][-1] == ',':
                    mRNAline += 1
                else:
                    exon_ranges.update({ gene_name : ''.join(raw_exons).split(',') })

    if getall:
        for key, anExon in exon_ranges.items():
            corrected_exon = []
            for pos in anExon:
                pos = pos.split('..')
                corrected_exon.append([int(pos[0])-left, int(pos[1])-left])
            exon_ranges[key] = corrected_exon
        seqline = seqline[left:right]
        
        return seqline, exon_ranges
    else:
        return seqline

def rev_comp(seq):
    complement = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    return ''.join(complement[base] for base in re.sub('[^ATCG]', '', seq[::-1]))

def check_substr(str1, str2, k):
    if len(str1) < len(str2):
        a, b = str1, str2
    else:
        b, a = str1, str2
    
    rev_strand = False
    while True:
        for substr in (a[i:i+k] for i in range(len(a)-k+1)):
            if substr in b:
                return True
        if rev_strand:
            return False
        else:
            rev_strand = True
            a = rev_comp(a)
 
def match_seq(ref, seq):
    subseq = [ref[i:i+100] for i in range(0, len(ref), 100)]
    correct = False
    for bin in subseq:
        if 'N' in bin:
            continue
        if bin in seq:
            correct = True
            break
    
    if not correct:
        correct = check_substr(ref, seq, 50)
     
    return correct

def write_fasta(gene, alleles = {}):
    file_path = os.path.dirname('RBG/fasta/')
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    filename = file_path + '/%s_gen.fasta' % gene
    ofile = open(filename, 'w')
    
    for allele, seq in alleles.items():
        ofile.write('>' + allele + ' ' + str(len(seq)) + 'bp' + '\n')
        for subseq in (seq[i:i+60] for i in range(0, len(seq), 60)):
            ofile.write(subseq + '\n')

    ofile.close()

def extract_RBC():
    # RBC Database URL
    rbc_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Alleles/dbRBC_allelevFINAL.xml"

    raw_data = get_xml(rbc_url) 
    
    locus = {}
    genSeq = {}
    accession = {}

    for allele in raw_data.iter('allele'):
        # Update gene loci
        alleleGene = allele[1].text.replace('_', '*').split('*')[0]
        alleleNam = allele[1].text.replace(alleleGene, '')[1:]
         
        switcher = {
            'H' : 'FUT1',
            'SE' : 'FUT2',
            'LE' : 'FUT3',
            'K' : 'KEL',
            'DARC' : 'ACKR1',
            'LU' : 'BCAM',
            'YT' : 'ACHE',
            'C4B' : 'C4A',
            'RHCE' : 'RHD'}

        if alleleGene in switcher:
            if allele.find('locus').find('name').text == 'RG': # switch one C4B to C4A
                alleleGene = switcher[alleleGene]
            elif allele.find('locus').find('name').text == 'RHD': # switch a few RHCE to RHD
                alleleGene = switcher[alleleGene]
            else:
                alleleGene = switcher[alleleGene]
        elif alleleGene == 'SE32??':
            alleleNam = allele[2].text.split('*')[1]
            alleleGene = 'FUT2'

        alleleNam = alleleGene + '*' + alleleNam.replace('?', '').replace('*','_')
        
        if alleleGene not in locus:
            locus.update({ alleleGene : [] })
        
        if alleleNam not in locus[alleleGene]:
            locus[alleleGene].append(alleleNam)
            genSeq.update({ alleleNam : '' })
            accession.update({ alleleNam : [] })
            
            # Get Accession Numbers
            if 'NM_' in allele[1].text or 'XM_' in allele[1].text:
                accession[alleleNam].append('_'.join(allele[1].text.split('_')[-2:]))
            elif '_' in allele[1].text:            
                accession[alleleNam].append(allele[1].text.split('_')[-1])

            for info in allele.iter('allele_info'):
                if info.find('accession') is not None:
                    accession[alleleNam].append(info.find('accession').text.replace('ABO_', ''))
            
            # Get Nucleotide sequence
            for seq in allele.find('blocks').iter('block'):
                seqReg = seq.find('sequence').text.replace('*','').replace('.','')
                if not seqReg:
                    continue

                gen_left, gen_right = len(genSeq[alleleNam])+1, len(genSeq[alleleNam])+len(seqReg)
                assert gen_left <= gen_right
                genSeq[alleleNam] += seqReg
    
    # Remove Genes with less than 5 alleles
    for gene, allelelist in locus.items():
        if len(allelelist) < 5:
            del locus[gene]
            for name in allelelist:
                del genSeq[name]
                del accession[name]

    # Check and correct sequences in loaded data
    refGeneSeq = {}
    refGeneExon = {}
    refseq = get_geneRefSeq(locus)
    for key, value in refseq.items():
        if not value:
            continue
        try:
            seq_ofRef, exons_ofRef = get_seqbyRef(value[0].split('.')[0], key, getall = True)
        except ValueError:
            del refseq[key]
            continue
        refGeneSeq.update({ key : seq_ofRef })
        refGeneExon.update(exons_ofRef)
    
    for key, seq in genSeq.items():
        for loci, alleles in locus.items():
            if key in alleles:
                gene_name = loci
        
        assert gene_name in refGeneSeq
        
        try:
            correct = match_seq(refGeneSeq[gene_name], seq)
        except ValueError:
            print "Problem with %s" % key

        if not correct:            
            for accNumber in accession[key]:
                reffasta = get_seqbyRef(accNumber)
                try:
                    match = match_seq(reffasta, seq)    
                except ValueError:
                    print "Problem with %s" % key
        
                if not match:
                    print >> sys.stdout, 'Correcting %s' % key
                    genSeq.update({ key : reffasta })
    
    refgene_indata = False
    for gene, seq in refGeneSeq.items():
        for name in locus[gene]:
            if (seq in genSeq[name]) and ('01' in name):
                refGeneSeq[name] = refGeneSeq.pop(gene)
                refGeneExon[name] = refGeneExon.pop(gene)
                refgene_indata = True
                break
        if refgene_indata:
            refgene_indata = False
            continue
        else:
            newgene = gene + '*refSeq.01'
            refGeneSeq[newgene] = refGeneSeq.pop(gene)
            refGeneExon[newgene] = refGeneExon.pop(gene)
            genSeq.update({ newgene : refGeneSeq[newgene] })

    for gene, allelelist in locus.items():
        alleleseq = {}
        for allele in allelelist:
            alleleseq.update({ allele : genSeq[allele] })
        write_fasta(gene, alleleseq)

if __name__ == '__main__':

    extract_RBC()
