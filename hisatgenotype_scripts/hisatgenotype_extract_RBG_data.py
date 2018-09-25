#!/usr/bin/env python

import sys, os, subprocess, re 
import inspect, random 
import math 
import glob
import urllib2
import string
import xml.etree.ElementTree as etree
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align.Applications import ClustalwCommandline
from argparse import ArgumentParser, FileType

"""
Download RBC data from website 
"""

# Symbol for empty nucleotide/sequence
emptySeq = '*'

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
        print >> sys.stderr, 'No data available for accession: %s' % access
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
                gene_range = line.split(' ')[-1].replace('>','').replace('<','').split('..')
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

    if getall and gene_found:
        for key, anExon in exon_ranges.items():
            corrected_exon = []
            for pos in anExon:
                pos = pos.split('..') if not any(i in pos for i in r"<>") else pos.replace('>','').replace('<','').split('..')
                corrected_exon.append([int(pos[0])-left, int(pos[1])-left])
            exon_ranges[key] = corrected_exon
        seqline = seqline[left:right]
        
        return seqline, exon_ranges
    elif getall and not gene_found:
        print 'No exons found for gene %s in accession %s' % (gene_name, access)
        return seqline, exon_ranges.update({ gene_name : [-1, -1] })
    else:
        return seqline

def rev_comp(seq):
    complement = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    return ''.join(complement[base] for base in re.sub('[^ATCG]', '', seq[::-1]))

def match_seq(ref, seq):
    subseq = [ref[i:i+100] for i in range(0, len(ref), 100)]
    correct = False
    for bin in subseq:
        if 'N' in bin:
            continue
        if bin in seq:
            correct = True
            break

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

    if not correct:
        correct = check_substr(ref, seq, 50)
     
    return correct

def collapse_alleles(alleles = {}):
    remove = []
    for allele_i, seq_i in alleles.items():
        seq_i_strip = seq_i.replace(emptySeq, '').replace('.','')
        
        for allele_j, seq_j in alleles.items():
            seq_j_strip = seq_j.replace(emptySeq, '').replace('.','')
            if allele_i == allele_j:
                continue
            
            if seq_i == seq_j and allele_i > allele_j:
                if len(allele_i) <= len(allele_j):
                    print '\t\t %s is %s : Removing' % (allele_i, allele_j)
                    remove.append(allele_i)
                else:
                    print '\t\t %s is %s : Removing' % (allele_j, allele_i)
                    remove.append(allele_j)
                break

            if len(seq_i_strip) < len(seq_j_strip):
                if seq_i_strip in seq_j_strip:
                    if allele_i not in remove:
                        if 'HG19.ref' in allele_i:
                            print '\t\t Collapsing %s into %s' % (allele_i, allele_j)
                            remove.append(allele_i)
                        elif ('refSeq' in allele_j) or (('refSeq' in allele_i) and ('.' not in allele_j)):
                            print '\t\t Collapsing %s into %s' % (allele_j, allele_i)
                            remove.append(allele_j)
                        elif 'exon' in allele_i:
                            print '\t\t Collapsing %s into %s' % (allele_i, allele_j)
                            remove.append(allele_i)                            
                        #else:
                        #    print '\t\t Collapsing %s into %s' % (allele_i, allele_j)
                        #    remove.append(allele_i)
                    break
    
    for dup in remove:
        if dup in alleles:
            del alleles[dup]

    return alleles

def write_fasta(gene, alleles = {}, loc = 'gen'):
    file_path = os.path.dirname('RBG/fasta/')
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    filename = file_path + '/%s_%s.fasta' % (gene, loc)
    ofile = open(filename, 'w')
    
    for allele, seq in alleles.items():
        seq = seq.replace(emptySeq, '').replace('.','')
        ofile.write('>' + allele + ' ' + str(len(seq)) + 'bp' + '\n')
        for subseq in (seq[i:i+60] for i in range(0, len(seq), 60)):
            ofile.write(subseq + '\n')

    ofile.close()
    return filename

def write_msf(filename, gene, alleles = {}):
    file_path = os.path.dirname('RBG/msf/')
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    # Sanity checking
    seq_len = 0
    for allele, msf_seq in alleles.items():
        if seq_len == 0:
            seq_len = len(msf_seq)
        else:
            assert seq_len == len(msf_seq), '%s gene lengths are %d and %d' % (gene, seq_len, len(msf_seq))
    assert seq_len > 0

    file_out = file_path + '/' + filename.split('/')[-1].replace('.fasta', '.msf')
    ofile = open(file_out, 'w')

    # Follow MSF style of IMGT/HLA database
    for i in range(0, seq_len, 50):
        for allele, msf_seq in alleles.items():
            output = "%12s" % allele
            for j in range(i, i+50, 10):
                if j >= seq_len:
                    break
                if j == i:
                    output += "\t"
                else:
                    output += " "
                output += msf_seq[j:j+10]
            ofile.write(output + '\n')
        ofile.write('\n')

    ofile.close()

def extract_RBC():
    # RBC Database URL
    print >> sys.stdout, 'Loading dbRBC from NCBI'
    rbc_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Alleles/dbRBC_allelevFINAL.xml"

    raw_data = get_xml(rbc_url) 
    
    locus = {}
    fullGene = {}
    mixed_partials = []

    geneSegment = {}
    genSeq = {}
    genExon = {}

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

        if alleleGene in switcher.keys():
            if allele.find('locus').find('name').text == 'RG': # switch one C4B to C4A
                alleleGene = switcher[alleleGene]
            elif allele.find('locus').find('name').text == 'RHD': # switch a few RHCE to RHD
                alleleGene = switcher[alleleGene]
            elif alleleGene != 'RHCE' and alleleGene != 'C4B':
                alleleGene = switcher[alleleGene]
        elif alleleGene == 'SE32??':
            alleleNam = allele[2].text.split('*')[1]
            alleleGene = 'FUT2'

        alleleNam = alleleGene + '*' + alleleNam.replace('?', '').replace('*','_')
        
        if alleleGene not in locus:
            locus.update({ alleleGene : [] })
        
        if alleleNam not in locus[alleleGene]:
            locus[alleleGene].append(alleleNam)
            fullGene.update({ alleleNam : True })
            
            geneSegment.update({ alleleNam : [] })

            accession.update({ alleleNam : [] })
            
            # Get Accession Numbers
            if 'NM_' in allele[1].text or 'XM_' in allele[1].text:
                accession[alleleNam].append('_'.join(allele[1].text.split('_')[-2:]))
            elif '_' in allele[1].text:            
                accession[alleleNam].append(allele[1].text.split('_')[-1])

            for info in allele.iter('allele_info'):
                if info.find('accession') is not None:
                    accession[alleleNam].append(info.find('accession').text.replace('ABO_', ''))
            
            intron_count, empty_intron, exon_count, empty_exon = 0, 0, 0, 0
            # Get Nucleotide sequence
            for seqnode in allele.find('blocks').iter('block'):
                seq = seqnode.find('sequence').text
                region = seqnode.find('type').text

                if region.startswith(('Intron')):
                    intron_count += 1
                    if (seq.count('*')/len(seq)) > 0.5:
                        empty_intron +=1

                if region.startswith(('Exon')):
                    exon_count += 1
                    if (seq.count('*')/len(seq)) > 0.5:
                        empty_exon +=1

                seq = seq.replace(emptySeq, '').replace('.', '')
                alt_nuc = 'RYKMSWBDHV*' # Make Translation Table 
                seqReg = seq.translate(string.maketrans(alt_nuc, emptySeq * len(alt_nuc)))
                geneSegment[alleleNam].append([region , seqReg])

            if (empty_intron == intron_count) and (intron_count != 0):
                fullGene[alleleNam] = False
            elif empty_exon > 0 and empty_intron > 0:
                mixed_partials.append(alleleNam)

    # Remove Genes with less than 5 alleles
    for gene, allelelist in locus.items():
        if len(allelelist) < 5:
            del locus[gene]
            for name in allelelist:
                del geneSegment[name]
                del accession[name]

    # Loading and adding RefSeq gene to data
    print >> sys.stdout, 'Loading RefSeq of RBG Genes from NCBI'
    refGeneSeq = {}
    refGeneExon = {}
    refseq = get_geneRefSeq(locus)
    for key, value in refseq.items():
        if not value:
            continue
        try:
            seq_ofRef, exons_ofRef = get_seqbyRef(value[0].split('.')[0], key, getall = True)
        except ValueError:
            print >> sys.stderr, 'Key error %s not found' % key
            del refseq[key]
            continue
        refGeneSeq.update({ key : seq_ofRef })
        refGeneExon.update(exons_ofRef)
    
    refgene_indata = False
    refGeneKey = {}
    for gene, seq in refGeneSeq.items():
        for allele in locus[gene]:
            gen_seq = ''.join([reg[1] for reg in geneSegment[allele]])
            if (seq in gen_seq) and ('01' in allele):
                refGeneKey.update({ gene : allele })
                refgene_indata = True
                break

        if refgene_indata:
            refgene_indata = False
            continue
        else:
            newgene = gene + '*refSeq.01'
            refGeneKey.update({ gene : newgene })
            
            print >> sys.stdout, '\tAdding %s to data' % newgene
            genSeq.update({ newgene : refGeneSeq[gene] })
            genExon.update({ newgene : refGeneExon[gene] })
            fullGene.update({ newgene : True })
            locus[gene].append(newgene)
    
    # extract HG38 genes
    print >> sys.stdout, 'Extracting Genes from HG38'
    cigar_re = re.compile('\d+\w')

    for gene, ref in refGeneSeq.items():
        if gene in ['CR1', 'C4B', 'C4A']:
            continue
        
        hg_allele = gene + '*HG38.ref.01'

        genSeq.update({ hg_allele : '' })
        genExon.update({ hg_allele : [] })
        fullGene.update({ hg_allele : True })
        locus[gene].append(hg_allele)

        genestart = sys.maxint
        geneend = 0
        genechr = ''
        exons = []

        print '\tExtracting %s from Genome' % gene
        strand_reverse = False     
        for pos in refGeneExon[gene]:
            cmd_aligner = ['hisat2', '-x', 'grch38/genome', '-f', '-c']
            exleft, exright = pos
            exon = ref[exleft - 1 : exright]
            
            cmd_aligner += ['%s' % exon]
            align_proc = subprocess.Popen(cmd_aligner,
                                          stdout=subprocess.PIPE,
                                          stderr=open("/dev/null", 'w'))

            best_chr, best_left, best_right, best_AS, best_strand = "", -1, -1, -sys.maxint, ''
            for line in align_proc.stdout:
                if line.startswith('@'):
                    continue
                line = line.strip()
                cols = line.split()
                allele_id, flag, chr, left, mapQ, cigar_str = cols[:6]
                flag = int(flag)
                if flag & 0x10 != 0:
                    strand_reverse = True
                left = int(left) - 1
                right = left

                cigars = cigar_re.findall(cigar_str)
                cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
                if len(cigars) > 1 or cigars[0][0] != 'M':
                    continue
                for i in range(len(cigars)):
                    cigar_op, length = cigars[i]
                    if cigar_op in "MND":
                        right += length
                AS = ""
                for i in range(11, len(cols)):
                    col = cols[i]
                    if col.startswith("AS"):
                        AS = col[5:]
                assert AS != ""
                AS = int(AS)
                if AS > best_AS and (not genechr or chr == genechr):
                    best_chr, best_left, best_right, best_AS = chr, left, right, AS
 
                if best_left < genestart:
                    genestart = best_left
                if best_right > geneend:
                    geneend = best_right

                if not genechr:
                    genechr = best_chr
                assert genechr == best_chr, "Error: exons in Gene %s aligning to chr %s and %s" % (gene, genechr, chr)

                exons.append([best_left, best_right])
            align_proc.communicate()
 
        cmd_extractGene = ['samtools', 'faidx', 'genome.fa', '%s:%d-%d' % (genechr, genestart, geneend)]
        extract_proc = subprocess.Popen(cmd_extractGene,
                                        stdout=subprocess.PIPE,
                                        stderr=open("/dev/null", 'w'))
        for line in extract_proc.stdout:
            line = line.strip()
            if line.startswith('>'):
                continue
            genSeq[hg_allele] += line
        extract_proc.communicate()

        exon_pos = range(len(genSeq[hg_allele])+1)[::-1]
        for exon in exons:
            left = (exon[0] - genestart) + 1
            right = (exon[1] - genestart) + 1
            if strand_reverse:
                nleft = exon_pos[right-1]
                nright = exon_pos[left-1]
            else:
                nleft = left
                nright = right
            genExon[hg_allele].append([nleft, nright])
        
        genExon[hg_allele] = sorted(genExon[hg_allele], key = lambda x: x[0])
        print genExon[hg_allele]
        
        if strand_reverse:            
            genSeq[hg_allele] = rev_comp(genSeq[hg_allele])

    refgeneSegments = {}
    for gene, exons in genExon.items():
        if 'HG38' not in gene:
            continue
        refgeneSegments.update({ gene : [] })
        boundries = [bound for exon in exons for bound in exon]
        prev = 0
        for bound in boundries:
            refgeneSegments[gene].append(genSeq[gene][prev:bound-1])
            prev = bound-1
        print '%s\t:\t%d' % (gene, len(refgeneSegments[gene])) 

    for allele, regions in geneSegment.items():
        if not fullGene[allele]:
            continue
        seq_started = False
        seg_count = 0
        gene, _ = allele.split('*')
        if gene in ['CR1', 'C4B', 'C4A']:
            continue

        print '%s\t:\t%d' % (gene, len(regions))
        '''
        gene += '*HG38.ref.01'
        # assert len(regions) == len(refgeneSegments[gene])+1, 'Allele %s and Gene %s have different lengths' % (allele, gene)
        for region in regions:
            if region[1]:
                seq_started = True
            if seq_started and not region[1]:
                if seg_count > len(refgeneSegments[gene]):
                    break
                print 'Filling in %s from %s' % (allele, gene)
                geneSegment[allele][seg_count][1] = refgeneSegments[gene][seg_count]

            seg_count += 1
        '''

    exit(1)
    # Check and correct sequences in loaded data
    print >> sys.stdout, 'Validating Alleles'    
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
            """ 
            TODO: Implement allele handling by comparison with other alleles of gene
            Right now there are ~5 alleles that don't match the reference. Individual BLAST shows all are annotated correctly.
            These are not checked or handled here. If the database changes this could be a problem
            RHCE*AM295503; GYPB*AY509883; RHCE*AM295499; RHD*AJ867388; CD55*S72858
            """
            for accNumber in accession[key]:
                reffasta, refexon = get_seqbyRef(accNumber, gene_name, getall = True)
                try:
                    match = match_seq(reffasta, seq)    
                except ValueError:
                    print "Problem with %s" % key
        
                if not match and gene_name in refGeneKey:
                    print >> sys.stdout, 'Correcting %s' % key
                    genSeq[key] = reffasta
                    genExon[key] = refexon[gene_name]

                    sum_seqlen = 0
                    for i in genExon[key]:
                        sum_seqlen += (int(i[2])-int(i[1])-1)
                    if sum_seqlen < len(reffasta):
                        fullGene[key] = True
                    else:
                        fullGene[key] = False

    # Remove any miss aligned data by length and write genes out
    for gene, allelelist in locus.items():
        print >> sys.stdout, 'Processing %s:' % gene
        alleleseq = {}

        for allele in allelelist:
            alleleseq.update({ allele : genSeq[allele].replace(emptySeq, '').replace('.', '') })
             
        print >> sys.stdout, '\t Checking %s for redundancy' % gene
        alleleseq = collapse_alleles(alleleseq)
        
        full_allele = {}
        partial_allele = {}
        for allele, seq in alleleseq.items():
            partial_allele.update({ allele : '' })

            if gene == 'ABO':
                print genExon[allele]

            for exon in genExon[allele]:
                left, right = exon
                partial_allele[allele] += seq[left-1 : right]

            if fullGene[allele]:
                full_allele.update({ allele : seq })         
                
 
        #print >> sys.stdout, '\t Writing %s Fasta' % gene
        #fastaname_gen = write_fasta(gene, full_allele)
        #fastaname_nuc = write_fasta(gene, partial_allele, loc = 'nuc')

        # print >> sys.stdout, '\t Writing %s MSF' % gene
        # write_msf(fastaname, gene, alleleseq)   

if __name__ == '__main__':

    extract_RBC()
