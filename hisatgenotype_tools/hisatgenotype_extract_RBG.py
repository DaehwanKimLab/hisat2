#!/usr/bin/env python

import sys, os, subprocess, re
import urllib2
import string
import xml.etree.ElementTree as etree
from multiprocessing import Pool
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common

"""
Download RBC data from website 
"""

# Symbol for empty nucleotide/sequence
emptySeq = '*'
skip_seqs = ['AQP3', 'XG', 'KEL'] 
    # In here due to error with only one sequence in tmp fasta file during msf building

def get_xml(url):
    file = urllib2.urlopen(url)
    data = file.read()
    file.close()

    tree = etree.fromstring(data)
    return tree

def get_website(url, verbose = False):
    tries = 0
    while tries < 4:
        tries += 1
        try:
            #if verbose: print >> sys.stderr, 'Attempt %d of 3 to connect to URL' % tries
            resp = urllib2.urlopen(url)
        except IOError:
            #if verbose: print >> sys.stderr, 'No Connection to URL'
            if tries == 3:
                raise ValueError()
        else:
            #if verbose: print >> sys.stderr, 'Connection successful!'
            file = resp.read()
            return file.splitlines(True)

def get_RefSeqID(locus = {}, verbose = False):
    refseq = {}
    for gene in locus:
        refseq.update({ gene : [] })

    webaddress = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene'
    try:
        website = get_website(webaddress, verbose)
    except ValueError:
        if verbose: print >> sys.stderr, 'Cannot access refseq database at this time'
        raise ValueError()

    for line in website:
        line = line.strip().split('\t')
        if line[2] in refseq:
            refseq[line[2]].append(line[3])
      
    return refseq

def get_seqbyRef(access, gene_name = '', getall = False, verbose = False):
    webaddress = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=gb&sendto=on&id=' + access
    try:
        website = get_website(webaddress, verbose)
    except ValueError:
        if verbose: print >> sys.stderr, 'No data available for accession: %s' % access
        raise ValueError()

    seqline = ''
    exon_ranges, exon_numbers, exon_hit = {}, [], False
    gene_hit, gene_found, gene_end, cds_found = False, False, False, False
    start_seq = False
    mRNAline = 1

    for line in website:
        if line.startswith('//'): # end of fine
            break
        
        # get sequence
        line = line.strip()
        if start_seq:
            line = line.replace(' ','')
            seqline += re.sub(r'\d+', '', line.upper())
            continue
        if line.startswith('ORIGIN'):
            start_seq = True
            continue

        # Get exons
        if getall and not gene_end:
            assert gene_name != ''
               
            if line.startswith('gene'):
                if gene_found:
                    gene_end = True
                    continue
                gene_range = line.split(' ')[-1].replace('>','').replace('<','').split('..')
                gene_hit = True
            elif line.startswith('/gene') and gene_hit:
                if gene_name in line.replace('"', '').split('='):
                    left, right = int(gene_range[0]) - 1, int(gene_range[1])
                    gene_found = True
            elif (line.startswith('CDS') or line.startswith('mRNA')) and \
                  gene_found and not cds_found:
                if line.startswith('CDS'):
                    cds_found = True
                raw_exons = re.findall('\(([^)]+)', line)
                if raw_exons and raw_exons[0][-1] == ',':
                    mRNAline += 1
                if mRNAline == 1:
                    exon_ranges.update({ gene_name : ''.join(raw_exons).split(',') })
                continue
            elif line.startswith('exon') and gene_found:
                exon_hit = True
            elif line.find('number') != -1 and exon_hit:
                exon_numbers.append(line[-1])
                exon_hit = False
                
            if mRNAline > 1:
                raw_exons.append(line.replace(')', ''))
                if raw_exons[-1][-1] == ',':
                    mRNAline += 1
                else:
                    exon_ranges.update({ gene_name : ''.join(raw_exons).split(',') })
                    mRNAline = 1        

    if getall and gene_found:
        use_given = True if exon_numbers else False
        for key, anExon in exon_ranges.items():
            if not anExon:
                raise ValueError()
            corrected_exon = []
            for itr in range(len(anExon)):
                if len(exon_numbers) != len(anExon) and gene_name != 'ABO':
                    use_given = False
                exon_num = exon_numbers[itr] if use_given else itr+1
                pos = anExon[itr].replace('>','').replace('<','').split('..')
                corrected_exon.append([exon_num, int(pos[0])-left, int(pos[1])-left])
            exon_ranges[key] = corrected_exon
        seqline = seqline[left:right]
        
        return seqline, exon_ranges
    elif getall and not gene_found:
        if verbose: print 'No exons found for gene %s in accession %s' % (gene_name, access)
        return seqline, exon_ranges.update({ gene_name : [-1, -1] })
    else:
        return seqline

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
                a = typing_common.reverse_complement(a)

    if not correct:
        correct = check_substr(ref, seq, 50)
     
    return correct

def write_fasta(out_dir, gene, alleles = {}, loc = 'gen'):
    file_path = os.path.dirname('%s/fasta/' % out_dir)
    if not os.path.exists(file_path):
        os.makedirs(file_path)
 
    filename = file_path + '/%s_%s.fasta' % (gene, loc)
    ofile = open(filename, 'w')
    
    for allele, seq in alleles.items():
        seq = seq.replace(emptySeq, '').replace('.','').replace('N','')
        ofile.write('>' + allele + ' ' + str(len(seq)) + 'bp' + '\n')
        for subseq in (seq[i:i+60] for i in range(0, len(seq), 60)):
            ofile.write(subseq + '\n')

    ofile.close()
    return filename

def run_clustalo(in_file):
    out_file = in_file.replace('.fasta', '.msf')
    cmd = 'clustalo --in=%s --full --outfmt=msf --out=%s --threads=1 --force' % (in_file, out_file)
    os.system(cmd)

def run_clustalo_pile(in_file):
    out_file = in_file.replace('.fasta', '.msf')
    cmd = 'clustalo --in=%s --pileup --outfmt=msf --out=%s --threads=1 --force' % (in_file, out_file)
    # cmd = 'clustalw2 -infile=%s -output=GDE -outfile=%s -quiet -case=upper > /dev/null' % (in_file, out_file)
    os.system(cmd)

def write_msf(threads, out_dir, genfilename, nucfilename, gene, allele_key, alleles = {}, fullGene = [], exon_seg = [], verbose = False):
    filepath = '%s/msf/' % out_dir
    file_path = os.path.dirname(filepath)
    if not os.path.exists(file_path):
        os.makedirs(file_path)
 
    numsegs = 0
    files = []
    msf_aligned_gen_list = {}
    msf_aligned_gen = {}
    msf_aligned_nuc = {}
    for allele, segs in alleles.items():
        msf_aligned_nuc.update({ allele : '' })
        if fullGene[allele]:
            msf_aligned_gen_list.update({ allele : [] })
            msf_aligned_gen.update({ allele : '' })

        if numsegs == 0:
            numsegs = len(segs)
        else:
            assert numsegs == len(segs), 'Allele %s of length %d when %d is needed' % (allele, len(segs), numsegs)

        for itr in range(len(segs)):
            if not fullGene[allele] and itr not in exon_seg:
                continue

            seq = segs[itr].replace('.', '') # re.sub(r'((?:[ACTG]))\.((?:[ACTG]))', r'\1N\2', segs[itr]).replace('.', 'N')
            parts = [s for s in seq.split(emptySeq) if s]

            if not seq.replace(emptySeq, ''):
                continue
            
            file_name = filepath + 'tmp_%s_gen_%d.fasta' % (gene, itr)
            if file_name not in files and itr == (len(files)):
                files.append(file_name)

            if os.path.exists(file_name):
                append_write = 'a'
            else:
                append_write = 'w'

            f = open(file_name, append_write)
            if len(parts) > 1:
                for itr in range(len(parts)):
                    part =  '%s_pt%d' % (allele, itr)
                    allele_key.update({ re.sub('[^0-9a-zA-Z]+', '_', part) : allele })
                    f.write('> %s\n' % (part))
                    f.write('%s\n' % parts[itr])
            else:
                f.write('> %s\n' % allele)
                f.write('%s\n' % parts[0])
            f.close()

    pool = Pool(int(threads))
    for itr in range(len(files)):
        in_file = files[itr]
        if itr in exon_seg:
            pool.apply_async(run_clustalo_pile, args=(in_file,))
        else:
            pool.apply_async(run_clustalo, args=(in_file,))
    pool.close()
    pool.join()
    
    for itr in range(len(files)):
        file_ = files[itr]
        parts, fulls, filealleles = {}, {}, []
        seq = ''
        header = True
        for line in open(file_.replace('.fasta', '.msf'), 'r'):
            line = line.strip()
            if line.startswith('//'):
                header = False
                continue
            if not line or not line[0].isalpha():
                continue

            if header:
                if line.startswith('Name'):
                    _, msfname, _, length = line.split()[:4]
                    
                    if 'pt' in msfname:
                        parts.update({ msfname : '' })
                    else:
                        fulls.update({ allele_key[msfname] : '' })
                    
                    filealleles.append(allele_key[msfname])
                continue

            msfallele = line.split()[0]
            seq = ''.join(line.split()[1:])

            assert msfallele in allele_key, "%s not in msf2name" % msfallele

            if 'pt' in msfallele:
                parts[msfallele] += seq
            else:
                fulls[allele_key[msfallele]] += seq

        parts_ = {}
        for allele, seq in parts.items():
            if allele_key[allele] not in parts_:
                parts_.update({ allele_key[allele] : [] })
            
            parts_[allele_key[allele]].append(seq)

        for allele, seqs in parts_.items():
            merg = ['~'] * len(seqs[0])
            for seq in seqs:
                assert len(seq) == len(merg)
                for jtr in range(len(seq)):
                    if merg[jtr] == '~':
                        merg[jtr] = seq[jtr]
                    elif merg[jtr] != seq[jtr]:
                        assert seq[jtr] == '~'
            fulls.update({ allele : ''.join(merg) })
        
        # QC Check
        for allele in filealleles:
            assert allele in fulls, '%s allele not in fulls' % allele
        
        # Parsing Sequences based on full length bool
        length = 0
        for allele, seq in fulls.items():
            if length == 0:
                length = len(seq)
            else:           
                assert length == len(seq), '%s : %s : %d : %d' % (itr, allele, length, len(seq))

            if fullGene[allele]:
                msf_aligned_gen_list[allele].append(seq)
            if itr in exon_seg:
                msf_aligned_nuc[allele] += seq

        for allele in msf_aligned_gen_list:
            if allele in filealleles:
                continue
            msf_aligned_gen_list[allele].append('~'*int(length))

        if itr in exon_seg:
            for allele in msf_aligned_nuc:
                if allele in filealleles:
                    continue
                msf_aligned_nuc[allele] += ('~'*int(length))

        # Correct Boundries
        if itr != 0:
            boundries = {'nts' : 0,
                        '~' : 0,
                        'len': 0 }
            for allele, seqs in msf_aligned_gen_list.items():
                left, right = seqs[itr-1][-1], seqs[itr][0]
                if left == '~' and right == '~':
                    continue

                boundries['len'] += 1
                if left in ['A', 'C', 'G', 'T', 'N']:
                    boundries['nts'] += 1
                else:
                    assert left == '~'
                    boundries['~'] += 1 
            
                if right in ['A', 'C', 'G', 'T', 'N']:
                    boundries['nts'] += 1
                else:
                    assert right == '~'
                    boundries['~'] += 1 
        
            if boundries['nts'] == boundries['~'] == boundries['len']:
                for allele, seqs in msf_aligned_gen_list.items():
                    left, right = seqs[itr-1][-1], seqs[itr][0]
                    if left == '~':
                        seqs[itr-1] = seqs[itr-1][:-1]
                    elif right == '~':
                        seqs[itr] = seqs[itr][1:]
                    else:
                        if verbose: print 'Error in boundry correction'
                        exit(1)
            
    for allele, seqs in msf_aligned_gen_list.items():
        msf_aligned_gen[allele] = ''.join(seqs)
   
    for file_ in files:
        try:
            os.remove(file_)
            os.remove(file_.replace('.fasta', '.msf'))
            os.remove(file_.replace('.fasta', '.dnd'))
        except OSError:
            pass

    def write_MSF_file(filename, msf_aligned):
        # Sanity checking and deletion check
        seq_len = 0  
        for allele, msf_seq in msf_aligned.items():
            if seq_len == 0:
                seq_len = len(msf_seq)
            else:
                assert seq_len == len(msf_seq), '%s allele lengths are %d and %d' % (allele, seq_len, len(msf_seq))
        assert seq_len > 0

        file_out = file_path + '/' + filename.split('/')[-1].replace('.fasta', '.msf')
        ofile = open(file_out, 'w')

        # Follow MSF style of IMGT/HLA database
        for i in range(0, seq_len, 50):
            for allele, msf_seq in msf_aligned.items():
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

    write_MSF_file(genfilename, msf_aligned_gen)
    write_MSF_file(nucfilename, msf_aligned_nuc)

def extract_RBC(out_dir,
                threads,
                gene_list,
                verbose):
    # RBC Database URL
    if verbose: 
        print >> sys.stdout, 'Loading dbRBC from NCBI'
    else:
        sys.stdout.write("Processing: [")
        sys.stdout.flush()

    rbc_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Alleles/dbRBC_allelevFINAL.xml"

    raw_data = get_xml(rbc_url) 
    
    locus = {}
    fullGene = {}
    mixed_partials = []

    geneSegment = {}
    genSeq = {}
    genExon = {}

    accession = {}
    exon_counts = {}

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

        if gene_list:
            if alleleGene not in gene_list:
                continue

        if alleleNam in skip_seqs or alleleGene in skip_seqs:
            continue

        alleleNam = alleleGene + '*' + alleleNam.replace('?', '').replace('*','_')        
        
        if alleleGene not in locus:
            locus.update({ alleleGene : [] })
            exon_counts.update({ alleleGene : 0 })

        if alleleNam not in locus[alleleGene]:
            genesegs = []
            allele_access = []
            
            # Get Accession Numbers
            if 'NM_' in allele[1].text or 'XM_' in allele[1].text:
                allele_access.append('_'.join(allele[1].text.split('_')[-2:]))
            elif '_' in allele[1].text:            
                allele_access.append(allele[1].text.split('_')[-1])

            for info in allele.iter('allele_info'):
                if info.find('accession') is not None:
                    allele_access.append(info.find('accession').text.split('_')[-1])
            
            intron_count, empty_intron, exon_count, empty_exon = 0, 0, 0, 0
            # Get Nucleotide sequence
            for seqnode in allele.find('blocks').iter('block'):
                seq = seqnode.find('sequence').text.replace('*', emptySeq)
                region = seqnode.find('type').text

                alt_nuc = 'RYKMSWBDHV' # Make Translation Table 
                seqReg = seq.translate(string.maketrans(alt_nuc, '.' * len(alt_nuc)))                
                raw_seq = seqReg.replace(emptySeq, '').replace('.', '')
                if len(raw_seq) < 5:
                    raw_seq = ''
                    seqReg = re.sub(r'[ACGT]' , emptySeq, seqReg)              

                if region.startswith(('Intron')):
                    intron_count += 1
                    if seq.count(emptySeq) != 0:
                        empty_intron +=1

                if region.startswith(('Exon')):
                    exon_count += 1
                    if seq.count(emptySeq) == len(seq):
                        empty_exon +=1
                    
                genesegs.append([region , seqReg])

            if exon_counts[alleleGene] == 0:
                exon_counts[alleleGene] = exon_count
            elif exon_counts[alleleGene] != exon_count:
                if verbose: print "\t%s is inconsistent with other alleles: Removing" % alleleNam
                continue
            
            fullGene.update({ alleleNam : True })
            assert intron_count >= empty_intron
            if empty_intron == intron_count:
                fullGene[alleleNam] = False
            elif empty_exon > 0 and empty_intron > 0:
                mixed_partials.append(alleleNam)

            locus[alleleGene].append(alleleNam)
            geneSegment.update({ alleleNam : genesegs })
            genSeq.update({ alleleNam : '' })
            genExon.update({ alleleNam : [] })
            accession.update({ alleleNam : allele_access })
 
    # Loading and adding RefSeq gene to data
    if verbose: 
        print >> sys.stdout, 'Loading RefSeq of RBG Genes from NCBI'
    else:
        sys.stdout.write("====")
        sys.stdout.flush()

    refGeneSeq = {}
    refGeneExon = {}
    refseq = get_RefSeqID(locus, verbose)
    for gene, value in refseq.items():
        exon_num = 0

        if gene == 'HBD': # Temporary until refSeq fixed
            value = ['NG_000007.2']

        if not value:
            continue
        try:
            seq_ofRef, exons_ofRef = get_seqbyRef(value[0], gene, True, verbose)
        except ValueError:
            if verbose: print >> sys.stderr, '\tKey error %s not found: Removing' % gene
            del refseq[gene]
            for name in locus[gene]:
                del geneSegment[name]
                del accession[name]
            del locus[gene]
            continue

        # Remove genes with non-matching exon numbers
        exon_num = len(exons_ofRef[gene])
        if exon_num != exon_counts[gene]:
            if verbose: print "\tRemoving %s from data: Non-matching exon count to reference" % gene
            for name in locus[gene]:
                del geneSegment[name]
                del accession[name]
            del locus[gene]
            continue

        refGeneSeq.update({ gene : seq_ofRef })
        refGeneExon.update(exons_ofRef)
    
    for gene, seq in refGeneSeq.items():
        newgene = gene + '*refSeq.01'
        if verbose: print >> sys.stdout, '\tAdding %s to data' % newgene

        genSeq.update({ newgene : refGeneSeq[gene] })
        genExon.update({ newgene : refGeneExon[gene] })
        fullGene.update({ newgene : True })
        locus[gene].append(newgene)

    # Remove Genes with less than 5 alleles and trim 5' end
    for gene, allelelist in locus.items():
        if verbose: print "Trimming %s 5' end" % gene
        if len(allelelist) < 5:
            del locus[gene]
            for name in allelelist:
                if name not in geneSegment or name not in accession:
                    continue
                del geneSegment[name]
                del accession[name]
            continue
        
        """ 5' trimming in "promoter" region SAVE for later
        base_counts = []
        for allele in allelelist:
            five_prime = geneSegment[allele][0][1]
            if not base_counts:
                base_counts = [0 for i in range(len(five_prime))]
            else:
                assert len(five_prime) == len(base_counts)
            
            for itr in range(len(five_prime)):
                if five_prime[itr] != emptySeq:
                    base_counts[itr] += 1

        for itr in range(len(base_counts)):
            if base_counts[itr] == 1:
                continue
            cut_pos = itr
            break
        """

        cut_pos = int(refGeneExon[gene][0][1]) - 1

        for allele in allelelist:
            if allele not in geneSegment:
                continue
            
            geneSegment[allele][0][1] = geneSegment[allele][0][1][-cut_pos:]
            exon_count = 0
            for itr in range(len(geneSegment[allele])): 
                raw_seq = geneSegment[allele][itr][1].replace(emptySeq, '').replace('.', '')
                if not raw_seq:
                    continue
                
                left, right = len(genSeq[allele]) + 1, len(genSeq[allele]) + len(raw_seq)
                genSeq[allele] += raw_seq
                
                if geneSegment[allele][itr][0] == 'Exon':
                    exon_count += 1
                    genExon[allele].append([exon_count, left, right])
                     
    # extract HG38 genes
    if verbose: 
        print >> sys.stdout, 'Extracting Genes from HG38'
    else:
        sys.stdout.write("====")
        sys.stdout.flush()

    cigar_re = re.compile('\d+\w')

    for gene, ref in refGeneSeq.items():
        if gene in ['CR1', 'C4B', 'C4A', 'RHCE', 'RHD', 'GCNT2']:
            continue
        if gene not in locus:
            continue
        
        hg_allele = gene + '*HG38.ref.01'

        skip = False
        genestart = sys.maxint
        geneend = 0
        genechr = ''
        exons = []

        if verbose: print '\tExtracting %s from Genome' % gene
        cmd_aligner = ['hisat2', '-x', 'grch38/genome', '-f', '-c']
        
        cds_start, cds_end, upstream, downstream = sys.maxint, 0, 0, 0
        for pos in refGeneExon[gene]:
            _, exleft, exright = pos
            
            if cds_start > exleft:
                cds_start = exleft
            if cds_end < exright:
                cds_end = exright
            
            exon = ref[exleft - 1 : exright]
            exons.append(exon)
            
        cmd_aligner += ['%s' % ','.join(exons)]
        upstream = cds_start - 1
        downstream = len(ref) - cds_end

        align_proc = subprocess.Popen(cmd_aligner,
                                      stdout=subprocess.PIPE,
                                      stderr=open("/dev/null", 'w'))

        strand_reverse = False
        exons = []
        for line in align_proc.stdout:
            if line.startswith('@'):
                continue
            line = line.strip()
            cols = line.split()
            _, flag, chr, left, _, cigar_str = cols[:6]
            flag = int(flag)
            if flag & 0x10 != 0:
                strand_reverse = True
            if flag & 0x100 != 0:
                if verbose: print '\tWarning!: Gene %s with flag %d' % (gene, flag)
                continue
            left = int(left) - 1
            right = left

            cigars = cigar_re.findall(cigar_str)
            cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
            for i in range(len(cigars)):
                cigar_op, length = cigars[i]
                if cigar_op in "MND":
                    right += length
            
            if left < genestart:
                genestart = left
            if right > geneend:
                geneend = right

            if not genechr:
                genechr = chr
            elif genechr != chr and not skip:
                if verbose: print "\tError: exons in Gene %s aligning to chr %s and %s: Not adding" % (gene, genechr, chr)
                skip = True
            
            exons.append([left, right])

        align_proc.communicate()

        if skip:
            continue

        genSeq.update({ hg_allele : '' })
        genExon.update({ hg_allele : [] })
        fullGene.update({ hg_allele : True })
        locus[gene].append(hg_allele)

        assert downstream not in [1, -1]
        if strand_reverse:
            genestart -= downstream
            geneend += upstream
        else:
            genestart -= upstream
            geneend += downstream

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

        exon_pos = range(len(genSeq[hg_allele]))[::-1]
        hgexons = []
        for exon in exons:
            left = (exon[0] - genestart)
            right = (exon[1] - genestart)
            if strand_reverse:
                nleft = exon_pos[right] + 1
                nright = exon_pos[left]
            else:
                nleft = left + 1
                nright = right
            hgexons.append([nleft, nright])
        
        hgexons = sorted(hgexons, key = lambda x: x[0])

        for itr in range(len(hgexons)):
            genExon[hg_allele].append([itr+1, hgexons[itr][0], hgexons[itr][1]])
        
        if strand_reverse:            
            genSeq[hg_allele] = typing_common.reverse_complement(genSeq[hg_allele])

    # Check and correct sequences in loaded data
    if verbose: 
        print >> sys.stdout, 'Validating Alleles' 
    else:
        sys.stdout.write("====")
        sys.stdout.flush()

    for allele, seq in genSeq.items():
        for loci, alleles in locus.items():
            if allele in alleles:
                gene_name = loci
        
        if gene_name not in refGeneSeq:
            continue
        
        try:
            correct = match_seq(refGeneSeq[gene_name], seq)
        except ValueError:
            if verbose: print "Problem with %s" % allele

        if not correct and allele in accession:
            for accNumber in accession[allele]:
                try:
                    reffasta, refexon = get_seqbyRef(accNumber, gene_name, True, verbose)                  
                    match = match_seq(reffasta, seq)    
                except ValueError:
                    if verbose: print "Problem with %s" % allele
                    continue
       
                if refexon is None or gene_name not in refexon:
                    continue
 
                if not match:
                    if verbose: print >> sys.stdout, 'Correcting %s' % allele
                    genSeq[allele] = reffasta
                    genExon[allele] = refexon[gene_name]

                    boundries = [bound for exon in genExon[allele] for bound in exon[1:]]
                    startex, startintron = genExon[allele][0][0], False
                    
                    if boundries[0] != 1:
                        startintron = True 
                    else:
                        boundries.pop(0)

                    prev, segs = 0, []
                    for bound in boundries:
                        segs.append(['', reffasta[prev:bound-1]])
                        prev = bound-1
                    
                    exon_count = 0
                    for itr in range(len(geneSegment[allele])):
                        if 'Exon' in geneSegment[allele][itr]:
                            exon_count += 1
                        if exon_count == int(startex):
                            if startintron:
                                geneSegment[allele][itr-1] = segs[0]
                                segs.pop(0)
                            for seg in segs:
                                geneSegment[allele][itr] = seg
                                itr += 1
                            break

                    sum_seqlen = 0
                    for i in genExon[allele]:
                        sum_seqlen += (int(i[2])-int(i[1])-1)
                    if sum_seqlen < len(reffasta):
                        fullGene[allele] = True
                    else:
                        fullGene[allele] = False

    filepath = '%s/RBG/' % out_dir
    file_path = os.path.dirname(filepath)
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    ofile = open('%s/rbg.dat' % file_path, 'w')
    ofile.write('Red Blood Group gene allele database: V.1.0\n')
    ofile.close

    # Remove any miss aligned data by length and write genes out
    if not verbose:
        sys.stdout.write("====")
        sys.stdout.flush()

    for gene, allelelist in locus.items():
        if verbose: print >> sys.stdout, 'Processing %s:' % gene
        allele_index = {}
        msf_key = {}
        seq_list = []
        msf_alleleseq = {}
        exon_seg = []
        
        for itr in range(len(allelelist)):
            allele = allelelist[itr]
            msf_key.update({ re.sub('[^0-9a-zA-Z]+', '_', allele) : allele })
            allele_index.update({ allele : itr })
            seq_list.append(genSeq[allele])

        if verbose: print >> sys.stdout, '\t Checking %s for redundancy' % gene
        allele_index, seq_list = typing_common.collapse_alleles(allele_index, seq_list, emptySeq, False, verbose)

        full_allele = {}
        partial_allele = {}
        for allele, ind in allele_index.items():
            seq = seq_list[ind]
            ofile = open('RBG/rbg.dat', 'a')
            ofile.write('DE\t%s\n' % allele)
            for exon in genExon[allele]:
                ofile.write('FT\texon\t%s\n' % '..'.join(str(x) for x in exon[1:]))
                ofile.write('FT\t\t/number=%s\n' % exon[0])
            ofile.write('\\\\\n')
            ofile.close()

            partial_allele.update({ allele : '' })
            
            for exon in genExon[allele]:
                _, left, right = exon
                partial_allele[allele] += seq[left-1 : right]

            msf_alleleseq.update({ allele : [] })

            if allele in geneSegment:
                for itr in range(len(geneSegment[allele])):
                    sequ = geneSegment[allele][itr]
                    msf_alleleseq[allele].append(sequ[1])
                    if 'Exon' in sequ[0] and itr not in exon_seg:
                        exon_seg.append(itr)

            else:
                boundries = []
                for exon in genExon[allele]:
                    boundries.append(exon[1])
                    boundries.append(exon[2] + 1)

                prev = 0
                for bound in boundries:
                    if bound == 1:
                        msf_alleleseq[allele].append('')
                        continue
                    msf_alleleseq[allele].append(genSeq[allele][prev:bound-1])
                    prev = bound-1
                msf_alleleseq[allele].append(genSeq[allele][prev:len(genSeq[allele])]) 

            if fullGene[allele]:
                full_allele.update({ allele : seq })

        if verbose: print >> sys.stdout, '\t Writing %s Fasta' % gene
        fastaname_gen = write_fasta(filepath, gene, full_allele)
        fastaname_nuc = write_fasta(filepath, gene, partial_allele, loc = 'nuc')

        if verbose: print >> sys.stdout, '\t Writing %s MSF' % gene
         
        write_msf(threads,
                  filepath,
                  fastaname_gen, 
                  fastaname_nuc, 
                  gene,
                  msf_key, 
                  msf_alleleseq, 
                  fullGene, 
                  exon_seg,
                  verbose) 

    if not verbose:
        sys.stdout.write("====] Complete!\n")
        sys.stdout.flush()
         
if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract Red Blood Group Genes')
    parser.add_argument("--out-dir",
                        dest="out_dir",
                        type=str,
                        default=".",
                        help="Directory name for extracted read files")
    parser.add_argument("--gene-list",
                        dest="gene_list",
                        type=str,
                        default="",
                        help="A comma-separated list of genes (default: empty)")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="Verbose Output")

    args = parser.parse_args()

    if args.gene_list != '':
        gene_list = args.gene_list.upper().split(',')

    extract_RBC(args.out_dir,
                args.threads,
                gene_list,
                args.verbose)
