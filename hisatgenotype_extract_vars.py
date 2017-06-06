#!/usr/bin/env python

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#


import os, sys, subprocess, re
import inspect
import glob
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common, hisatgenotype_gene_typing as gene_typing


"""
Mapping from base pair to a location in MSF format
"""
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


"""
"""
def create_consensus_seq(seqs,
                         seq_len,
                         min_var_freq,
                         remove_empty = True):
    consensus_freq = [[0, 0, 0, 0, 0] for i in range(seq_len)]
    for i in range(len(seqs)):                
        seq = seqs[i]
        if len(seq) != seq_len:
            continue                    
        for j in range(seq_len):
            nt = seq[j]
            assert nt in "ACGT.E"
            if nt == 'A':
                consensus_freq[j][0] += 1
            elif nt == 'C':
                consensus_freq[j][1] += 1
            elif nt == 'G':
                consensus_freq[j][2] += 1
            elif nt == 'T':
                consensus_freq[j][3] += 1
            else:
                assert nt in ".E"
                consensus_freq[j][4] += 1

    for j in range(len(consensus_freq)):
        for k in range(len(consensus_freq[j])):
            consensus_freq[j][k] /= float(len(seqs))
            consensus_freq[j][k] *= 100.0

    consensus_seq = ""
    has_empty = False
    for c in range(len(consensus_freq)):
        freq = consensus_freq[c]
        A, C, G, T, E = freq
        # No alleles have bases at this particular location
        if E >= 100.0:
            has_empty = True
            consensus_seq += 'E'
            continue
        if E >= 100.0 - min_var_freq:
            idx = 4
        else:
            idx = freq.index(max(freq[:4]))
        assert idx < 5
        consensus_seq += "ACGT."[idx]
    consensus_seq = ''.join(consensus_seq)

    # Remove dots (deletions)
    skip_pos = set()
    if has_empty and remove_empty:
        for seq_i in range(len(seqs)):
            seqs[seq_i] = list(seqs[seq_i])
        for i in range(len(consensus_seq)):
            if consensus_seq[i] != 'E':
                continue
            skip_pos.add(i)
            for seq_i in range(len(seqs)):
                if i >= len(seqs[seq_i]):
                    continue
                seqs[seq_i][i] = 'E'
        for seq_i in range(len(seqs)):
            seqs[seq_i] = ''.join(seqs[seq_i])
            seqs[seq_i] = seqs[seq_i].replace('E', '')
        consensus_seq = consensus_seq.replace('E', '')

    # Convert a list form of consensus_freq to a dictionary form
    temp_freq = []
    for j in range(len(consensus_freq)):
        if j in skip_pos:
            continue
        freq_dic = {}
        for k in range(len(consensus_freq[j])):
            freq = consensus_freq[j][k]
            if freq <= 0.0:
                continue
            nt = "ACGT."[k]                    
            freq_dic[nt] = freq
        temp_freq.append(freq_dic)
    consensus_freq = temp_freq

    assert len(consensus_seq) == len(consensus_freq)                
    return consensus_seq, consensus_freq



"""
Left-shift deletions if poissble
"""
def leftshift_deletions(backbone_seq, seq, debug = False):
    if len(seq) != len(backbone_seq):
        return seq
    seq = list(seq)
    seq_len = len(seq)
    bp_i = 0
    # Skip the first deletion
    while bp_i < seq_len:
        if seq[bp_i] in "ACGT":
            break
        bp_i += 1

    while bp_i < seq_len:
        bp = seq[bp_i]
        if bp != '.':
            bp_i += 1
            continue
        bp_j = bp_i + 1
        while bp_j < seq_len:
            bp2 = seq[bp_j]
            if bp2 != '.':
                break
            else:
                bp_j += 1

        if bp_j >= seq_len:
            bp_i = bp_j
            break

        if debug:
            print >> sys.stderr, bp_i, bp_j, backbone_seq[bp_i-10:bp_i], backbone_seq[bp_i:bp_j], backbone_seq[bp_j:bp_j+10]
            print >> sys.stderr, bp_i, bp_j, ''.join(seq[bp_i-10:bp_i]), ''.join(seq[bp_i:bp_j]), ''.join(seq[bp_j:bp_j+10])
        prev_i, prev_j = bp_i, bp_j

        while bp_i > 0 and seq[bp_i-1] in "ACGT" and backbone_seq[bp_j-1] in "ACGT":
            if seq[bp_i-1] != backbone_seq[bp_j-1]:
                break
            seq[bp_j-1] = seq[bp_i-1]
            seq[bp_i-1] = '.'
            bp_i -= 1
            bp_j -= 1
        bp_i = bp_j
        while bp_i < seq_len:
            if seq[bp_i] in "ACGT":
                break
            bp_i += 1

        # DK - debugging purposes
        if debug:
            print prev_i, prev_j, ''.join(seq[prev_i-10:prev_i]), ''.join(seq[prev_i:prev_j]), ''.join(seq[prev_j:prev_j+10])

    return ''.join(seq)


"""
"""
def extract_vars(base_fname,
                 base_dname,
                 locus_list,
                 inter_gap,
                 intra_gap,
                 whole_haplotype,
                 min_var_freq,
                 ext_seq_len,
                 leftshift,
                 partial,
                 verbose):
    base_fullpath_name = base_fname
    if base_dname != "" and not os.path.exists(base_dname):
        os.mkdir(base_dname)
        base_fullpath_name = "%s/%s" % (base_dname, base_fname)

    # Download human genome and HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]

    if not typing_common.check_files(HISAT2_fnames):
        typing_common.download_genome_and_index()
    
    # Corresponding genomic loci found by HISAT2 (reference is GRCh38)
    #   e.g. hisat2 --no-unal --score-min C,0 -x grch38/genome -f hisatgenotype_db/HLA/fasta/A_gen.fasta
    locus_file = open(base_fullpath_name + ".locus", 'w')
    left_ext_seq_dic, right_ext_seq_dic = {}, {}
    genes, gene_strand = {}, {}

    # Clone a git repository, hisatgenotype_db
    if not os.path.exists("hisatgenotype_db"):
        typing_common.clone_hisatgenotype_database()
    fasta_dname = "hisatgenotype_db/%s/fasta" % base_fname.upper()

    # Check HLA genes
    gene_names = []
    if base_fname == "hla":
        fasta_fnames = glob.glob("%s/*_gen.fasta" % fasta_dname)
    else:
        assert base_fname in ["codis", "cyp"]
        fasta_fnames = glob.glob("%s/*.fasta" % fasta_dname)
    for gen_fname in fasta_fnames:
        gene_name = gen_fname.split('/')[-1].split('_')[0]
        if gene_name == "hla":
            continue
        gene_names.append(gene_name)

    if locus_list == []:
        locus_list = gene_names

    cigar_re = re.compile('\d+\w')
    remove_locus_list = []
    for gene in locus_list:
        aligner_cmd = ["hisat2"]
        if base_fname in ["hla", "coids"]:
            aligner_cmd += ["--score-min", "C,0"]
        aligner_cmd += ["--no-unal",
                        "-x", "grch38/genome",
                        "-f", "%s/%s_gen.fasta" % (fasta_dname, gene)]
        align_proc = subprocess.Popen(aligner_cmd,
                                      stdout=subprocess.PIPE,
                                      stderr=open("/dev/null", 'w'))
        allele_id = ""
        best_chr, best_left, best_right, best_AS, best_strand = "", -1, -1, -sys.maxint, ''
        for line in align_proc.stdout:
            if line.startswith('@'):
                continue
            line = line.strip()
            cols = line.split()
            temp_allele_id, flag, chr, left, _, cigar_str = cols[:6]
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

            flag = int(flag)
            strand = '-' if flag & 0x10 else '+'
            AS = ""
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith("AS"):
                    AS = col[5:]
            assert AS != ""
            AS = int(AS)
            if AS > best_AS:
                allele_id = temp_allele_id
                best_chr, best_left, best_right, best_AS, best_strand = chr, left, right, AS, strand

        chr, left, right, strand = best_chr, best_left, best_right, best_strand
        align_proc.communicate()
        if allele_id == "":
            remove_locus_list.append(gene)
            continue
        if base_fname == "hla":
            allele_name = ""
            for line in open("%s/%s_gen.fasta" % (fasta_dname, gene)):
                line = line.strip()
                if not line.startswith('>'):
                    continue
                tmp_allele_id, tmp_allele_name = line[1:].split()[:2]
                if allele_id == tmp_allele_id:
                    allele_name = tmp_allele_name
                    break
        else:
            allele_name = allele_id
        assert allele_name != "" and strand != ''
        genes[gene] = allele_name
        gene_strand[gene] = strand
        print >> sys.stderr, "%s-%s's reference allele is %s on '%s' strand of chromosome %s" % \
            (base_fname.upper(), gene, allele_name, strand, chr)

        assert chr != "" and left >= 0 and right > left
        if ext_seq_len > 0:
            left_ext_seq, right_ext_seq = "", ""
            left1, left2 = max(1, left - ext_seq_len), max(1, left - 1)
            if left2 > 0:
                extract_seq_cmd = ["samtools", "faidx", "genome.fa", "%s:%d-%d" % (chr, left1, left2)]
                extract_seq_proc = subprocess.Popen(extract_seq_cmd,
                                                    stdout=subprocess.PIPE,
                                                    stderr=open("/dev/null", 'w'))
                for line in extract_seq_proc.stdout:
                    if line.startswith('>'):
                        continue
                    line = line.strip()
                    left_ext_seq += line
            extract_seq_cmd = ["samtools", "faidx", "genome.fa", "%s:%d-%d" % (chr, right, right + ext_seq_len - 1)]
            extract_seq_proc = subprocess.Popen(extract_seq_cmd,
                                                stdout=subprocess.PIPE,
                                                stderr=open("/dev/null", 'w'))
            for line in extract_seq_proc.stdout:
                if line.startswith('>'):
                    continue
                line = line.strip()
                right_ext_seq += line

            if strand == '-':
                left_ext_seq, right_ext_seq = typing_common.reverse_complement(right_ext_seq), typing_common.reverse_complement(left_ext_seq)
            left_ext_seq_dic[gene], right_ext_seq_dic[gene] = left_ext_seq, right_ext_seq
            

    # Extract exon information from hla.data
    gene_exons = {}
    if base_fname == "hla":        
        skip = False
        for line in open("hisatgenotype_db/%s/hla.dat" % base_fname.upper()):
            if line.startswith("DE"):
                allele_name = line.split()[1][:-1]
                if allele_name.startswith("HLA-"):
                    allele_name = allele_name[4:]
                gene = allele_name.split('*')[0]
                if line.find("partial") != -1 or \
                        not gene in genes or \
                        allele_name != genes[gene]:
                    skip = True
                    continue
                skip = False
            elif not skip:
                if not line.startswith("FT"):
                    continue
                if line.find("exon") != -1:
                    exon_range = line.split()[2].split("..")
                    exon_left, exon_right = int(exon_range[0]) - 1, int(exon_range[1]) - 1
                    assert exon_left >= 0
                    assert exon_left < exon_right
                    if not gene in gene_exons:
                        gene_exons[gene] = []
                    if gene in left_ext_seq_dic:
                        left_ext_seq_len = len(left_ext_seq_dic[gene])
                    else:
                        left_ext_seq_len = 0
                    gene_exons[gene].append([exon_left + left_ext_seq_len, exon_right + left_ext_seq_len])

    tmp_locus_list = []
    for gene in locus_list:
        if gene in remove_locus_list:
            continue
        if base_fname == "hla" and gene not in gene_exons:
            continue
        tmp_locus_list.append(gene)
    locus_list = tmp_locus_list
    for key in genes.keys():
        if key in locus_list:
            continue
        del genes[key]
        del gene_strand[key]

    # Write the backbone sequences into a fasta file
    backbone_file = open(base_fullpath_name + "_backbone.fa", 'w')        
    # variants w.r.t the backbone sequences into a SNP file
    var_file = open(base_fullpath_name + ".snp", 'w')
    var_index_file = open(base_fullpath_name + ".index.snp", 'w')
    # variant frequence
    var_freq_file = open(base_fullpath_name + ".snp.freq", 'w')
    # haplotypes
    haplotype_file = open(base_fullpath_name + ".haplotype", 'w')
    # pairs of a variant and the corresponding HLA allels into a LINK file    
    link_file = open(base_fullpath_name + ".link", 'w')
    # Write all the sequences with dots removed into a file
    input_file = open(base_fullpath_name + "_sequences.fa", 'w')
    num_vars, num_haplotypes = 0, 0
    full_alleles = {}
    for gene, ref_gene in genes.items():
        strand = gene_strand[gene]
        left_ext_seq, right_ext_seq = "", ""
        if gene in left_ext_seq_dic:
            left_ext_seq, right_ext_seq = left_ext_seq_dic[gene], right_ext_seq_dic[gene]

        def read_MSF_file(fname, left_ext_seq = "", right_ext_seq = ""):
            names = {} # HLA allele names to numeric IDs
            seqs = []  # HLA multiple alignment sequences
            for line in open(fname):
                line = line.strip()
                if not line or \
                        not line[0].isalnum():
                    continue

                if line.startswith("MSF"):
                    continue

                if line.startswith("Name"):
                    try:
                        name = line.split('\t')[0]
                        name = name.split()[1]
                    except ValueError:
                        continue

                    if name in names:
                        print >> sys.stderr, "Warning: %s is found more than once in Names" % (name)
                        continue

                    names[name] = len(names)
                else:
                    if len(seqs) == 0:
                        seqs = [left_ext_seq for i in range(len(names))]
                    try:
                        cols = line.split()
                        name = cols[0]
                        fives = cols[1:]
                        assert len(fives) > 0
                    except ValueError:
                        continue

                    if name not in names:
                        names[name] = len(names)

                    id = names[name]
                    if id >= len(seqs):
                        assert id == len(seqs)
                        seqs.append(left_ext_seq)
                        
                    seqs[id] += ''.join(fives)

                    # Add sub-names of the allele
                    sub_name = ""
                    for group in name.split(':')[:-1]:
                        if sub_name != "":
                            sub_name += ":"
                        sub_name += group
                        if sub_name not in full_alleles:
                            full_alleles[sub_name] = [name]
                        else:
                            full_alleles[sub_name].append(name)

            if len(right_ext_seq) > 0:
                for i_ in range(len(seqs)):
                    seqs[i_] += right_ext_seq

            return names, seqs

        if base_fname == "hla":
            MSA_fname = "hisatgenotype_db/%s/msf/%s_gen.msf" % (base_fname.upper(), gene)
        else:
            MSA_fname = "hisatgenotype_db/%s/msf/%s_gen.msf" % (base_fname.upper(), gene)
            
        if not os.path.exists(MSA_fname):
            print >> sys.stderr, "Warning: %s does not exist" % MSA_fname
            continue

        names, seqs = read_MSF_file(MSA_fname, left_ext_seq, right_ext_seq)

        # Identify a consensus sequence
        assert len(seqs) > 0

        # Check sequences are of equal length
        def find_seq_len(seqs):
            seq_lens = {}
            for s in range(len(seqs)):
                seq_len = len(seqs[s])
                if seq_len not in seq_lens:
                    seq_lens[seq_len] = 1
                else:
                    seq_lens[seq_len] += 1

            max_seq_count = 0
            for tmp_seq_len, tmp_seq_count in seq_lens.items():
                if tmp_seq_count > max_seq_count:
                    seq_len = tmp_seq_len
                    max_seq_count = tmp_seq_count
            return seq_len

        seq_len = find_seq_len(seqs)        
        backbone_name = "%s*BACKBONE" % gene
        backbone_seq, backbone_freq = create_consensus_seq(seqs,
                                                           seq_len,
                                                           min_var_freq,
                                                           not partial) # Remove empty sequences?
        # Allele sequences can shrink, so readjust the sequence length
        if not partial:
            seq_len = find_seq_len(seqs)

        if partial and base_fname == "hla":
            partial_MSA_fname = "hisatgenotype_db/HLA/msf/%s_nuc.msf" % gene
            if not os.path.exists(partial_MSA_fname):
                print >> sys.stderr, "Warning: %s does not exist" % partial_MSA_fname
                continue
            partial_names, partial_seqs = read_MSF_file(partial_MSA_fname)

            # DK - debugging purposes
            # Partial alleles vs. Full alleles
            """
            counts = [0, 0, 0, 0]
            for partial_name in partial_names.keys():
                if partial_name in names:
                    continue
                name_group = partial_name.split(':')
                for group_i in [3, 2, 1, 0]:
                    if group_i == 0:
                        counts[group_i] += 1
                    if group_i > len(name_group):
                        continue
                    sub_name = ':'.join(name_group[:group_i])
                    if sub_name in full_alleles:
                        print partial_name, sub_name, full_alleles[sub_name][:5]
                        counts[group_i] += 1
                        break
            print "DK: counts:", counts
            sys.exit(1)
            """
                
            ref_seq = seqs[names[ref_gene]]
            ref_seq_map = create_map(ref_seq)
            ref_partial_seq = partial_seqs[partial_names[ref_gene]]
            ref_partial_seq_map = create_map(ref_partial_seq)
            exons = gene_exons[gene]
            exon_len = 0
            ref_exons = [] # converted exons to MSF file (e.g. A_gen.msf)
            ref_partial_exons = [] # converted exons to MSF file (e.g. A_nuc.msf)

            complete = True
            for exon in exons:
                left, right = exon
                ref_exons.append([ref_seq_map[left], ref_seq_map[right]])
                next_exon_len = right - left + exon_len
                if next_exon_len >= len(ref_partial_seq_map):
                    print >> sys.stderr, "Warning: partial sequences (%s) seem to be incomplete" % gene
                    complete = False
                    break
                ref_partial_exons.append([ref_partial_seq_map[exon_len], ref_partial_seq_map[next_exon_len]])
                exon_len += (right - left + 1)
                # Make sure two MSF files (e.g. A_gen.msf and A_nuc.msf) share the same MSF lengths in the exonic sequences
                ref_exon_len = ref_exons[-1][1] - ref_exons[-1][0] + 1
                ref_partial_exon_len = ref_partial_exons[-1][1] - ref_partial_exons[-1][0] + 1
                assert ref_exon_len == ref_partial_exon_len

            if complete:
                partial_seq_len = find_seq_len(partial_seqs)
                partial_backbone_seq, partial_backbone_freq = create_consensus_seq(partial_seqs,
                                                                                   partial_seq_len,
                                                                                   min_var_freq,
                                                                                   False) # Remove empty sequences?
                for name, seq_id in partial_names.items():
                    if name in names:
                        continue
                    seq = partial_seqs[seq_id]
                    new_seq = ""
                    right = 0
                    for e in range(len(exons)):
                        ref_exon = ref_exons[e]
                        ref_partial_exon = ref_partial_exons[e]
                        new_seq += backbone_seq[right:ref_exon[0]]
                        exon_seq = seq[ref_partial_exon[0]:ref_partial_exon[1] + 1]
                        nt_exon_seq = exon_seq.replace('.', '')
                        if len(nt_exon_seq) == 0:
                            exon_seq = partial_backbone_seq[ref_partial_exon[0]:ref_partial_exon[1] + 1]
                        new_seq += exon_seq
                        right = ref_exon[1] + 1
                    new_seq += backbone_seq[right:]
                    names[name] = len(seqs)
                    seqs.append(new_seq)

                backbone_seq, backbone_freq = create_consensus_seq(seqs,
                                                                   seq_len,
                                                                   min_var_freq,
                                                                   True) # Remove empty sequences?
                seq_len = find_seq_len(seqs)
                
        if min_var_freq <= 0.0:
            assert '.' not in backbone_seq and 'E' not in backbone_seq
        
        # Reverse complement MSF if this gene is on '-' strand
        if strand == '-':
            # Reverse exons
            ref_seq = seqs[names[ref_gene]]
            ref_seq = ref_seq.replace('.', '')
            ref_seq_len = len(ref_seq)
            if base_fname == "hla":
                exons = []
                for left, right in reversed(gene_exons[gene]):
                    left, right = ref_seq_len - right - 1, ref_seq_len - left - 1
                    exons.append([left, right])
                gene_exons[gene] = exons

            for i in range(len(seqs)):
                seqs[i] = typing_common.reverse_complement(seqs[i])
            backbone_seq, backbone_freq = create_consensus_seq(seqs, seq_len, min_var_freq, True)

        if leftshift:
            for seq_i in range(len(seqs)):
                seqs[seq_i] = leftshift_deletions(backbone_seq, seqs[seq_i])
            backbone_seq, backbone_freq = create_consensus_seq(seqs, seq_len, min_var_freq, True)
            seq_len = find_seq_len(seqs)

        print >> sys.stderr, "%s: number of HLA alleles is %d." % (gene, len(names))

        Vars = {}
        for cmp_name, id in names.items():
            if cmp_name == backbone_name:
                continue
            assert id < len(seqs)
            cmp_seq = seqs[id]
            if len(cmp_seq) != seq_len:
                print >> sys.stderr, "Warning: the length of %s (%d) is different from %d" % \
                    (cmp_name, len(cmp_seq), seq_len)
                continue

            # DK - debugging purposes
            """
            if cmp_name == "A*03:01:07":
                print cmp_name
                cmp_seq2 = seqs[names["A*32:29"]]
                for s in range(0, seq_len, 100):
                    print s, backbone_seq[s:s+100]
                    print s, cmp_seq2[s:s+100]
                    print s, cmp_seq[s:s+100]
                # sys.exit(1)
            """
            def insertVar(type, info):
                pos, backbone_pos, data = info
                if type in "MI":
                    varKey = "%d-%s-%s" % (pos, type, data)
                else:
                    varKey = "%d-%s-%d" % (pos, type, data)

                if varKey not in Vars:
                    if type == 'M':
                        assert backbone_pos < backbone_freq
                        assert data in backbone_freq[backbone_pos]
                        freq = backbone_freq[backbone_pos][data]
                    elif type == 'D':
                        del_len = int(data)
                        freq = 100.0
                        assert backbone_pos + del_len <= backbone_freq
                        for d in range(del_len):
                            assert '.' in backbone_freq[backbone_pos + d]
                            freq2 = backbone_freq[backbone_pos + d]['.']
                            if freq2 < freq:
                                freq = freq2
                    else:
                        assert type == 'I'
                        ins_len = len(data)
                        freq = 100.0
                        assert backbone_pos + ins_len <= backbone_freq
                        for i in range(ins_len):
                            nt = data[i]
                            assert nt in backbone_freq[backbone_pos + i]
                            freq2 = backbone_freq[backbone_pos + i][nt]
                            if freq2 < freq:
                                freq = freq2
                        assert freq <= min_var_freq
                    
                    Vars[varKey] = [freq, [cmp_name]]
                else:
                    Vars[varKey][1].append(cmp_name)

            insertion, deletion = [], []
            ndots = 0
            for s in range(seq_len):
                assert not (insertion and deletion)
                bc = backbone_seq[s]
                cc = cmp_seq[s]
                if bc != '.' and cc != '.':
                    if insertion:
                        insertVar('I', insertion)
                        insertion = []
                    elif deletion:
                        insertVar('D', deletion)
                        deletion = []
                    if bc != cc:
                        mismatch = [s - ndots, s, cc]
                        insertVar('M', mismatch)
                elif bc == '.' and cc != '.':
                    if deletion:
                        insertVar('D', deletion)
                        deletion = []
                    if insertion:
                        insertion[2] += cc
                    else:
                        insertion = [s - ndots, s, cc]
                elif bc != '.' and cc == '.':
                    if insertion:
                        insertVar('I', insertion)
                        insertion = []
                    if deletion:
                        deletion[2] += 1
                    else:
                        deletion = [s - ndots, s, 1]

                if bc == '.':
                    ndots += 1

                """
                if backbone_seq[s] != cmp_seq[s]:
                    print "%s is different %s at %d: %s vs. %s" % \
                        (backbone_name, cmp_name, s+1, backbone_seq[s], cmp_seq[s])
                """

            if insertion:
                insertVar('I', insertion)
            elif deletion:
                insertVar('D', deletion)


        print >> sys.stderr, "Number of variants is %d." % (len(Vars.keys()))

        # Compare variants
        def cmp_varKey(a, b):
            a_locus, a_type, a_data = a.split('-')
            b_locus, b_type, b_data = b.split('-')
            a_locus, b_locus = int(a_locus), int(b_locus)
            if a_locus != b_locus:
                return a_locus - b_locus
            if a_type != b_type:
                if a_type == 'I':
                    return -1
                elif b_type == 'I':
                    return 1
                elif a_type == 'M':
                    return -1
                else:
                    assert b_type == 'M'
                    return 1
            assert a_data != b_data
            if a_type in "MI":
                if a_data < b_data:
                    return -1
                else:
                    return 1
            else:
                assert a_type == 'D'
                return int(a_data) - int(b_data)            

        Vars_ = {}
        for key, values in Vars.items():
            freq, names_ = values
            for name in names_:
                if not name in Vars_:
                    Vars_[name] = [key]
                else:
                    Vars_[name].append(key)
        for name, vars in Vars_.items():
            Vars_[name] = sorted(vars, cmp=cmp_varKey)

        # Sanity check -
        #    (1) Reconstruct the other sequences from the backbone sequence and variants and
        #    (2) Confirm these constructed sequences are the same as those input sequences.
        for cmp_name, id in names.items():
            if cmp_name == backbone_name:
                continue

            constr_seq = backbone_seq.replace('.', '')
            constr_seq = list(constr_seq)
            locus_diff = 0

            if cmp_name not in Vars_:
                continue
            
            for var in Vars_[cmp_name]:
                try:
                    locus, type, data = var.split('-')
                    locus = int(locus)
                except ValueError:
                    continue

                if type == 'M':
                    assert len(data) == 1
                    constr_seq[locus + locus_diff] = data[0]
                elif type == 'I':
                    assert locus + locus_diff >= 0
                    assert locus + locus_diff <= len(constr_seq)
                    constr_seq = constr_seq[:locus + locus_diff] + list(data) + constr_seq[locus + locus_diff:]
                    locus_diff += len(data)
                else:
                    assert type == 'D'
                    assert locus + locus_diff + len(data) <= len(constr_seq)
                    assert locus + locus_diff >= 0
                    del_len = int(data)
                    constr_seq = constr_seq[:locus + locus_diff] + constr_seq[locus + locus_diff + del_len:]
                    locus_diff -= del_len

            constr_seq = "".join(constr_seq)
            assert id < len(seqs)
            cmp_seq = seqs[id].replace('.', '')
            if len(constr_seq) != len(cmp_seq):
                print >> sys.stderr, "Error: reconstruction fails (%s)! Lengths different: %d vs. %d" % \
                    (cmp_name, len(constr_seq), len(cmp_seq))
                assert False

            # Sanity check
            for s in range(len(constr_seq)):
                if constr_seq[s] != cmp_seq[s]:
                    print >> sys.stderr, "Differ at %d: %s vs. %s (reconstruction vs. original)" % \
                        (s, constr_seq[s], cmp_seq[s])
                    print "%s:%s vs. %s:%s" % \
                        (constr_seq[s-10:s], constr_seq[s:s+10], cmp_seq[s-10:s], cmp_seq[s:s+10])

            if constr_seq != cmp_seq.replace('.', ''):
                print >> sys.stderr, "Error: reconstruction fails for %s" % (cmp_name)
                assert False

        # Write the backbone sequences into a fasta file
        print >> backbone_file, ">%s" % (backbone_name)
        backbone_seq_ = backbone_seq.replace('.', '')
        for s in range(0, len(backbone_seq_), 60):
            print >> backbone_file, backbone_seq_[s:s+60]

        # Remap the backbone allele, which is sometimes slighly different from
        #   fasta version
        ref_backbone_id = names[ref_gene]
        ref_backbone_seq = seqs[ref_backbone_id]
        aligner_cmd = ["hisat2"]
        if base_fname == "hla":
            aligner_cmd += ["--score-min", "C,0"]
        aligner_cmd += ["--no-unal",
                        "-x", "grch38/genome",
                        "-f", 
                        "-c", "%s" % ref_backbone_seq.replace('.', '')]
        align_proc = subprocess.Popen(aligner_cmd,
                                      stdout=subprocess.PIPE,
                                      stderr=open("/dev/null", 'w'))
        best_chr, best_left, best_right, best_AS = "", 0, 0, -sys.maxint
        for line in align_proc.stdout:
            if line.startswith('@'):
                continue
            line = line.strip()
            cols = line.split()
            allele_id, flag, chr, left, mapQ, cigar_str = cols[:6]
            flag = int(flag)
            assert flag & 0x10 == 0
            left = int(left) - 1
            right = left
            AS = ""
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith("AS"):
                    AS = col[5:]
            AS = int(AS)
            cigars = cigar_re.findall(cigar_str)
            cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
            for i in range(len(cigars)):
                cigar_op, length = cigars[i]
                if cigar_op in "MND":
                    right += length
            if AS > best_AS:
                best_chr, best_left, best_right, best_AS = chr, left, right, AS

        chr, left, right = best_chr, best_left, best_right
        align_proc.communicate()
        if left == right:
            print >> sys.stderr, "Warning: %s (%s) is not remapped" % (gene, ref_gene)
            continue
        assert left < right

        base_locus = 0                
        ref_seq = seqs[names[ref_gene]]
        ref_seq_map = create_map(ref_seq)

        del_count = []
        for nt in backbone_seq:
            assert nt in "ACGT."
            add = 1 if nt == '.' else 0
            if len(del_count) == 0:
                del_count.append(add)
            else:
                del_count.append(del_count[-1] + add)
        
        if base_fname == "hla":
            exon_str = ""
            for exon_left, exon_right in gene_exons[gene]:
                exon_left, exon_right = ref_seq_map[exon_left], ref_seq_map[exon_right]
                exon_left -= del_count[exon_left]
                exon_right -= del_count[exon_right]
                if exon_str != "":
                    exon_str += ','
                exon_str += ("%d-%d" % (exon_left, exon_right))

            # Sanity check for exonic sequence
            sanity_check = True
            if sanity_check and \
               os.path.exists("hisatgenotype_db/HLA/fasta/%s_nuc.fasta" % gene):
                exons_ = []
                for exon in exon_str.split(','):
                    exon_left, exon_right = exon.split('-')
                    exon_left, exon_right = int(exon_left), int(exon_right)
                    exons_.append([exon_left, exon_right])

                backbone_seq_ = backbone_seq.replace('.', '')
                if ref_gene in Vars_:
                    vars_ = Vars_[ref_gene]
                else:
                    vars_ = []
                seq_ = list(backbone_seq_)
                has_insertion = False
                for var_ in vars_:
                    var_pos, var_type, var_data = var_.split('-')
                    var_pos = int(var_pos)
                    assert var_pos >= 0 and var_pos < len(backbone_seq_)
                    if var_type == 'M':
                        seq_[var_pos] = var_data
                    elif var_type == 'D':
                        del_len = int(var_data)
                        assert var_pos + del_len <= len(ref_seq)
                        seq_[var_pos:var_pos + del_len] = ['.'] * del_len
                    else:
                        assert var_type == 'I'
                        has_insertion = True

                seq_ = ''.join(seq_)
                exon_seq_ = ""
                for exon_left, exon_right in exons_:
                    exon_seq_ += seq_[exon_left:exon_right+1]
                exon_seq_ = exon_seq_.replace('.', '')
                if gene_strand[gene] == '-':
                    exon_seq_ = typing_common.reverse_complement(exon_seq_)

                cmp_exon_seq_, allele_name_ = "", ""
                for line in open("hisatgenotype_db/HLA/fasta/%s_nuc.fasta" % gene):
                    if line.startswith(">"):
                        if allele_name_ == ref_gene:
                            break
                        allele_name_ = line.strip().split()[1]
                        cmp_exon_seq_ = ""
                    else:
                        cmp_exon_seq_ += line.strip()
                """
                print "Has insertions:", has_insertion
                print "constructed:", len(exon_seq_)
                for p in range(0, len(exon_seq_), 60):
                    print exon_seq_[p:p+60]
                print "true:", len(cmp_exon_seq_)
                for p in range(0, len(cmp_exon_seq_), 60):
                    print cmp_exon_seq_[p:p+60]
                """
                if exon_seq_ != cmp_exon_seq_:
                    print >> sys.stderr, "Waring: exonic sequences do not match (%s)" % gene
        else:
            exon_str = "%d-%d" % (left, right - 1)

        print >> locus_file, "%s\t%s\t%d\t%d\t%d\t%s\t%s" % \
            (backbone_name, chr, left, right - 1, len(backbone_seq.replace('.', '')), exon_str, gene_strand[gene])

        # Write
        #       (1) variants w.r.t the backbone sequences into a SNP file
        #       (2) pairs of a variant and the corresponding HLA allels into a LINK file    
        keys = sorted(Vars.keys(), cmp=cmp_varKey)
        var2ID = {}
        for k in range(len(keys)):
            locus, type, data = keys[k].split('-')
            locus = int(locus)
            if type == 'M':
                type_str = "single"
            elif type == 'I':
                type_str = "insertion"
            else:
                assert type == 'D'
                type_str = "deletion"

            freq, names_ = Vars[keys[k]]
            names_ = sorted(names_)            
            varID = "hv%d" % (num_vars)
            tmp_backbone_name = backbone_name
            print >> var_file, "%s\t%s\t%s\t%d\t%s" % \
                (varID, type_str, tmp_backbone_name, base_locus + locus, data)
            if freq >= min_var_freq:
                print >> var_index_file, "%s\t%s\t%s\t%d\t%s" % \
                    (varID, type_str, tmp_backbone_name, base_locus + locus, data)
            print >> var_freq_file, "%s\t%.2f" % (varID, freq)
            print >> link_file, "%s\t%s" % (varID, ' '.join(names_))
            var2ID[keys[k]] = num_vars
            num_vars += 1

        add_seq_len = 0
        # Write haplotypes
        excluded_vars = set()
        var_leftmost, var_rightmost = sys.maxint, -1
        for k in range(len(keys)):
            key = keys[k]
            if Vars[key][0] < min_var_freq:
                excluded_vars.add(key)

            # Update leftmost and rightmost of Vars
            locus, type, data = key.split('-')
            left = right = int(locus)
            if type == 'D':
                right = left + int(data) - 1
            if k == 0:
                var_leftmost = left
            if var_rightmost < right:
                var_rightmost = right

        i = 0
        while i < len(keys):
            key_i = keys[i]
            locus, type, data = key_i.split('-')
            locus = int(locus)
            if type == 'D':
                locus += (int(data) - 1)
            prev_locus = locus
            if whole_haplotype:
                j = len(keys)
            else:
                j = i + 1
                while j < len(keys):
                    key_j = keys[j]
                    locus2, type2, data2 = key_j.split('-')
                    locus2 = int(locus2)
                    if prev_locus + inter_gap < locus2:
                        break
                    prev_locus = locus2
                    if type == 'D':
                        prev_locus += (int(data) - 1)
                    j += 1

            alleles = set()
            for k in range(i, j):
                key_k = keys[k]
                freq, names_ = Vars[key_k]
                if freq < min_var_freq:
                    continue
                add_alleles = set(names_)
                alleles |= add_alleles

            haplotypes = set()
            cur_vars = set(keys[i:j]) - excluded_vars
            for allele in alleles:
                allele_vars = set(Vars_[allele]) - excluded_vars
                allele_cur_vars = '#'.join(sorted(list(cur_vars & allele_vars), cmp=cmp_varKey))
                haplotypes.add(allele_cur_vars)

            # Split some haplotypes that include large gaps inside
            def split_haplotypes(haplotypes):
                split_haplotypes = set()
                for haplotype in haplotypes:
                    haplotype = haplotype.split('#')
                    assert len(haplotype) > 0
                    if len(haplotype) == 1:
                        split_haplotypes.add(haplotype[0])
                        continue
                    prev_s, s = 0, 1
                    while s < len(haplotype):
                        prev_locus, prev_type, prev_data = haplotype[s-1].split('-')
                        locus, type, data = haplotype[s].split('-')
                        prev_locus, locus = int(prev_locus), int(locus)
                        if prev_type == 'D':
                            prev_locus += (int(prev_data) - 1)
                        if prev_locus + intra_gap < locus:
                            split_haplotypes.add('#'.join(haplotype[prev_s:s]))
                            prev_s = s
                        s += 1
                        if s == len(haplotype):
                            split_haplotypes.add('#'.join(haplotype[prev_s:s]))
                return split_haplotypes

            if not whole_haplotype:
                haplotypes = split_haplotypes(haplotypes)

            def cmp_haplotype(a, b):
                a = a.split('#')
                a1_locus, _, _ = a[0].split('-')
                a2_locus, a2_type, a2_data = a[-1].split('-')
                a_begin, a_end = int(a1_locus), int(a2_locus)
                if a2_type == 'D':
                    a_end += (int(a2_data) - 1)
                b = b.split('#')
                b1_locus, _, _ = b[0].split('-')
                b2_locus, b2_type, b2_data = b[-1].split('-')
                b_begin, b_end = int(b1_locus), int(b2_locus)
                if b2_type == 'D':
                    b_end += (int(b2_data) - 1)
                if a_begin != b_begin:
                    return a_begin - b_begin
                return a_end - b_end

            haplotypes = sorted(list(haplotypes), cmp=cmp_haplotype)
            
            # DK - for debugging purposes
            """
            dis = prev_locus - locus
            print "\n[%d, %d]: %d haplotypes" % (i, j, len(haplotypes)), dis
            if len(cur_vars) in range(0, 1000):
                # print "vars:", sorted(list(cur_vars), cmp=cmp_varKey
                print "num:", len(haplotypes)
                for haplotype in haplotypes:
                    print haplotype.split('#')
                print "\nnum:", len(haplotypes2)
                for haplotype in haplotypes2:
                    print haplotype.split('#')
            """

            # Write haplotypes
            sanity_vars = set()
            for h_i in range(len(haplotypes)):
                h = haplotypes[h_i].split('#')
                varIDs = []
                for var in h:
                    varIDs.append("hv%s" % var2ID[var])
                    # DK - for debugging purposes
                    # varIDs.append(var)
                    sanity_vars.add(var2ID[var])
                if whole_haplotype:
                    h_begin, h_end = var_leftmost, var_rightmost
                else:
                    h1_locus, _, _ = h[0].split('-')
                    h2_locus, h2_type, h2_data = h[-1].split('-')
                    h_begin, h_end = int(h1_locus), int(h2_locus)
                    if h2_type == 'D':
                        h_end += (int(h2_data) - 1)
                    assert h_begin <= h_end
                    h_new_begin = h_begin
                    for h_j in reversed(range(0, h_i)):
                        hc = haplotypes[h_j].split('#')
                        hc_begin, hc_type, hc_data = hc[-1].split('-')
                        hc_begin = int(hc_begin)
                        hc_end = hc_begin
                        if hc_type == 'D':
                            hc_end += (int(hc_data) - 1)
                        if hc_end + inter_gap < h_begin:
                            break
                        if h_new_begin > hc_end:
                            h_new_begin = hc_end
                    assert h_new_begin <= h_begin
                    h_begin = h_new_begin
                tmp_backbone_name = backbone_name
                print >> haplotype_file, "ht%d\t%s\t%d\t%d\t%s" % \
                    (num_haplotypes, tmp_backbone_name, base_locus + h_begin, base_locus + h_end, ','.join(varIDs))
                num_haplotypes += 1
                add_seq_len += (h_end - h_begin + 1)
            assert len(sanity_vars) == len(cur_vars)
                    
            i = j

        print >> sys.stderr, "Length of additional sequences for haplotypes:", add_seq_len
                    
        # Write all the sequences with dots removed into a file
        for name, ID in names.items():
            print >> input_file, ">%s" % (name)
            assert ID < len(seqs)
            seq = seqs[ID].replace('.', '')
            for s in range(0, len(seq), 60):
                print >> input_file, seq[s:s+60]

    backbone_file.close()
    locus_file.close()
    var_file.close()
    var_index_file.close()
    var_freq_file.close()
    haplotype_file.close()
    link_file.close()
    input_file.close()

    # Read partial alleles from hla.data, and write them into a file
    partial_allele_list = []
    if base_fname == "hla":
        for line in open("hisatgenotype_db/HLA/hla.dat"):
            if not line.startswith("DE"):
                continue
            allele_name = line.split()[1][:-1]
            if allele_name.startswith("HLA-"):
                allele_name = allele_name[4:]
            gene = allele_name.split('*')[0]
            if line.find("partial") != -1:
                partial_allele_list.append(allele_name)

    partial_file = open("%s.partial" % base_fullpath_name, 'w')
    for partial_allele in partial_allele_list:
        print >> partial_file, partial_allele
    partial_file.close()
   
    
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Extract variants from multiple sequence alignments")
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="hla",
                        help="base filename for backbone sequence, variants, and linking info (Default: hla)")
    parser.add_argument("--locus-list",
                        dest="locus_list",
                        type=str,
                        default="",
                        help="A comma-separated list of gene names (default: empty, all genes)")
    parser.add_argument("--inter-gap",
                        dest="inter_gap",
                        type=int,
                        default=30,
                        help="Maximum distance for variants to be in the same haplotype (default: 30)")
    parser.add_argument("--intra-gap",
                        dest="intra_gap",
                        type=int,
                        default=50,
                        help="Break a haplotype into several haplotypes (default: 50)")
    parser.add_argument("--whole-haplotype",
                        dest="whole_haplotype",
                        action="store_true",
                        help="Include partial alleles (e.g. A_nuc.fasta)")
    parser.add_argument("--min-var-freq",
                        dest="min_var_freq",
                        type=float,
                        default=0.0,
                        help="Exclude variants whose freq is below than this value in percentage (Default: 0.0)")    
    parser.add_argument("--ext-seq",
                        dest="ext_seq_len",
                        type=int,
                        default=0,
                        help="Length of extra sequences flanking backbone sequences (Default: 0)")
    parser.add_argument("--leftshift",
                        dest="leftshift",
                        action="store_true",
                        help="Shift deletions to the leftmost")
    parser.add_argument("--no-partial",
                        dest="partial",
                        action="store_false",
                        help="Exclude partial alleles, exon-only sequences in HLA")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.locus_list == "":
        locus_list = []
    else:
        locus_list = args.locus_list.split(',')
    if args.inter_gap > args.intra_gap:
        print >> sys.stderr, "Error: --inter-gap (%d) must be smaller than --intra-gap (%d)" % (args.inter_gap, args.intra_gap)
        sys.exit(1)
             
    if args.base_fname.find('/') != -1:
        elems = args.base_fname.split('/')
        base_fname = elems[-1]
        base_dname = '/'.join(elems[:-1])
    else:
        base_fname = args.base_fname
        base_dname = ""
        
    extract_vars(base_fname,
                 base_dname,
                 locus_list,
                 args.inter_gap,
                 args.intra_gap,
                 args.whole_haplotype,
                 args.min_var_freq,
                 args.ext_seq_len,
                 args.leftshift,
                 args.partial,
                 args.verbose)

