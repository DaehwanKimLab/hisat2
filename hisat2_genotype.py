#!/usr/bin/env python

#
# Copyright 2016, Daehwan Kim <infphilo@gmail.com>
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


import sys, os, subprocess, re
import inspect, random
import math
from argparse import ArgumentParser, FileType


"""
"""
def read_genome(genome_file):
    chr_dic, chr_names = {}, []
    chr_name, sequence = "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
                chr_names.append(chr_name)
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
        chr_names.append(chr_name)
    return chr_dic, chr_names


"""
"""
def genotype(reference_type,
             base_fname,
             fastq,
             read_fnames,
             threads,
             simulate_interval,
             num_mismatch,
             verbose,
             daehwan_debug):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(genotype))
    ex_path = os.path.dirname(curr_script)

    # Load genomic sequences
    chr_dic, chr_names = read_genome(open("%s.fa" % base_fname))

    # variants, backbone sequence, and other sequeces
    genotype_fnames = ["%s.fa" % base_fname,
                       "%s.gene" % base_fname,
                       "%s.snp" % base_fname,
                       "%s.haplotype" % base_fname,
                       "%s.link" % base_fname,
                       "%s.coord" % base_fname,
                       "%s.clnsig" % base_fname]
    # hisat2 graph index files
    genotype_fnames += ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    if not check_files(genotype_fnames):
        print >> sys.stderr, "Error: some files are missing!"
        sys.exit(1)

    # Align reads, and sort the alignments into a BAM file
    hisat2 = os.path.join(ex_path, "hisat2")
    aligner_cmd = [hisat2,
                   "--no-unal",
                   "-p", str(threads)]
    # aligner_cmd += ["--mm"]
    aligner_cmd += ["-x", "%s" % base_fname]

    assert len(read_fnames) > 0
    if not fastq:
        aligner_cmd += ["-f"]
    single = len(read_fnames) == 1
    if single:
        aligner_cmd += [read_fnames[0]]
    else:
        aligner_cmd += ["-1", read_fnames[0],
                        "-2", read_fnames[1]]

    if verbose:
        print >> sys.stderr, ' '.join(aligner_cmd)

    align_proc = subprocess.Popen(aligner_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=open("/dev/null", 'w'))

    sambam_cmd = ["samtools",
                  "view",
                  "-bS",
                  "-"]
    sambam_proc = subprocess.Popen(sambam_cmd,
                                   stdin=align_proc.stdout,
                                   stdout=open("hla_input_unsorted.bam", 'w'),
                                   stderr=open("/dev/null", 'w'))
    sambam_proc.communicate()
    bamsort_cmd = ["samtools",
                   "sort",
                   "--threads", str(threads),
                   "hla_input_unsorted.bam"]
    bamsort_proc = subprocess.Popen(bamsort_cmd,
                                    stdout=open("hla_input.bam", 'w'),
                                    stderr=open("/dev/null", 'w'))
    bamsort_proc.communicate()

    bamindex_cmd = ["samtools",
                    "index",
                    "hla_input.bam"]
    bamindex_proc = subprocess.Popen(bamindex_cmd,
                                     stderr=open("/dev/null", 'w'))
    bamindex_proc.communicate()

    os.system("rm hla_input_unsorted.bam")

    # Read partial alleles from hla.data (temporary)
    partial_alleles = set()
    """
    for line in open("IMGTHLA/hla.dat"):
        if not line.startswith("DE"):
            continue
        allele_name = line.split()[1][4:-1]
        gene = allele_name.split('*')[0]
        if line.find("partial") != -1:
            partial_alleles.add(allele_name)
    """

    # Read HLA alleles (names and sequences)
    genes, gene_loci, gene_seqs = {}, {}, {}
    for line in open("%s.gene" % base_fname):
        family, allele_name, chr, left, right = line.strip().split()
        gene_name = "%s-%s" % (family, allele_name.split('*')[0])
        assert gene_name not in genes
        genes[gene_name] = allele_name
        left, right = int(left), int(right)
        """
        exons = []
        for exon in exon_str.split(','):
            exon_left, exon_right = exon.split('-')
            exons.append([int(exon_left), int(exon_right)])
        """
        gene_loci[gene_name] = [allele_name, chr, left, right]
        assert chr in chr_dic
        chr_seq = chr_dic[chr]
        assert left < right
        assert right < len(chr_seq)
        gene_seqs[gene_name] = chr_dic[chr][left:right+1]

    # Read link information
    Links, var_genes, allele_vars = {}, {}, {}
    for line in open("%s.link" % base_fname):
        var_id, alleles = line.strip().split('\t')
        alleles = alleles.split()
        assert not var_id in Links
        Links[var_id] = alleles
        for allele in alleles:
            if allele not in allele_vars:
                allele_vars[allele] = set()
            allele_vars[allele].add(var_id)
            gene_name = "HLA-%s" % (allele.split('*')[0])
            var_genes[var_id] = gene_name

    # gene alleles
    allele_names = {}
    for gene_name in genes.keys():
        if gene_name not in allele_names:
            allele_names[gene_name] = []
        gene_name2 = gene_name.split('-')[1]
        for allele_name in allele_vars.keys():
            allele_name1 = allele_name.split('*')[0]
            if gene_name2 == allele_name1:
                allele_names[gene_name].append(allele_name)


    # Read HLA variants, and link information
    Vars, Var_list = {}, {}
    for line in open("%s.snp" % base_fname):
        var_id, var_type, chr, pos, data = line.strip().split('\t')
        pos = int(pos)

        # daehwan - for debugging purposes
        if var_id not in var_genes:
            continue
        
        assert var_id in var_genes
        gene_name = var_genes[var_id]
        """
        if reference_type != "gene":
            allele, dist = None, 0
            for tmp_gene, values in refHLA_loci.items():
                allele_name, chr, left, right, exons = values
                if allele == None:
                    allele = allele_name
                    dist = abs(pos - left)
                else:
                    if dist > abs(pos - left):
                        allele = allele_name
                        dist = abs(pos - left)
        """
            
        if not gene_name in Vars:
            Vars[gene_name] = {}
            assert not gene_name in Var_list
            Var_list[gene_name] = []
            
        assert not var_id in Vars[gene_name]
        """
        left = 0
        if reference_type != "gene":
            _, _, left, _, _ = refHLA_loci[gene]
        """
        Vars[gene_name][var_id] = [var_type, pos, data]
        Var_list[gene_name].append([pos, var_id])

    for gene_name, in_var_list in Var_list.items():
        Var_list[gene_name] = sorted(in_var_list)
    def lower_bound(Var_list, pos):
        low, high = 0, len(Var_list)
        while low < high:
            m = (low + high) / 2
            m_pos = Var_list[m][0]
            if m_pos < pos:
                low = m + 1
            elif m_pos > pos:
                high = m
            else:
                assert m_pos == pos
                while m > 0:
                    if Var_list[m-1][0] < pos:
                        break
                    m -= 1
                return m
        return low       
           
    
    # HLA gene allele lengths
    """
    HLA_lengths = {}
    for HLA_gene, HLA_alleles in HLAs.items():
        HLA_lengths[HLA_gene] = {}
        for allele_name, seq in HLA_alleles.items():
            HLA_lengths[HLA_gene][allele_name] = len(seq)
    """

    # Cigar regular expression
    cigar_re = re.compile('\d+\w')

    test_list = [[sorted(genes.keys())]]
    for test_i in range(len(test_list)):
        test_HLA_list = test_list[test_i]
        for test_HLA_names in test_HLA_list:
            print >> sys.stderr, "\t%s" % (test_HLA_names)
            for gene in test_HLA_names:
                ref_allele = genes[gene]
                ref_seq = gene_seqs[gene]
                # ref_exons = refHLA_loci[gene][-1]

                # Read alignments
                alignview_cmd = ["samtools",
                                 "view"]
                alignview_cmd += ["hla_input.bam"]
                base_locus = 0
                _, chr, left, right = gene_loci[gene]
                base_locus = left
                alignview_cmd += ["%s:%d-%d" % (chr, left + 1, right + 1)]

                bamview_proc = subprocess.Popen(alignview_cmd,
                                                stdout=subprocess.PIPE,
                                                stderr=open("/dev/null", 'w'))

                sort_read_cmd = ["sort", "-k", "1", "-n"]
                alignview_proc = subprocess.Popen(sort_read_cmd,
                                                  stdin=bamview_proc.stdout,
                                                  stdout=subprocess.PIPE,
                                                  stderr=open("/dev/null", 'w'))

                # Count alleles
                HLA_counts, HLA_cmpt = {}, {}
                coverage = [0 for i in range(len(ref_seq) + 1)]
                num_reads, total_read_len = 0, 0
                prev_read_id = None
                prev_exon = False
                for line in alignview_proc.stdout:
                    cols = line.strip().split()
                    read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
                    read_seq, qual = cols[9], cols[10]
                    num_reads += 1
                    total_read_len += len(read_seq)
                    flag, pos = int(flag), int(pos)
                    pos -= 1
                    if pos < 0:
                        continue

                    if flag & 0x4 != 0:
                        continue

                    NM, Zs, MD = "", "", ""
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("Zs"):
                            Zs = col[5:]
                        elif col.startswith("MD"):
                            MD = col[5:]
                        elif col.startswith("NM"):
                            NM = int(col[5:])

                    if NM > num_mismatch:
                        continue

                    # daehwan - for debugging purposes
                    debug = False
                    if read_id in ["2339"] and False:
                        debug = True
                        print "read_id: %s)" % read_id, pos, cigar_str, "NM:", NM, MD, Zs
                        print "            ", read_seq

                    vars = []
                    if Zs:
                        vars = Zs.split(',')

                    assert MD != ""
                    MD_str_pos, MD_len = 0, 0
                    read_pos, left_pos = 0, pos
                    right_pos = left_pos
                    cigars = cigar_re.findall(cigar_str)
                    cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
                    cmp_list = []
                    for i in range(len(cigars)):
                        cigar_op, length = cigars[i]
                        if cigar_op == 'M':
                            first = True
                            MD_len_used = 0
                            while True:
                                if not first or MD_len == 0:
                                    if MD[MD_str_pos].isdigit():
                                        num = int(MD[MD_str_pos])
                                        MD_str_pos += 1
                                        while MD_str_pos < len(MD):
                                            if MD[MD_str_pos].isdigit():
                                                num = num * 10 + int(MD[MD_str_pos])
                                                MD_str_pos += 1
                                            else:
                                                break
                                        MD_len += num
                                # Insertion or full match followed
                                if MD_len >= length:
                                    MD_len -= length
                                    cmp_list.append(["match", right_pos + MD_len_used, length - MD_len_used])
                                    break
                                first = False
                                read_base = read_seq[read_pos + MD_len]
                                MD_ref_base = MD[MD_str_pos]
                                MD_str_pos += 1
                                assert MD_ref_base in "ACGT"
                                cmp_list.append(["match", right_pos + MD_len_used, MD_len - MD_len_used])
                                cmp_list.append(["mismatch", right_pos + MD_len, 1])
                                MD_len_used = MD_len + 1
                                MD_len += 1
                                # Full match
                                if MD_len == length:
                                    MD_len = 0
                                    break
                        elif cigar_op == 'I':
                            cmp_list.append(["insertion", right_pos, length])
                        elif cigar_op == 'D':
                            if MD[MD_str_pos] == '0':
                                MD_str_pos += 1
                            assert MD[MD_str_pos] == '^'
                            MD_str_pos += 1
                            while MD_str_pos < len(MD):
                                if not MD[MD_str_pos] in "ACGT":
                                    break
                                MD_str_pos += 1
                            cmp_list.append(["deletion", right_pos, length])
                        elif cigar_op == 'S':
                            cmp_list.append(["soft", right_pos, length])
                        else:                    
                            assert cigar_op == 'N'
                            cmp_list.append(["intron", right_pos, length])

                        if cigar_op in "MND":
                            right_pos += length

                        if cigar_op in "MIS":
                            read_pos += length

                    """
                    exon = False
                    for exon in ref_exons:
                        exon_left, exon_right = exon
                        if right_pos <= exon_left or pos > exon_right:
                            continue
                        else:
                            exon = True
                            break
                    """

                    if left_pos < base_locus or \
                            right_pos - base_locus > len(ref_seq):
                        continue
                
                    def add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, exon = True):
                        max_count = max(HLA_count_per_read.values())
                        cur_cmpt = set()
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue
                            """
                            if allele in exclude_allele_list:
                                continue
                            """
                            cur_cmpt.add(allele)                    
                            if not allele in HLA_counts:
                                HLA_counts[allele] = 1
                            else:
                                HLA_counts[allele] += 1

                        if len(cur_cmpt) == 0:
                            return

                        # daehwan - for debugging purposes                            
                        alleles = ["", ""]
                        # alleles = ["B*40:304", "B*40:02:01"]
                        allele1_found, allele2_found = False, False
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue
                            if allele == alleles[0]:
                                allele1_found = True
                            elif allele == alleles[1]:
                                allele2_found = True
                        if allele1_found != allele2_found:
                            print alleles[0], HLA_count_per_read[alleles[0]]
                            print alleles[1], HLA_count_per_read[alleles[1]]
                            if allele1_found:
                                print ("%s\tread_id %s - %d vs. %d]" % (alleles[0], prev_read_id, max_count, HLA_count_per_read[alleles[1]]))
                            else:
                                print ("%s\tread_id %s - %d vs. %d]" % (alleles[1], prev_read_id, max_count, HLA_count_per_read[alleles[0]]))
                            print read_seq

                        cur_cmpt = sorted(list(cur_cmpt))
                        cur_cmpt = '-'.join(cur_cmpt)
                        add = 1
                        """
                        if partial and not exon:
                            add *= 0.2
                        """
                        if not cur_cmpt in HLA_cmpt:
                            HLA_cmpt[cur_cmpt] = add
                        else:
                            HLA_cmpt[cur_cmpt] += add

                    if read_id != prev_read_id:
                        if prev_read_id != None:
                            add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, prev_exon)

                        HLA_count_per_read = {}
                        for HLA_name in allele_names[gene]:
                            if HLA_name.find("BACKBONE") != -1:
                                continue
                            HLA_count_per_read[HLA_name] = 0

                    def add_count(var_id, add):
                        assert var_id in Links
                        alleles = Links[var_id]
                        for allele in alleles:
                            if allele.find("BACKBONE") != -1:
                                continue
                            HLA_count_per_read[allele] += add
                            # daehwan - for debugging purposes
                            if debug:
                                if allele in ["DQA1*05:05:01:01", "DQA1*05:05:01:02"]:
                                    print allele, add, var_id

                    # Decide which allele(s) a read most likely came from
                    # also sanity check - read length, cigar string, and MD string
                    for var_id, data in Vars[gene].items():
                        var_type, var_pos, var_data = data
                        if var_type != "deletion":
                            continue
                        if left_pos >= var_pos and right_pos <= var_pos + int(var_data):
                            add_count(var_id, -1)                            
                    ref_pos, read_pos, cmp_cigar_str, cmp_MD = left_pos, 0, "", ""
                    cigar_match_len, MD_match_len = 0, 0            
                    for cmp in cmp_list:
                        type = cmp[0]
                        length = cmp[2]
                        if type == "match":
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos + length <= var_pos:
                                    break
                                if ref_pos <= var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":
                                        if ref_pos < var_pos and ref_pos + length > var_pos + len(var_data):
                                            add_count(var_id, -1)
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print cmp, var_id, Links[var_id]
                                    elif var_type == "deletion":
                                        del_len = int(var_data)
                                        if ref_pos < var_pos and ref_pos + length > var_pos + del_len:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print cmp, var_id, Links[var_id], -1, Vars[gene][var_id]
                                            # Check if this might be one of the two tandem repeats (the same left coordinate)
                                            cmp_left, cmp_right = cmp[1], cmp[1] + cmp[2]
                                            test1_seq1 = ref_seq[cmp_left-base_locus:cmp_right-base_locus]
                                            test1_seq2 = ref_seq[cmp_left-base_locus:var_pos-base_locus] + ref_seq[var_pos + del_len - base_locus:cmp_right + del_len - base_locus]
                                            # Check if this happens due to small repeats (the same right coordinate - e.g. 19 times of TTTC in DQA1*05:05:01:02)
                                            cmp_left -= read_pos
                                            cmp_right += (len(read_seq) - read_pos - cmp[2])
                                            test2_seq1 = ref_seq[cmp_left+int(var_data)-base_locus:cmp_right-base_locus]
                                            test2_seq2 = ref_seq[cmp_left-base_locus:var_pos-base_locus] + ref_seq[var_pos+int(var_data)-base_locus:cmp_right-base_locus]
                                            if test1_seq1 != test1_seq2 and test2_seq1 != test2_seq2:
                                                add_count(var_id, -1)
                                    else:
                                        if debug:
                                            print cmp, var_id, Links[var_id], -1
                                        add_count(var_id, -1)
                                var_idx += 1

                            read_pos += length
                            ref_pos += length
                            cigar_match_len += length
                            MD_match_len += length
                        elif type == "mismatch":
                            read_base = read_seq[read_pos]
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "single":
                                        if var_data == read_base:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print cmp, var_id, 1, var_data, read_base, Links[var_id]

                                            # daehwan - for debugging purposes
                                            if False:
                                                read_qual = ord(qual[read_pos])
                                                add_count(var_id, (read_qual - 60) / 60.0)
                                            else:
                                                add_count(var_id, 1)
                                        # daehwan - check out if this routine is appropriate
                                        # else:
                                        #    add_count(var_id, -1)
                                var_idx += 1
                            cmp_MD += ("%d%s" % (MD_match_len, ref_seq[ref_pos-base_locus]))
                            MD_match_len = 0
                            cigar_match_len += 1
                            read_pos += 1
                            ref_pos += 1
                        elif type == "insertion":
                            ins_seq = read_seq[read_pos:read_pos+length]
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            # daehwan - for debugging purposes
                            if debug:
                                print left_pos, cigar_str, MD, vars
                                print ref_pos, ins_seq, Var_list[gene][var_idx], Vars[gene][Var_list[gene][var_idx][1]]
                                # sys.exit(1)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":                                
                                        if var_data == ins_seq:
                                            # daehwan - for debugging purposes
                                            if debug:
                                                print cmp, var_id, 1, Links[var_id]
                                            add_count(var_id, 1)
                                var_idx += 1

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dI" % length)
                        elif type == "deletion":
                            del_len = length
                            # Deletions can be shifted bidirectionally
                            temp_ref_pos = ref_pos
                            while temp_ref_pos > 0:
                                last_bp = ref_seq[temp_ref_pos + del_len - 1 - base_locus]
                                prev_bp = ref_seq[temp_ref_pos - 1 - base_locus]
                                if last_bp != prev_bp:
                                    break
                                temp_ref_pos -= 1
                            var_idx = lower_bound(Var_list[gene], temp_ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if temp_ref_pos < var_pos:
                                    first_bp = ref_seq[temp_ref_pos - base_locus]
                                    next_bp = ref_seq[temp_ref_pos + del_len - base_locus]
                                    if first_bp == next_bp:
                                        temp_ref_pos += 1
                                        continue
                                    else:
                                        break
                                if temp_ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "deletion":
                                        var_len = int(var_data)
                                        if var_len == length:
                                            if debug:
                                                print cmp, var_id, 1, Links[var_id]
                                                print ref_seq[var_pos - 10-base_locus:var_pos-base_locus], ref_seq[var_pos-base_locus:var_pos+int(var_data)-base_locus], ref_seq[var_pos+int(var_data)-base_locus:var_pos+int(var_data)+10-base_locus]
                                            add_count(var_id, 1)
                                var_idx += 1

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_MD += ("%d" % MD_match_len)
                            MD_match_len = 0
                            cmp_cigar_str += ("%dD" % length)
                            cmp_MD += ("^%s" % ref_seq[ref_pos-base_locus:ref_pos+length-base_locus])
                            ref_pos += length
                        elif type == "soft":
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dS" % length)
                        else:
                            assert type == "intron"
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_cigar_str += ("%dN" % length)
                            ref_pos += length                    
                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                    cmp_MD += ("%d" % MD_match_len)
                    if read_pos != len(read_seq) or \
                            cmp_cigar_str != cigar_str or \
                            cmp_MD != MD:
                        print >> sys.stderr, "Error:", cigar_str, MD
                        print >> sys.stderr, "\tcomputed:", cmp_cigar_str, cmp_MD
                        print >> sys.stderr, "\tcmp list:", cmp_list
                        assert False            

                    prev_read_id = read_id
                    # prev_exon = exon

                if num_reads <= 0:
                    continue

                if prev_read_id != None:
                    add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read)

                HLA_counts = [[allele, count] for allele, count in HLA_counts.items()]
                def HLA_count_cmp(a, b):
                    if a[1] != b[1]:
                        return b[1] - a[1]
                    assert a[0] != b[0]
                    if a[0] < b[0]:
                        return -1
                    else:
                        return 1
                HLA_counts = sorted(HLA_counts, cmp=HLA_count_cmp)
                for count_i in range(len(HLA_counts)):
                    count = HLA_counts[count_i]
                    print >> sys.stderr, "\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1])
                    if count_i >= 9:
                        break
                print >> sys.stderr

                def normalize(prob):
                    total = sum(prob.values())
                    for allele, mass in prob.items():
                        prob[allele] = mass / total

                def normalize2(prob, length):
                    total = 0
                    for allele, mass in prob.items():
                        assert allele in length
                        total += (mass / length[allele])
                    for allele, mass in prob.items():
                        assert allele in length
                        prob[allele] = mass / length[allele] / total

                def prob_diff(prob1, prob2):
                    diff = 0.0
                    for allele in prob1.keys():
                        if allele in prob2:
                            diff += abs(prob1[allele] - prob2[allele])
                        else:
                            diff += prob1[allele]
                    return diff

                def HLA_prob_cmp(a, b):
                    if a[1] != b[1]:
                        if a[1] < b[1]:
                            return 1
                        else:
                            return -1
                    assert a[0] != b[0]
                    if a[0] < b[0]:
                        return -1
                    else:
                        return 1

                HLA_prob, HLA_prob_next = {}, {}
                for cmpt, count in HLA_cmpt.items():
                    alleles = cmpt.split('-')
                    for allele in alleles:
                        if allele not in HLA_prob:
                            HLA_prob[allele] = 0.0
                        HLA_prob[allele] += (float(count) / len(alleles))

                """
                assert gene in HLA_lengths
                HLA_length = HLA_lengths[gene]
                """
                HLA_length = {}
                
                # normalize2(HLA_prob, HLA_length)
                normalize(HLA_prob)
                def next_prob(HLA_cmpt, HLA_prob, HLA_length):
                    HLA_prob_next = {}
                    for cmpt, count in HLA_cmpt.items():
                        alleles = cmpt.split('-')
                        alleles_prob = 0.0
                        for allele in alleles:
                            assert allele in HLA_prob
                            alleles_prob += HLA_prob[allele]
                        for allele in alleles:
                            if allele not in HLA_prob_next:
                                HLA_prob_next[allele] = 0.0
                            HLA_prob_next[allele] += (float(count) * HLA_prob[allele] / alleles_prob)
                    # normalize2(HLA_prob_next, HLA_length)
                    normalize(HLA_prob_next)
                    return HLA_prob_next

                diff, iter = 1.0, 0
                while diff > 0.0001 and iter < 1000:
                    HLA_prob_next = next_prob(HLA_cmpt, HLA_prob, HLA_length)
                    diff = prob_diff(HLA_prob, HLA_prob_next)
                    HLA_prob = HLA_prob_next
                    iter += 1

                """
                for allele, prob in HLA_prob.items():
                    allele_len = len(HLAs[gene][allele])
                    HLA_prob[allele] /= float(allele_len)
                normalize(HLA_prob)
                """
                HLA_prob = [[allele, prob] for allele, prob in HLA_prob.items()]

                HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
                success = [False for i in range(len(test_HLA_names))]
                found_list = [False for i in range(len(test_HLA_names))]
                for prob_i in range(len(HLA_prob)):
                    prob = HLA_prob[prob_i]
                    print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, prob[0], prob[1] * 100.0)
                    if prob_i >= 9:
                        break
                print >> sys.stderr

                """
                if len(test_HLA_names) == 2:
                    HLA_prob, HLA_prob_next = {}, {}
                    for cmpt, count in HLA_cmpt.items():
                        alleles = cmpt.split('-')
                        for allele1 in alleles:
                            for allele2 in HLA_names[gene]:
                                if allele1 < allele2:
                                    allele_pair = "%s-%s" % (allele1, allele2)
                                else:
                                    allele_pair = "%s-%s" % (allele2, allele1)
                                if not allele_pair in HLA_prob:
                                    HLA_prob[allele_pair] = 0.0
                                HLA_prob[allele_pair] += (float(count) / len(alleles))

                    if len(HLA_prob) <= 0:
                        continue

                    # Choose top allele pairs
                    def choose_top_alleles(HLA_prob):
                        HLA_prob_list = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
                        HLA_prob_list = sorted(HLA_prob_list, cmp=HLA_prob_cmp)
                        HLA_prob = {}
                        best_prob = HLA_prob_list[0][1]
                        for i in range(len(HLA_prob_list)):
                            allele_pair, prob = HLA_prob_list[i]
                            if prob * 2 <= best_prob:
                                break                        
                            HLA_prob[allele_pair] = prob
                        normalize(HLA_prob)
                        return HLA_prob
                    HLA_prob = choose_top_alleles(HLA_prob)

                    def next_prob(HLA_cmpt, HLA_prob):
                        HLA_prob_next = {}
                        for cmpt, count in HLA_cmpt.items():
                            alleles = cmpt.split('-')
                            prob = 0.0
                            for allele in alleles:
                                for allele_pair in HLA_prob.keys():
                                    if allele in allele_pair:
                                        prob += HLA_prob[allele_pair]
                            for allele in alleles:
                                for allele_pair in HLA_prob.keys():
                                    if not allele in allele_pair:
                                        continue
                                    if allele_pair not in HLA_prob_next:
                                        HLA_prob_next[allele_pair] = 0.0
                                    HLA_prob_next[allele_pair] += (float(count) * HLA_prob[allele_pair] / prob)
                        normalize(HLA_prob_next)
                        return HLA_prob_next

                    diff, iter = 1.0, 0
                    while diff > 0.0001 and iter < 1000:
                        HLA_prob_next = next_prob(HLA_cmpt, HLA_prob)
                        diff = prob_diff(HLA_prob, HLA_prob_next)
                        HLA_prob = HLA_prob_next
                        HLA_prob = choose_top_alleles(HLA_prob)
                        iter += 1

                    HLA_prob = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
                    HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)

                    success = [False]
                    for prob_i in range(len(HLA_prob)):
                        allele_pair, prob = HLA_prob[prob_i]
                        allele1, allele2 = allele_pair.split('-')
                        if best_alleles and prob_i < 1:
                            print >> sys.stdout, "PairModel %s (abundance: %.2f%%)" % (allele_pair, prob * 100.0)
                        if simulation:
                            if allele1 in test_HLA_names and allele2 in test_HLA_names:
                                rank_i = prob_i
                                while rank_i > 0:
                                    if HLA_prob[rank_i-1][1] == prob:                                        
                                        rank_i -= 1
                                    else:
                                        break
                                print >> sys.stderr, "\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (rank_i + 1, allele_pair, prob * 100.0)
                                if rank_i == 0:
                                    success[0] = True
                                break
                        print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, allele_pair, prob * 100.0)
                        if not simulation and prob_i >= 9:
                            break
                    print >> sys.stderr
                """

    # Read variants with clinical significance
    clnsigs = {}
    for line in open("%s.clnsig" % base_fname):
        var_id, var_gene, var_clnsig = line.strip().split('\t')
        clnsigs[var_id] = [var_gene, var_clnsig]

    vars, Var_list = {}, {}
    for line in open("%s.snp" % base_fname):
        var_id, type, chr, left, data = line.strip().split()
        if var_id not in clnsigs:
            continue
        left = int(left)
        if type == "deletion":
            data = int(data)
        vars[var_id] = [chr, left, type, data]
        if chr not in Var_list:
            Var_list[chr] = []
        Var_list[chr].append([left, var_id])

    var_counts = {}

    # Read alignments
    alignview_cmd = ["samtools",
                     "view",
                     "hla_input.bam"]
    bamview_proc = subprocess.Popen(alignview_cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=open("/dev/null", 'w'))

    for line in bamview_proc.stdout:
        cols = line.strip().split()
        read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
        read_seq, qual = cols[9], cols[10]
        flag, pos = int(flag), int(pos)
        pos -= 1
        if pos < 0:
            continue

        if flag & 0x4 != 0:
            continue

        if chr not in Var_list:
            continue

        assert chr in chr_dic
        chr_seq = chr_dic[chr]

        NM, Zs, MD, NH = "", "", "", ""
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith("Zs"):
                Zs = col[5:]
            elif col.startswith("MD"):
                MD = col[5:]
            elif col.startswith("NM"):
                NM = int(col[5:])
            elif col.startswith("NH"):
                NH = int(col[5:])

        assert NH != ""
        NH = int(NH)
        if NH > 1:
            continue

        if NM > num_mismatch:
            continue

        read_vars = []
        if Zs:
            read_vars = Zs.split(',')
        for read_var in read_vars:
            _, _, var_id = read_var.split('|')
            if var_id not in clnsigs:
                continue
            if var_id not in var_counts:
                var_counts[var_id] = [1, 0]
            else:
                var_counts[var_id][0] += 1

        assert MD != ""
        MD_str_pos, MD_len = 0, 0
        read_pos, left_pos = 0, pos
        right_pos = left_pos
        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
        cmp_list = []
        for i in range(len(cigars)):
            cigar_op, length = cigars[i]
            if cigar_op == 'M':
                chr_var_list = Var_list[chr]
                var_idx = lower_bound(chr_var_list, right_pos)
                while var_idx < len(chr_var_list):
                    var_pos, var_id = chr_var_list[var_idx]
                    if var_pos >= right_pos + length:
                        break
                    if var_pos >= right_pos:
                        assert var_id in vars
                        _, _, var_type, var_data = vars[var_id]
                        contradict = False
                        if var_type == "single":
                            contradict = (read_seq[read_pos + var_pos - right_pos] == chr_seq[var_pos])
                        elif var_type == "insertion":
                            contradict = (right_pos < var_pos)
                        else:
                            contradict = True
                        if contradict:
                            if var_id not in var_counts:
                                var_counts[var_id] = [0, 1]
                            else:
                                var_counts[var_id][1] += 1
                    
                    var_idx += 1
                    
            if cigar_op in "MND":
                right_pos += length

            if cigar_op in "MIS":
                read_pos += length

    for var_id, counts in var_counts.items():
        if counts[0] < 2: # or counts[0] * 3 < counts[1]:
            continue
        assert var_id in vars
        var_chr, var_left, var_type, var_data = vars[var_id]
        assert var_id in clnsigs
        var_gene, var_clnsig = clnsigs[var_id]
        print >> sys.stderr, "\t\t\t%s %s: %s:%d %s %s (%s): %d-%d" % \
                (var_gene, var_id, var_chr, var_left, var_type, var_data, var_clnsig, counts[0], counts[1])


                
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='HISAT2 genotyping')
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--base-name",
                        dest="base_fname",
                        type=str,
                        default="genotype_genome",
                        help="base filename for genotype genome")
    parser.add_argument('-f',
                        dest='fastq',
                        action='store_false',
                        help='FASTA file')    
    parser.add_argument("-U",
                        dest="read_fname_U",
                        type=str,
                        default="",
                        help="filename for single-end reads")
    parser.add_argument("-1",
                        dest="read_fname_1",
                        type=str,
                        default="",
                        help="filename for paired-end reads")
    parser.add_argument("-2",
                        dest="read_fname_2",
                        type=str,
                        default="",
                        help="filename for paired-end reads")    
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--simulate-interval",
                        dest="simulate_interval",
                        type=int,
                        default=1,
                        help="Reads simulated at every these base pairs (default: 1)")
    parser.add_argument("--num-mismatch",
                        dest="num_mismatch",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument("--daehwan-debug",
                        dest="daehwan_debug",
                        type=str,
                        default="",
                        help="e.g., test_id:10,read_id:10000,basic_test")

    args = parser.parse_args()
    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    daehwan_debug = {}
    if args.daehwan_debug != "":
        for item in args.daehwan_debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                daehwan_debug[key] = value
            else:
                daehwan_debug[item] = 1

    if args.read_fname_U != "":
        read_fnames = [args.read_fname_U]
    else:
        if args.read_fname_1 == "" or args.read_fname_2 == "":
            print >> sys.stderr, "Error: please specify read file names correctly: -U or -1 and -2"
            sys.exit(1)
        read_fnames = [args.read_fname_1, args.read_fname_2] 

    random.seed(1)
    genotype(args.reference_type,
             args.base_fname,
             args.fastq,
             read_fnames,
             args.threads,
             args.simulate_interval,
             args.num_mismatch,
             args.verbose,
             daehwan_debug)
