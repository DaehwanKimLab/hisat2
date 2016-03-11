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
import inspect
import random
from argparse import ArgumentParser, FileType


"""
"""
def test_BRCA_genotyping(reference_type,
                         brca_list,
                         aligners,
                         read_fname,
                         alignment_fname,
                         threads,
                         simulate_interval,
                         enable_coverage,
                         num_mismatch,
                         verbose,
                         daehwan_debug):
    # Is it simulation?
    simulation = (not read_fname and not alignment_fname)
    
    # File location for ClinVar
    clinvar_url_base = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
    clinvar_fname = "clinvar_20160203.vcf.gz"

    if not os.path.exists(clinvar_fname):
        os.system("wget %s/%s" % (clinvar_url_base, clinvar_fname))
        assert os.path.exists(clinvar_fname)
    
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_BRCA_genotyping))
    ex_path = os.path.dirname(curr_script)

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")

    # Check if the pre-existing files (hla*) are compatible with the current parameter setting
    if os.path.exists("brca.ref"):
        left = 0
        BRCA_genes = set()
        for line in open("brca.ref"):
            BRCA_name, chr, left, _ = line.strip().split()
            BRCA_gene = BRCA_name.split('*')[0]
            BRCA_genes.add(BRCA_gene)
            left = int(left)
        delete_brca_files = False
        if reference_type == "gene":
            if left > 0:
                delete_hla_files = True
        elif reference_type == "chromosome":
            if left == 0:
                delete_hla_files = True
        else:
            assert reference_type == "genome"
        if not set(brca_list).issubset(BRCA_genes):
            delete_brca_files = True
        if delete_brca_files:
            os.system("rm brca*")

    # Extract BRCA variants, backbone sequence, and other sequeces
    BRCA_fnames = ["brca_backbone.fa",
                   "brca.ref",
                   "brca.snp",
                   "brca.haplotype",
                   "brca.clnsig"]

    if not check_files(BRCA_fnames):
        extract_brca_script = os.path.join(ex_path, "hisat2_extract_snps_haplotypes_VCF.py")
        extract_cmd = [extract_brca_script,
                       "genome.fa",
                       "--base", "brca",
                       "--reference-type", "gene",
                       "--genotype-vcf", clinvar_fname,
                       "--genotype-gene-list", ','.join(brca_list),
                       "--extra-files"]
        extract_cmd += ["--inter-gap", "30",
                        "--intra-gap", "50"]

        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        if not check_files(BRCA_fnames):
            print >> sys.stderr, "Error: extract_BRCA_vars failed!"
            sys.exit(1)

    # Build HISAT2 graph indexes based on the above information
    HLA_hisat2_graph_index_fnames = ["brca.graph.%d.ht2" % (i+1) for i in range(8)]
    if not check_files(HLA_hisat2_graph_index_fnames):
        hisat2_build = os.path.join(ex_path, "hisat2-build")
        build_cmd = [hisat2_build,
                     "-p", str(threads),
                     "--snp", "brca.snp",
                     "--haplotype", "brca.haplotype",
                     "brca_backbone.fa",
                     "brca.graph"]
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_hisat2_graph_index_fnames):
            print >> sys.stderr, "Error: indexing BRCA genes failed!  Perhaps, you may have forgotten to build hisat2 executables?"
            sys.exit(1)

    """
    # Build HISAT2 linear indexes based on the above information
    HLA_hisat2_linear_index_fnames = ["hla.linear.%d.ht2" % (i+1) for i in range(8)]
    if reference_type == "gene" and not check_files(HLA_hisat2_linear_index_fnames):
        hisat2_build = os.path.join(ex_path, "hisat2-build")
        build_cmd = [hisat2_build,
                     "hla_backbone.fa,hla_sequences.fa",
                     "hla.linear"]
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_hisat2_graph_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!"
            sys.exit(1)

    # Build Bowtie2 indexes based on the above information
    HLA_bowtie2_index_fnames = ["hla.%d.bt2" % (i+1) for i in range(4)]
    HLA_bowtie2_index_fnames += ["hla.rev.%d.bt2" % (i+1) for i in range(2)]
    if reference_type == "gene" and not check_files(HLA_bowtie2_index_fnames):
        build_cmd = ["bowtie2-build",
                     "hla_backbone.fa,hla_sequences.fa",
                     "hla"]
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_bowtie2_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!"
            sys.exit(1)
    """

    # Read BRCA variants
    Vars, Var_list = {}, {}
    for line in open("brca.snp"):
        var_id, var_type, allele, pos, data = line.strip().split('\t')
        pos = int(pos)
        if reference_type != "gene":
            allele, dist = None, 0
            for tmp_gene, values in refHLA_loci.items():
                allele_name, chr, left, right = values
                if allele == None:
                    allele = allele_name
                    dist = abs(pos - left)
                else:
                    if dist > abs(pos - left):
                        allele = allele_name
                        dist = abs(pos - left)
            
        gene = allele
        if not gene in Vars:
            Vars[gene] = {}
            assert not gene in Var_list
            Var_list[gene] = []
            
        assert not var_id in Vars[gene]
        left = 0
        if reference_type != "gene":
            _, _, left, _ = refHLA_loci[gene]
        Vars[gene][var_id] = [var_type, pos - left, data]
        Var_list[gene].append([pos - left, var_id])

    for gene, in_var_list in Var_list.items():
        Var_list[gene] = sorted(in_var_list)
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

    # Read BRCA variants' clinical significance
    Vars_CLNSIG = {}
    for line in open("brca.clnsig"):
        var_id, CLNSIG = line.strip().split('\t')
        assert var_id not in Vars_CLNSIG
        Vars_CLNSIG[var_id] = CLNSIG

    BRCAs = {}
    def read_BRCA_alleles(fname, BRCAs):
        for line in open(fname):
            if line.startswith(">"):
                BRCA_gene = line.strip().split()[0][1:]
                if not BRCA_gene in BRCAs:
                    BRCAs[BRCA_gene] = ""
            else:
                BRCAs[BRCA_gene] += line.strip()
        return BRCAs
    if reference_type == "gene":
        read_BRCA_alleles("brca_backbone.fa", BRCAs)
    # read_BRCA_alleles("brca_sequences.fa", BRCAs)

    # Test BRCA genotyping
    test_list = []
    if simulation:
        test_passed = {}
        test_list = []
        genes = list(set(brca_list) & BRCA_genes)
        for gene in genes:
            for var_id in Vars[gene].keys():
                if var_id not in Vars_CLNSIG:
                    continue
                test_list.append([gene, [var_id]])
    else:
        test_list = [brca_list]

    for test_i in range(len(test_list)):
        if "test_id" in daehwan_debug:
            daehwan_test_ids = daehwan_debug["test_id"].split('-')
            if str(test_i + 1) not in daehwan_test_ids:
                continue

        print >> sys.stderr, "Test %d" % (test_i + 1)
        gene, test_var_ids = test_list[test_i]

        # daehwan - for debugging purposes
        # test_HLA_list = [["A*11:50Q", "A*11:01:01:01", "A*01:01:01:01"]]
        print >> sys.stderr, "\t%s" % (gene)
        for test_var_id in test_var_ids:
            assert test_var_id in Vars[gene]
            var_type, var_pos, var_data = Vars[gene][test_var_id]
            
            print >> sys.stderr, "\t\t%s:" % (test_var_id), var_type, var_pos, var_data, Vars_CLNSIG[test_var_id]

        var_true_counts = {}
        if simulation:
            BRCA_seq = BRCAs[gene]
            BRCA_reads_1, BRCA_reads_2 = [], []
            # Simulate reads from two HLA alleles
            def simulate_reads(seq, test_vars, simulate_interval = 1, frag_len = 250, read_len = 100):
                comp_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
                reads_1, reads_2 = [], []
                v = 0
                for i in range(0, len(seq) - frag_len + 1, simulate_interval):
                    while v < len(test_vars):
                        var_type, var_pos, var_data = test_vars[v]
                        if var_type == 'D':
                            var_pos += (int(var_data) - 1)
                        if var_pos >= i:
                            break
                        v += 1
                        
                    def simulate_read(left, right, test_var):
                        include_var = False
                        if test_var != None:
                            var_type, var_pos, var_data = test_var
                            var_pos2 = var_pos
                            if var_type == 'deletion':
                                var_pos2 += int(var_data)
                            if var_pos > left and var_pos2 < right:
                                include_var = True

                        if include_var:
                            if var_type == 'single':
                                return seq[left:var_pos] + var_data + seq[var_pos+1:right]
                            elif var_type == 'deletion':
                                return seq[left:var_pos] + seq[var_pos2:right]
                            else:
                                assert var_type == 'insertion'
                                return seq[left:var_pos] + var_data + seq[var_pos:right]
                        else:
                            return seq[left:right]

                    test_var = None
                    if v < len(test_vars):
                        test_var = test_vars[v]
                    reads_1.append(simulate_read(i, i+read_len, test_var))
                    tmp_read_2 = simulate_read(i+frag_len-read_len,i+frag_len, test_var)
                    tmp_read_2 = reversed(tmp_read_2)
                    read_2 = ""
                    for s in tmp_read_2:
                        if s in comp_table:
                            read_2 += comp_table[s]
                        else:
                            read_2 += s
                    reads_2.append(read_2)
                return reads_1, reads_2

            test_vars = []
            for test_var_id in test_var_ids:
                assert test_var_id in Vars[gene]
                test_vars.append(Vars[gene][test_var_id])

            tmp_reads_1, tmp_reads_2 = simulate_reads(BRCA_seq, test_vars, simulate_interval)
            BRCA_reads_1 += tmp_reads_1
            BRCA_reads_2 += tmp_reads_2

            # Write reads into a fasta read file
            def write_reads(reads, idx):
                read_file = open('brca_input_%d.fa' % idx, 'w')
                for read_i in range(len(reads)):
                    print >> read_file, ">%d" % (read_i + 1)
                    print >> read_file, reads[read_i]
                read_file.close()
            write_reads(BRCA_reads_1, 1)
            write_reads(BRCA_reads_2, 2)

        for aligner, index_type in aligners:
            if index_type == "graph":
                print >> sys.stderr, "\n\t\t%s %s on %s" % (aligner, index_type, reference_type)
            else:
                print >> sys.stderr, "\n\t\t%s %s" % (aligner, index_type)

            if alignment_fname == "":
                # Align reads, and sort the alignments into a BAM file
                if aligner == "hisat2":
                    hisat2 = os.path.join(ex_path, "hisat2")
                    aligner_cmd = [hisat2,
                                   "--no-unal",
                                   "--mm"]
                    if index_type == "linear":
                        aligner_cmd += ["-k", "10"]
                    aligner_cmd += ["-x", "brca.%s" % index_type]
                elif aligner == "bowtie2":
                    aligner_cmd = [aligner,
                                   "--no-unal",
                                   "-k", "10",
                                   "-x", "brca"]
                else:
                    assert False
                if simulation:
                    if "test_id" in daehwan_debug:
                        aligner_cmd += ["-f", "brca_input_1.fa"]
                    else:
                        aligner_cmd += ["-f",
                                        "-1", "brca_input_1.fa",
                                        "-2", "brca_input_2.fa"]
                else:
                    assert len(read_fname) in [1,2]
                    aligner_cmd += ["-p", str(threads)]
                    if len(read_fname) == 1:
                        aligner_cmd += [read_fname[0]]
                    else:
                        aligner_cmd += ["-1", "%s" % read_fname[0],
                                        "-2", "%s" % read_fname[1]]

                align_proc = subprocess.Popen(aligner_cmd,
                                              stdout=subprocess.PIPE,
                                              stderr=open("/dev/null", 'w'))

                sambam_cmd = ["samtools",
                              "view",
                              "-bS",
                              "-"]
                sambam_proc = subprocess.Popen(sambam_cmd,
                                               stdin=align_proc.stdout,
                                               stdout=open("brca_input_unsorted.bam", 'w'),
                                               stderr=open("/dev/null", 'w'))
                sambam_proc.communicate()
                if index_type == "graph":
                    bamsort_cmd = ["samtools",
                                   "sort",
                                   "brca_input_unsorted.bam",
                                   "brca_input"]
                    bamsort_proc = subprocess.Popen(bamsort_cmd,
                                                    stderr=open("/dev/null", 'w'))
                    bamsort_proc.communicate()

                    bamindex_cmd = ["samtools",
                                    "index",
                                    "brca_input.bam"]
                    bamindex_proc = subprocess.Popen(bamindex_cmd,
                                                     stderr=open("/dev/null", 'w'))
                    bamindex_proc.communicate()

                    os.system("rm brca_input_unsorted.bam")            
                else:
                    os.system("mv brca_input_unsorted.bam brca_input.bam")

            for tt in [0]:
                ref_seq = BRCAs[gene]

                # Read alignments
                alignview_cmd = ["samtools",
                                 "view"]
                if alignment_fname == "":
                    alignview_cmd += ["brca_input.bam"]
                else:
                    if not os.path.exists(alignment_fname + ".bai"):
                        os.system("samtools index %s" % alignment_fname)
                    alignview_cmd += [alignment_fname]
                base_locus = 0
                if index_type == "graph":
                    if reference_type == "gene":
                        alignview_cmd += ["%s" % gene]
                    else:
                        assert reference_type in ["chromosome", "genome"]
                        _, chr, left, right = refHLA_loci[gene]
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
                else:
                    alignview_proc = subprocess.Popen(alignview_cmd,
                                                 stdout=subprocess.PIPE,
                                                 stderr=open("/dev/null", 'w'))

                # Count alleles
                var_test_counts = {}
                num_reads, total_read_len = 0, 0
                prev_read_id = None
                if index_type == "graph":
                    # Cigar regular expression
                    cigar_re = re.compile('\d+\w')
                    for line in alignview_proc.stdout:
                        cols = line.strip().split()
                        read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
                        read_seq = cols[9]
                        num_reads += 1
                        total_read_len += len(read_seq)
                        flag, pos = int(flag), int(pos)
                        pos -= (base_locus + 1)
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
                                # Update coverage
                                if enable_coverage:
                                    if right_pos + length < len(coverage):
                                        coverage[right_pos] += 1
                                        coverage[right_pos + length] -= 1
                                    elif right_pos < len(coverage):
                                        coverage[right_pos] += 1
                                        coverage[-1] -= 1

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

                        if right_pos > len(ref_seq):
                            continue

                        # Decide which allele(s) a read most likely came from
                        # also sanity check - read length, cigar string, and MD string
                        def add_count(var_id, num):
                            if var_id not in var_test_counts:
                                var_test_counts[var_id] = 0
                            var_test_counts[var_id] += num
                            
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
                                        elif var_type == "deletion":
                                            del_len = int(var_data)
                                            if ref_pos < var_pos and ref_pos + length > var_pos + del_len:
                                                # Check if this might be one of the two tandem repeats (the same left coordinate)
                                                cmp_left, cmp_right = cmp[1], cmp[1] + cmp[2]
                                                test1_seq1 = ref_seq[cmp_left:cmp_right]
                                                test1_seq2 = ref_seq[cmp_left:var_pos] + ref_seq[var_pos + del_len:cmp_right + del_len]
                                                # Check if this happens due to small repeats (the same right coordinate - e.g. 19 times of TTTC in DQA1*05:05:01:02)
                                                cmp_left -= read_pos
                                                cmp_right += (len(read_seq) - read_pos - cmp[2])
                                                test2_seq1 = ref_seq[cmp_left+int(var_data):cmp_right]
                                                test2_seq2 = ref_seq[cmp_left:var_pos] + ref_seq[var_pos+int(var_data):cmp_right]
                                                if test1_seq1 != test1_seq2 and test2_seq1 != test2_seq2:
                                                    add_count(var_id, -1)
                                        else:
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
                                                add_count(var_id, 1)
                                    var_idx += 1

                                cmp_MD += ("%d%s" % (MD_match_len, ref_seq[ref_pos]))
                                MD_match_len = 0
                                cigar_match_len += 1
                                read_pos += 1
                                ref_pos += 1
                            elif type == "insertion":
                                ins_seq = read_seq[read_pos:read_pos+length]
                                var_idx = lower_bound(Var_list[gene], ref_pos)
                                while var_idx < len(Var_list[gene]):
                                    var_pos, var_id = Var_list[gene][var_idx]
                                    if ref_pos < var_pos:
                                        break
                                    if ref_pos == var_pos:
                                        var_type, _, var_data = Vars[gene][var_id]
                                        if var_type == "insertion":                                
                                            if var_data == ins_seq:
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
                                    last_bp = ref_seq[temp_ref_pos + del_len - 1]
                                    prev_bp = ref_seq[temp_ref_pos - 1]
                                    if last_bp != prev_bp:
                                        break
                                    temp_ref_pos -= 1
                                var_idx = lower_bound(Var_list[gene], temp_ref_pos)
                                while var_idx < len(Var_list[gene]):
                                    var_pos, var_id = Var_list[gene][var_idx]
                                    if temp_ref_pos < var_pos:
                                        first_bp = ref_seq[temp_ref_pos]
                                        next_bp = ref_seq[temp_ref_pos + del_len]
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
                                                add_count(var_id, 1)
                                    var_idx += 1
                                if cigar_match_len > 0:
                                    cmp_cigar_str += ("%dM" % cigar_match_len)
                                    cigar_match_len = 0
                                cmp_MD += ("%d" % MD_match_len)
                                MD_match_len = 0
                                cmp_cigar_str += ("%dD" % length)
                                cmp_MD += ("^%s" % ref_seq[ref_pos:ref_pos+length])
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

                    if num_reads <= 0:
                        continue

                for var_id, count in var_test_counts.items():
                    if count <= 0:
                        continue
                    print >> sys.stderr, "\t\t\t%s: %d" % (var_id, count)

    if simulation:
        for aligner_type, passed in test_passed.items():
            print >> sys.stderr, "%s:\t%d/%d passed (%.2f%%)" % (aligner_type, passed, len(test_list), passed * 100.0 / len(test_list))
    
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping')
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--brca-list",
                        dest="brca_list",
                        type=str,
                        default="BRCA1,BRCA2",
                        help="A comma-separated list of BRCA genes (default: BRCA1,BRCA2)")
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="hisat2.graph",
                        help="A comma-separated list of aligners (default: hisat2.graph)")
    parser.add_argument("--reads",
                        dest="read_fname",
                        type=str,
                        default="",
                        help="Fastq read file name")
    parser.add_argument("--alignment",
                        dest="alignment_fname",
                        type=str,
                        default="",
                        help="BAM file name")
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
    parser.add_argument("--coverage",
                        dest="coverage",
                        action='store_true',
                        help="Experimental purpose (assign reads based on coverage)")
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
    args.brca_list = args.brca_list.split(',')
    if args.aligners == "":
        print >> sys.stderr, "Error: --aligners must be non-empty."
        sys.exit(1)    
    args.aligners = args.aligners.split(',')
    for i in range(len(args.aligners)):
        args.aligners[i] = args.aligners[i].split('.')
    if args.read_fname:
        args.read_fname = args.read_fname.split(',')
    else:
        args.read_fname = []
    if args.alignment_fname != "" and \
            not os.path.exists(args.alignment_fname):
        print >> sys.stderr, "Error: %s doesn't exist." % args.alignment_fname
        sys.exit(1)
    daehwan_debug = {}
    if args.daehwan_debug != "":
        for item in args.daehwan_debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                daehwan_debug[key] = value
            else:
                daehwan_debug[item] = 1

    random.seed(1)
    test_BRCA_genotyping(args.reference_type,
                         args.brca_list,
                         args.aligners,
                         args.read_fname,
                         args.alignment_fname,
                         args.threads,
                         args.simulate_interval,
                         args.coverage,
                         args.num_mismatch,
                         args.verbose,
                         daehwan_debug)
