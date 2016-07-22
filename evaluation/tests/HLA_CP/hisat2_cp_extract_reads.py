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
import glob
from argparse import ArgumentParser, FileType


"""
"""
def reverse_complement(seq):
    comp_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    rc_seq = ""
    for s in reversed(seq):
        if s in comp_table:
            rc_seq += comp_table[s]
        else:
            rc_seq += s
    return rc_seq
            

"""
"""
def parallel_work(pids, 
                  work, 
                  ex_path,
                  fq_fname_base, 
                  fq_fname, 
                  fq_fname2, 
                  reference_type, 
                  ranges):
    child = -1
    for i in range(len(pids)):
        if pids[i] == 0:
            child = i
            break

    while child == -1:
        status = os.waitpid(0, 0)
        for i in range(len(pids)):
            if status[0] == pids[i]:
                child = i
                pids[i] = 0
                break

    child_id = os.fork()
    if child_id == 0:
        work(ex_path,
             fq_fname_base, 
             fq_fname, 
             fq_fname2, 
             reference_type, 
             ranges)
        os._exit(os.EX_OK)
    else:
        # print >> sys.stderr, '\t\t>> thread %d: %d' % (child, child_id)
        pids[child] = child_id

        
"""
"""
def wait_pids(pids):
    for pid in pids:
        if pid > 0:
            os.waitpid(pid, 0)
            

"""
"""
def extract_reads(base_fname,
                  reference_type,
                  hla_list,
                  partial,
                  DRB1_ref,
                  threads,
                  job_range,
                  verbose):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(extract_reads))
    ex_path = "../../.."

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        os.system("git clone https://github.com/jrob119/IMGTHLA.git")

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

    if reference_type == "gene":
        # Extract HLA variants, backbone sequence, and other sequeces
        HLA_fnames = ["hla_backbone.fa",
                      "hla_sequences.fa",
                      "hla.ref",
                      "hla.snp",
                      "hla.haplotype",
                      "hla.link"]

        if not check_files(HLA_fnames):
            extract_hla_script = os.path.join(ex_path, "hisat2_extract_HLA_vars.py")
            extract_cmd = [extract_hla_script]
            extract_cmd += ["--reference-type", reference_type,
                                "--hla-list", ','.join(hla_list)]
            if DRB1_ref:
                extract_cmd += ["--DRB1-REF"]
            if partial:
                extract_cmd += ["--partial"]
            extract_cmd += ["--inter-gap", "30",
                            "--intra-gap", "50"]
            print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
            proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()
            if not check_files(HLA_fnames):
                print >> sys.stderr, "Error: extract_HLA_vars failed!"
                sys.exit(1)

        # Build HISAT2 graph indexes based on the above information
        HLA_hisat2_graph_index_fnames = ["hla.graph.%d.ht2" % (i+1) for i in range(8)]
        if not check_files(HLA_hisat2_graph_index_fnames):
            hisat2_build = os.path.join(ex_path, "hisat2-build")
            build_cmd = [hisat2_build,
                         "--snp", "hla.snp",
                         "--haplotype", "hla.haplotype",
                         "hla_backbone.fa",
                         "hla.graph"]
            print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
            proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()        
            if not check_files(HLA_hisat2_graph_index_fnames):
                print >> sys.stderr, "Error: indexing HLA failed!  Perhaps, you may have forgotten to build hisat2 executables?"
                sys.exit(1)
    else:
        assert reference_type == "genome"
        genotype_fnames = ["%s.fa" % base_fname,
                           "%s.gene" % base_fname,
                           "%s.snp" % base_fname,
                           "%s.haplotype" % base_fname,
                           "%s.link" % base_fname,
                           "%s.coord" % base_fname,
                           "%s.clnsig" % base_fname]
        # hisat2 graph index files
        genotype_fnames += ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]
        if not check_files(genotype_fnames):        
            build_script = os.path.join(ex_path, "hisat2_build_genotype_genome.py")
            build_cmd = [build_script]
            if partial:
                build_cmd += ["--partial"]
            build_cmd += ["--inter-gap", "30",
                          "--intra-gap", "50"]
            build_cmd += ["--threads", "4"]
            build_cmd += ["--no-clinvar"]
            build_cmd += ["genome.fa", "genotype_genome"]
            print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
            proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()
            if not check_files(genotype_fnames):
                print >> sys.stderr, "Error: indexing genotype genome failed!  Perhaps, you may have forgotten to build hisat2 executables?"
                sys.exit(1)

    ranges = []
    genes, gene_loci = {}, {}
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

    if not os.path.exists("CP"):
        os.mkdir("CP")

    # Extract reads
    fq_fnames = glob.glob("genomes/*.1.fq.gz")
    count = 0
    pids = [0 for i in range(threads)]
    for file_i in range(len(fq_fnames)):
        fq_fname = fq_fnames[file_i]

        if job_range[1] > 1:
            if job_range[0] != (file_i % job_range[1]):
                continue

        fq_fname_base = fq_fname.split('/')[-1]
        fq_fname_base = fq_fname_base.split('.')[0]
        fq_fname2 = "genomes/%s.2.fq.gz" % fq_fname_base
        if not os.path.exists(fq_fname2):
            print >> sys.stderr, "%s does not exist." % fq_fname2
            continue
        if os.path.exists("CP/%s.extracted.fq.1.gz" % fq_fname_base):
            continue        
        count += 1

        print >> sys.stderr, "\t%d: Extracting reads from %s" % (count, fq_fname_base)
        def work(ex_path,
                 fq_fname_base,
                 fq_fname, 
                 fq_fname2, 
                 reference_type, 
                 ranges):
            hisat2 = os.path.join(ex_path, "hisat2")
            aligner_cmd = [hisat2]
            if reference_type == "gene":
                aligner_cmd += ["--al-conc-disc-gz", "CP/%s.extracted.fq.gz" % fq_fname_base]
                aligner_cmd += ["-x", "hla.graph"]
            else:
                aligner_cmd += ["-x", "genotype_genome"]
            aligner_cmd += ["-1", fq_fname,
                            "-2", fq_fname2]
            # print >> sys.stderr, "\t\trunning", ' '.join(aligner_cmd)
            if reference_type == "gene": 
                align_proc = subprocess.Popen(aligner_cmd,
                                              stdout=open("/dev/null", 'w'),
                                              stderr=open("/dev/null", 'w'))
                align_proc.communicate()
            else:
                assert reference_type == "genome"
                align_proc = subprocess.Popen(aligner_cmd,
                                              stdout=subprocess.PIPE,
                                              stderr=open("/dev/null", 'w'))
                # LP6005041-DNA_A01.extracted.fq.1.gz
                gzip1_proc = subprocess.Popen(["gzip"],
                                              stdin=subprocess.PIPE,
                                              stdout=open("CP/%s.extracted.fq.1.gz" % fq_fname_base, 'w'),
                                              stderr=open("/dev/null", 'w'))

                # LP6005041-DNA_A01.extracted.fq.2.gz
                gzip2_proc = subprocess.Popen(["gzip"],
                                              stdin=subprocess.PIPE,
                                              stdout=open("CP/%s.extracted.fq.2.gz" % fq_fname_base, 'w'),
                                              stderr=open("/dev/null", 'w'))
                prev_read_name, extract_read, read1, read2 = "", False, [], []
                for line in align_proc.stdout:
                    if line.startswith('@'):
                        continue
                    line = line.strip()
                    cols = line.split()
                    read_name, flag, chr, pos, mapQ, cigar, _, _, _, read, qual = cols[:11]
                    flag, pos = int(flag), int(pos)
                    """
                    strand = '-' if flag & 0x10 else '+'                   
                    AS = ""
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("AS"):
                            AS = int(col[5:])
                    """

                    if read_name != prev_read_name:
                        if extract_read:
                            gzip1_proc.stdin.write("@%s\n" % prev_read_name)
                            gzip1_proc.stdin.write("%s\n" % read1[0])
                            gzip1_proc.stdin.write("+\n")
                            gzip1_proc.stdin.write("%s\n" % read1[1])
                            gzip2_proc.stdin.write("@%s\n" % prev_read_name)
                            gzip2_proc.stdin.write("%s\n" % read2[0])
                            gzip2_proc.stdin.write("+\n")
                            gzip2_proc.stdin.write("%s\n" % read2[1])

                        prev_read_name, extract_read, read1, read2 = read_name, False, [], []
                        
                    if flag & 0x4 == 0:
                        for loci in gene_loci.values():
                            _, loci_chr, loci_left, loci_right = loci
                            if chr == loci_chr and pos >= loci_left and pos < loci_right:
                                extract_read = True
                                break

                    if flag & 0x40: # left read
                        if not read1:
                            if flag & 0x10: # reverse complement
                                read1 = [reverse_complement(read), qual[::-1]]
                            else:
                                read1 = [read, qual]
                    else:
                        assert flag & 0x80 # right read
                        if flag & 0x10: # reverse complement
                            read2 = [reverse_complement(read), qual[::-1]]
                        else:
                            read2 = [read, qual]

                if extract_read:
                    gzip1_proc.stdin.write("@%s\n" % prev_read_name)
                    gzip1_proc.stdin.write("%s\n" % read1[0])
                    gzip1_proc.stdin.write("+\n")
                    gzip1_proc.stdin.write("%s\n" % read1[1])
                    gzip2_proc.stdin.write("@%s\n" % prev_read_name)
                    gzip2_proc.stdin.write("%s\n" % read2[0])
                    gzip2_proc.stdin.write("+\n")
                    gzip2_proc.stdin.write("%s\n" % read2[1])                            
                        
                gzip1_proc.stdin.close()
                gzip2_proc.stdin.close()                        

        if threads <= 1:
            work(ex_path, 
                 fq_fname_base, 
                 fq_fname, 
                 fq_fname2, 
                 reference_type, 
                 ranges)
        else:
            parallel_work(pids, 
                          work, 
                          ex_path, 
                          fq_fname_base, 
                          fq_fname, 
                          fq_fname2, 
                          reference_type, 
                          ranges)

    if threads > 1:
        wait_pids(pids)


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping for Platinum Genomes')
    parser.add_argument("--base_fname",
                        dest="base_fname",
                        type=str,
                        default="genotype_genome",
                        help="base filename for genotype genome")
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="genome",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--hla-list",
                        dest="hla_list",
                        type=str,
                        default="A,B,C,DQA1,DQB1,DRB1",
                        help="A comma-separated list of HLA genes (default: A,B,C,DQA1,DQB1,DRB1)")
    parser.add_argument('--partial',
                        dest='partial',
                        action='store_true',
                        help='Include partial alleles (e.g. A_nuc.fasta)')
    parser.add_argument("--DRB1-REF",
                        dest="DRB1_ref",
                        action='store_true',
                        help="Some DRB1 alleles seem to include vector sequences, so use this option to avoid including them")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--job-range",
                        dest="job_range",
                        type=str,
                        default="0,1",
                        help="two numbers (e.g. 1,3)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()

    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    args.hla_list = args.hla_list.split(',')
    job_range = []
    for num in args.job_range.split(','):
        job_range.append(int(num))
    extract_reads(args.base_fname,
                  args.reference_type,
                  args.hla_list,
                  args.partial,
                  args.DRB1_ref,
                  args.threads,
                  job_range,
                  args.verbose)
