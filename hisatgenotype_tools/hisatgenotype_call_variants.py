#!/usr/bin/env python

import sys, os, subprocess, re
from multiprocessing import Pool
from argparse import ArgumentParser
import hisatgenotype_typing_common as typing_common
import hisatgenotype_args as arguments

def record_variants(aligner,
                    is_fastq,
                    reference,
                    use_graph,
                    read_fnames,
                    alignment_fname,
                    first_align,
                    keep_align,
                    threads):
    if first_align:
        if reference == "":
            print >> sys.stderr, "Error: --ref-genome or -x option not set"
            exit(1)

        fname_prefix = read_fnames[0].split('.')[0]
        bam_fname = fname_prefix + '.bam'

        cmd = [aligner]
        if aligner == "hisat2" or alinger == "bowtie2":
            cmd += ['-x', reference, '-p', threads]
            if len(read_fnames) > 1:
                cmd += ['-1', read_fnames[0], '-2', read_fnames[1]]
            else:
                cmd += ['-U', read_fnames[0]]
        elif alinger == "bwa":
            cmd += ["mem", reference]
            cmd += read_fnames
            
        cmd += ['|', 'samtools', 'view', '-Sb', '-', '>', bam_fname] 
        os.system(" ".join(cmd))
        alignment_fname.append(bam_fname)

    pull_align = ['samtools', 'view', alignment_fname[0]]
    proc = subprocess.Popen(pull_align,
                            stdout=subprocess.PIPE,
                            stderr=open("/dev/null", "w"))

    snp_dict = {}
    cigar_re = re.compile('\d+\w')
    for line in proc.stdout:
        line = line.strip()
        cols = line.split()
        _ , flag, chrom , pos, mapq, cigar_str, _ , _ , _ , seq, qual = cols[:11] # All standard fields from SAM
        flag, pos, mapq, read_len= int(flag), int(pos), int(mapq), len(seq)
        
        if flag & 0x4 != 0:
            continue

        md, zs, nh = '', '', ''
        for itr in range(11, len(cols)):
            opt, _ , value = cols[itr].split(':')
            if opt == "MD":
                md = value # MD field is string of missmatches
            elif opt == "Zs":
                zs = value # Zs field from Hisat2 contains known snps the read aligned to
            elif opt == "NH" or opt == "X0":
                nh = int(value) # NH in hisat2 and bowtie is number of mapped locations X0 in BWA is number of best locations

        assert md and nh

        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]

        # Skip all reads that don't have a variant
        if all(cigar[0] in "MSH" for cigar in cigars) and md.isdigit():
            md = int(md)
            for cigar in cigars:
                if cigar[0] == 'M' and cigar[1] == md:
                    continue

        # Generate lists of nucleotides and scores
        seq_qual = [[], [], []]
        pval_read = 1-(10**(-mapq/10))
        # for Substitution
        for itr in range(len(seq)):
            qual_score = float(ord(qual[itr])) - 33
            pval_base = 1-(10**(-qual_score/10))
            base_weight = (pval_base*pval_read)/nh

            seq_qual[0].append(seq[itr])
            seq_qual[1].append(base_weight)
            seq_qual[2].append("-")

        # for indel
        seq_qual[0].append('indel')
        seq_qual[1].append(pval_read/nh)

        md = re.split('(\d+)', md)
        md_str = []
        for field in md:
            if field.isdigit():
                for itr in range(int(field)):
                    md_str.append("-")
            else:
                md_str.append(field)

        var_str = []
        for itr in range(len(cigars)):
            cigar_op, lenth = cigars[itr]
            md_pos = 0
            if cigar_op in "SH":
                for jtr in range(length):
                    cigar_str.append(cigar_op)
            if cigar_op == "M":
                
 
if __name__ == '__main__':
    print "This script is not working yet"   
    exit(0)

    parser = ArgumentParser(
        description='HISAT-Genotype Call Variants')
    arguments.args_aligner_inputs(parser,
                                  keep=True) # Add option to keep alignments
    parser.add_argument('-x', '--ref-genome',
                        dest="reference",
                        type=str,
                        default = "",
                        help="Name of reference to use with aligner of choice")
    arguments.args_bamfile(parser)
    arguments.args_set_alinger(parser,
                               missmatch = False) # Turn off option for setting missmatch
    arguments.args_common(parser)

    args = parser.parse_args()

    args.aligner = args.aligner.lower()
    if args.aligner not in ['hisat2', 'bowtie2', 'bwa']:
        print >> sys.stderr, "Error: --alinger supports hisat2, bowtie2, or bwa"
        exit(1)

    read_fnames, bam_fname = [], []
    if args.read_fname_U:
        read_fnames = [args.read_fname_U]
    elif args.read_fname_1 or args.read_fname_2:
        if not args.read_fname_1 or not args.read_fname_2:
            print >> sys.stderr, "Error: Please specify both -1 and -2 options"
            exit(1)
        read_fnames = [args.read_fname_1,
                       args.read_fname_2]
    elif args.alignment_fname:
        bam_fname = [args.alignment_fname]
    else:
        print >> sys.stderr, "Error: Please provide file options"
        exit(1)

    assert (read_fnames or bam_fname) and not (read_fnames and bam_fname)

    record_variants(args.aligner,
                    args.fastq,
                    args.reference,
                    args.graph_index,
                    read_fnames,
                    bam_fname,
                    False if bam_fname else True,
                    args.keep_alignment,
                    args.threads)

