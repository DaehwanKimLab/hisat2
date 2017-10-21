#!/usr/bin/env python

import sys, os, subprocess
import multiprocessing
import string, re
import platform
from datetime import datetime, date, time
import copy
from argparse import ArgumentParser, FileType


"""
"""
def evaluate(read_fname,
             verbose):
    aligners = [
        ["hisat2", "", "", ""],
        ["hisat2", "", "snp", ""],
        ["bowtie2", "", "", ""],
        ]
    num_threads = 3

    cwd = os.getcwd()
    genome = "genome"
    align_stat = []
    # for paired in [False, True]:
    for paired in [False]:
        base_fname = "common_snp_reads"
        type_sam_fname = base_fname + ".sam"
        type_read1_fname = base_fname +  "_1.fa"
        type_read2_fname = base_fname +  "_2.fa"

        type_read1_fname = read_fname

        def get_aligner_version(aligner, version):
            version = ""
            if aligner == "hisat2" or \
                    aligner == "bowtie2":
                if version:
                    cmd = ["%s_%s/%s" % (aligner, version, aligner)]
                else:
                    cmd = ["%s/%s" % (aligner_bin_base, aligner)]
                cmd += ["--version"]                    
                cmd_process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                version = cmd_process.communicate()[0][:-1].split("\n")[0]
                version = version.split()[-1]
            elif aligner == "star":
                version = "2.4.2a"
            elif aligner == "gsnap":
                cmd = ["%s/gsnap" % (aligner_bin_base)]
                cmd_process = subprocess.Popen(cmd, stderr=subprocess.PIPE)
                version = cmd_process.communicate()[1][:-1].split("\n")[0]
                version = version.split()[2]
            elif aligner == "bwa":
                cmd = ["%s/bwa" % (aligner_bin_base)]
                cmd_process = subprocess.Popen(cmd, stderr=subprocess.PIPE)
                version = cmd_process.communicate()[1][:-1].split("\n")[2]
                version = version.split()[1]

            return version

        def get_aligner_cmd(aligner, type, index_type, version, read1_fname, read2_fname, out_fname, cmd_idx = 0):
            cmd = []
            if aligner == "hisat2":
                cmd = ["hisat2"]
                if num_threads > 1:
                    cmd += ["-p", str(num_threads)]
                cmd += ["-f"]
                cmd += ["--no-spliced-alignment"]
                if index_type:
                    index_cmd = "../grch38_snp_hisat2/genome_snp"
                else:
                    index_cmd = "../grch38_hisat2/genome"
                cmd += [index_cmd]
                if paired:
                    cmd += ["-1", read1_fname,
                            "-2", read2_fname]
                else:
                    cmd += [read1_fname]                        
            elif aligner == "star":
                cmd = ["%s/STAR" % (aligner_bin_base)]
                if num_threads > 1:
                    cmd += ["--runThreadN", str(num_threads)]
                cmd += ["--genomeDir"]
                if cmd_idx == 0:
                    if type == "gtf":
                        cmd += ["%s/STAR%s/gtf" % (index_base, index_add)]
                    else:
                        cmd += ["%s/STAR%s" % (index_base, index_add)]
                else:
                    assert cmd_idx == 1
                    cmd += ["."]

                if desktop:
                    cmd += ["--genomeLoad", "NoSharedMemory"]
                else:
                    cmd += ["--genomeLoad", "LoadAndKeep"]
                if type == "x2":
                    if cmd_idx == 1:
                        cmd += ["--alignSJDBoverhangMin", "1"]
                cmd += ["--readFilesIn",
                        read1_fname]
                if paired:
                    cmd += [read2_fname]
                if paired:
                    cmd += ["--outFilterMismatchNmax", "6"]
                else:
                    cmd += ["--outFilterMismatchNmax", "3"]
            elif aligner == "bowtie2":
                cmd = ["bowtie2"]
                if num_threads > 1:
                    cmd += ["-p", str(num_threads)]
                cmd += ["-f"]
                cmd += ["-x ../grch38_bowtie2/genome"]
                if paired:
                    cmd += ["-1", read1_fname,
                            "-2", read2_fname]
                else:
                    cmd += [read1_fname]
            elif aligner == "gsnap":
                cmd = ["%s/gsnap" % (aligner_bin_base),
                       "-A",
                       "sam"]
                if num_threads > 1:
                    cmd += ["-t", str(num_threads)]
                cmd += ["--max-mismatches=3",
                        "-D", "%s/GSNAP%s" % (index_base, index_add),
                        "-N", "1",
                        "-d", genome,
                        read1_fname]
                if paired:
                    cmd += [read2_fname]
            elif aligner == "bwa":
                cmd = ["%s/bwa" % (aligner_bin_base)]
                if type in ["mem", "aln"]:
                    cmd += [type]
                elif type == "sw":
                    cmd += ["bwa" + type]
                if num_threads > 1:
                    cmd += ["-t", str(num_threads)]
                cmd += ["%s/BWA%s/%s.fa" % (index_base, index_add, genome)]
                cmd += [read1_fname]
                if paired:
                    cmd += [read2_fname]
            else:
                assert False

            return cmd

        for aligner, type, index_type, version in aligners:
            aligner_name = aligner + type
            if version != "":
                aligner_name += ("_%s" % version)
            if aligner == "hisat2" and index_type != "":
                aligner_name += ("_" + index_type)
            two_step = (aligner == "tophat2" or type == "x2" or (aligner in ["hisat2", "hisat"] and type == ""))
            print >> sys.stderr, "\t%s\t%s" % (aligner_name, str(datetime.now()))
            if paired:
                aligner_dir = aligner_name + "_paired"
            else:
                aligner_dir = aligner_name + "_single"
            if not os.path.exists(aligner_dir):
                os.mkdir(aligner_dir)
            os.chdir(aligner_dir)

            out_fname = base_fname + ".sam"
            duration = -1.0

            # Align all reads
            aligner_cmd = get_aligner_cmd(aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname)
            start_time = datetime.now()
            if verbose:
                print >> sys.stderr, "\t", start_time, " ".join(aligner_cmd)
            if aligner in ["hisat2", "hisat", "bowtie", "bowtie2", "gsnap", "bwa"]:
                proc = subprocess.Popen(aligner_cmd, stdout=open(out_fname, "w"), stderr=subprocess.PIPE)
            else:
                proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc.communicate()
            finish_time = datetime.now()
            duration = finish_time - start_time
            duration = duration.total_seconds()
            if verbose:
                print >> sys.stderr, "\t", finish_time, "finished:", duration

            assert os.path.exists(out_fname)
            correct_reads, correct_multi_reads, num_reads = 0, 0, 0
            prev_read_id = None
            for line in open(out_fname):
                if line.startswith('@'):
                    continue
                read_id, flag, chr, pos, mapQ, cigar = line.split()[:6]
                if chr.startswith("chr"):
                    chr = chr[3:]
                pos = int(pos) - 1
                true_chr, true_pos, true_cigar = read_id.split('_')[1:4]
                true_pos = int(true_pos)

                if read_id != prev_read_id:
                    num_reads += 1

                if true_chr == chr and pos == true_pos and cigar == true_cigar:
                    correct_multi_reads += 1
                    if prev_read_id != read_id:
                        correct_reads += 1

                prev_read_id = read_id

            print >> sys.stderr, "\tfirst: %d / %d (%.2f%%)" % (correct_reads, num_reads, float(correct_reads)/num_reads*100)
            print >> sys.stderr, "\tall: %d / %d (%.2f%%)" % (correct_multi_reads, num_reads, float(correct_multi_reads)/num_reads*100)

            os.chdir("..")


"""
"""
if __name__ == "__main__":
    parser = ArgumentParser(
        description='test HISAT2, and compare HISAT2 with other popular aligners such as TopHat2, STAR, Bowtie1/2, GSNAP, BWA-mem, etc.')
    parser.add_argument('read_fname',
                        nargs='?',
                        type=str,
                        help='input read file')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    evaluate(args.read_fname,
             args.verbose)

    
