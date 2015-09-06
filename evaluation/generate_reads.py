#!/usr/bin/env python

import sys, os, random
from argparse import ArgumentParser, FileType

def shuffle_reads(read_fname, random_list):
    reads = []
    read_file = open(read_fname)
    for line in read_file:
        if line[0] == ">":
            reads.append([])
        reads[-1].append(line[:-1])        
    read_file.close()

    read_fname_out = read_fname + ".shuffle"
    read_file_out = open(read_fname_out, "w")
    assert len(random_list) == len(reads)
    for i in random_list:
        read = reads[random_list[i]]
        print >> read_file_out, "\n".join(read)
    read_file_out.close()


def shuffle_pairs(read1_fname, read2_fname):
    read1_file = open(read1_fname)
    num_reads = 0
    for line in read1_file:
        if line[0] == ">":
            num_reads += 1
    read1_file.close()

    random_list = [i for i in range(num_reads)]
    random.shuffle(random_list)

    shuffle_reads(read1_fname, random_list)
    shuffle_reads(read2_fname, random_list)


def simulate_reads():
    if not os.path.exists("reads"):
        os.mkdir("reads")
    os.chdir("reads")
    if not os.path.exists("simulation"):
        os.mkdir("simulation")
    os.chdir("simulation")

    _rna, _mismatch, _snp, _constant = True, True, True, True
    _dna = not _rna
    datasets = [
        ["22", 1000000, _dna, not _snp, not _mismatch, _constant],
        ["22", 1000000, _dna, not _snp, _mismatch, _constant],
        ["22", 1000000, _dna, _snp, not _mismatch, _constant],
        ["22", 1000000, _dna, _snp, _mismatch, _constant],
        ["22", 1000000, _rna, not _snp, not _mismatch, not _constant],
        ["22", 1000000, _rna, not _snp, not _mismatch, _constant],
        ["22", 1000000, _rna, not _snp, _mismatch, not _constant],
        ["22", 1000000, _rna, not _snp, _mismatch, _constant],
        ["22", 1000000, _rna, _snp, not _mismatch, not _constant],
        ["22", 1000000, _rna, _snp, not _mismatch, _constant],
        ["22", 1000000, _rna, _snp, _mismatch, not _constant],
        ["22", 1000000, _rna, _snp, _mismatch, _constant],
        ["22_20-21M", 1000000, _rna, not _snp, not _mismatch, not _constant],
        ["22_20-21M", 1000000, _rna, _snp, not _mismatch, _constant],
        ["genome", 20000000, _dna, not _snp, not _mismatch, _constant],
        ["genome", 20000000, _dna, _snp, not _mismatch, _constant],
        ["genome", 20000000, _rna, not _snp, not _mismatch, not _constant],
        ["genome", 20000000, _rna, _snp, not _mismatch, not _constant],
        ]

    data_dir_base = "../../../data"
    for genome, numreads, rna, snp, mismatch, constant in datasets:
        if rna:
            molecule = "RNA"
        else:
            molecule = "DNA"
        dirname = "%dM_%s" % (numreads / 1000000, molecule)
        if mismatch:
            dirname += "_mismatch"
        if snp:
            dirname += "_snp"
        if rna and constant:
            dirname += "_constant"
        dirname += "_reads"
        dirname += ("_" + genome)
        if os.path.exists(dirname):
            continue
        os.mkdir(dirname)
        os.chdir(dirname)
        genome_fname = data_dir_base + "/%s.fa" % (genome)

        if rna:
            gtf_fname = data_dir_base + "/%s.gtf" % (genome)
        else:
            gtf_fname = "/dev/null"

        if snp:
            snp_fname = data_dir_base + "/%s.snp" % (genome)
        else:
            snp_fname = "/dev/null"

        cmd_add = ""
        if not rna:
            cmd_add += "--dna "
        if mismatch:
            cmd_add += "--error-rate 0.5 "
        if rna and constant:
            cmd_add += "--expr-profile constant "
        cmd = "../../../aligners/bin/simulate_reads.py --sanity-check %s --num-fragment %d %s %s %s sim" % \
            (cmd_add, numreads, genome_fname, gtf_fname, snp_fname)
        print >> sys.stderr, cmd
        os.system(cmd)

        random.seed(0)
        print >> sys.stderr, "shuffle reads sim_1.fa and sim_2.fa"
        shuffle_pairs("sim_1.fa", "sim_2.fa")
        shuffle_reads_cmd = " mv sim_1.fa.shuffle sim_1.fa"
        shuffle_reads_cmd += "; mv sim_2.fa.shuffle sim_2.fa"
        os.system(shuffle_reads_cmd)

        os.chdir("..")

    os.chdir("..")
            
    
if __name__ == "__main__":
    parser = ArgumentParser(
        description='Generate reads using simulate_reads.py in HISAT2')
    args = parser.parse_args()
    simulate_reads()
