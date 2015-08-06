#!/usr/bin/env python

import sys, os
import random
use_message = '''
'''

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

    rna, mismatch, snp = True, True, True
    datasets = [
        ["genome", 1000000, "RNA", not mismatch, not snp],
        ["genome", 20000000, "RNA", not mismatch, not snp],
        ]

    for reference, numreads, rna, mismatch, snp in datasets:
        if mismatch:
            dirname = "%dM_%s_reads" % (numreads / 1000000, molecule_type)
        else:
            dirname = "%dM_%s_mismatch_reads" % (numreads / 1000000, molecule_type)
        if os.path.exists(dirname):
            continue
        os.mkdir(dirname)
        os.chdir(dirname)

        if reference == "genome":
            genome_fname = "../../data/genome.fa"
        elif reference == "":
            genome_fname = ""
        else:
            assert False

        if rna:
            gtf_fname = "../../data/genes.gtf"
        else:
            gtf_fname = "/dev/null"

        if snp:
            snp_fname = "../../data/genome.snp"
        else:
            snp_fname = "/dev/null"
            
        cmd = "../../aligners/bin/simulate_reads.py --sanity-check --num-fragment %d %s %s %s sim" % \
            (numreads, genome_fname, gtf_fname, snp_fname)
        print >> sys.stderr, cmd
        os.system(cmd)

        print >> sys.stderr, "shuffle reads sim_1.fa and sim_2.fa"
        shuffle_pairs("sim_1.fa", "sim_2.fa")
        shuffle_reads_cmd = " mv sim_1.fa.shuffle sim_1.fa"
        shuffle_reads_cmd += "; mv sim_2.fa.shuffle sim_2.fa"
        os.system(shuffle_reads_cmd)

        os.chdir("..")

    os.chdir("..")
            
    
if __name__ == "__main__":
    random.seed(0)
    simulate_reads()
