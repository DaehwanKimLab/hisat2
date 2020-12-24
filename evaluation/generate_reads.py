#!/usr/bin/env python2

import sys, os, random
from argparse import ArgumentParser, FileType
from multiprocessing import Process

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


def simulate_reads(input_genomes,
                   input_molecules,
                   numflag_list,
                   fresh,
                   v1):
    if not os.path.exists("reads"):
        os.mkdir("reads")
    os.chdir("reads")
    if not os.path.exists("simulation"):
        os.mkdir("simulation")
    os.chdir("simulation")

    if input_genomes:
        genomes = input_genomes
    else:
        genomes = ["22_20-21M", "22", "genome"]

    if input_molecules:
        molecules = input_molecules
    else:
        molecules = ["DNA", "RNA"]

    _mismatch, _snp, _constant, _allele = True, True, True, True
    datasets = [
        ["22",     10000,    "RNA", not _snp, not _mismatch,     _constant, not _allele],
        ["22",     10000,    "RNA",     _snp, not _mismatch,     _constant,     _allele],
        ["22",     1000000,  "DNA", not _snp, not _mismatch,     _constant, not _allele],
        ["22",     1000000,  "DNA", not _snp,     _mismatch,     _constant, not _allele],
        ["22",     1000000,  "DNA",     _snp, not _mismatch,     _constant, not _allele],
        ["22",     1000000,  "DNA",     _snp,     _mismatch,     _constant, not _allele],
        ["22",     1000000,  "RNA", not _snp, not _mismatch, not _constant, not _allele],
        ["22",     1000000,  "RNA", not _snp, not _mismatch,     _constant, not _allele],
        ["22",     1000000,  "RNA", not _snp,     _mismatch, not _constant, not _allele],
        ["22",     1000000,  "RNA", not _snp,     _mismatch,     _constant, not _allele],
        ["22",     1000000,  "RNA",     _snp, not _mismatch, not _constant, not _allele],
        ["22",     1000000,  "RNA",     _snp, not _mismatch,     _constant, not _allele],
        ["22",     1000000,  "RNA",     _snp,     _mismatch, not _constant, not _allele],
        ["22",     1000000,  "RNA",     _snp,     _mismatch,     _constant, not _allele],
        ["genome", 10000000, "DNA", not _snp, not _mismatch,     _constant, not _allele],
        ["genome", 10000000, "DNA",     _snp, not _mismatch,     _constant, not _allele],
        ["genome", 10000000, "DNA",     _snp,     _mismatch,     _constant, not _allele],
        ["genome", 10000000, "DNA", not _snp, not _mismatch, not _constant, not _allele],
        ["genome", 10000000, "RNA",     _snp, not _mismatch, not _constant, not _allele],
        ["genome", 10000000, "RNA",     _snp,     _mismatch, not _constant, not _allele],
    ]

    data_dir_base = "../../../data"

    def generate_reads(cmd):
        print >> sys.stderr, cmd
        os.system(cmd)

        random.seed(0)
        print >> sys.stderr, "shuffle reads sim_1.fa and sim_2.fa"
        shuffle_pairs("sim_1.fa", "sim_2.fa")
        shuffle_reads_cmd = " mv sim_1.fa.shuffle sim_1.fa"
        shuffle_reads_cmd += "; mv sim_2.fa.shuffle sim_2.fa"
        os.system(shuffle_reads_cmd)

    pid_list = []

    for genome, numfrags, molecule, snp, mismatch, constant, allele in datasets:
        if genome not in genomes:
            continue
        
        if molecule not in molecules:
            continue

        if numfrag_list and numfrags not in numfrag_list:
            continue
            
        if numfrags >= 1000000:
            dirname = "%dM_%s" % (numfrags / 1000000, molecule)
        else:
            dirname = "%dK_%s" % (numfrags / 1000, molecule)

        if mismatch:
            dirname += "_mismatch"
        if snp:
            dirname += "_snp"
        if allele:
            dirname += "_allele"
        if molecule == "RNA" and constant:
            dirname += "_constant"
        dirname += "_reads"
        dirname += ("_" + genome)

        if os.path.exists(dirname):
            if fresh:
                os.system("rm -rf %s" % (dirname))
            else:
                continue
        
        os.mkdir(dirname)
        os.chdir(dirname)
        genome_fname = data_dir_base + "/%s.fa" % (genome)

        if molecule == "RNA":
            gtf_fname = data_dir_base + "/%s.gtf" % (genome)
        else:
            gtf_fname = "/dev/null"

        if snp:
            snp_fname = data_dir_base + "/%s.snp" % (genome)
        else:
            snp_fname = "/dev/null"

        cmd_add = ""
        if molecule == "DNA":
            cmd_add += "--dna "
        if mismatch:
            cmd_add += "--error-rate 0.2 "
        if molecule == "RNA" and constant:
            cmd_add += "--expr-profile constant "

        if allele:
            cmd_add += "--allele-specific "
            
            # DK - debugging purposes
            cmd_add += "--num-tran 10 "

        if v1:
            simulator = "hisat2_simulate_reads_v1.py"
        else:
            simulator = "hisat2_simulate_reads.py"
        cmd = "../../../aligners/bin/%s --sanity-check %s --num-fragment %d %s %s %s sim" % \
            (simulator, cmd_add, numfrags, genome_fname, gtf_fname, snp_fname)

        p = Process(target=generate_reads, args=(cmd,))
        p.start()
        pid_list.append(p)

        os.chdir("..")

    os.chdir("..")

    # wait
    for p in pid_list:
        p.join()
            
    
if __name__ == "__main__":
    parser = ArgumentParser(
        description='Generate reads using simulate_reads.py in HISAT2')
    parser.add_argument('--genome-list',
                        dest='genome_list',
                        type=str,
                        default="",
                        help='comma-separated list of genomes (e.g. 22,genome')
    parser.add_argument('--molecule-list',
                        dest='molecule_list',
                        type=str,
                        default="",
                        help='comma-separated list of molecule types (e.g. DNA,RNA')
    parser.add_argument('--numfrag-list',
                        dest='numfrag_list',
                        type=str,
                        default="",
                        help='comma-separated list of fragment numbers')
    parser.add_argument('--fresh',
                        dest='fresh',
                        action='store_true',
                        help='delete existing directories')
    parser.add_argument('--v1',
                        dest='v1',
                        action='store_true',
                        help='use an original version of hisat2_simulate_reads_v1.py')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()

    genomes = []
    for genome in args.genome_list.split(','):
        if genome == "":
            continue
        genomes.append(genome)

    molecules = []
    for molecule in args.molecule_list.split(','):
        if molecule == "":
            continue
        molecules.append(molecule)

    numfrag_list = []
    for frags in args.numfrag_list.split(','):
        numfrag_list.append(int(frags))
    
    simulate_reads(genomes,
                   molecules,
                   numfrag_list,
                   args.fresh,
                   args.v1)

