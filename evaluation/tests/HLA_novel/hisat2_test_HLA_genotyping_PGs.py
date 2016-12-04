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


import sys, os, subprocess, re
import inspect
import random
from argparse import ArgumentParser, FileType

# Gold Standard (experimentally verified, a lot of literature, ...)
gold_allele_info = {
    "NA12877" : {"A" : ["03:01", "02:01"], "B" : ["15:01", "44:02"], "C" : ["05:01", "03:04"], "DQA1" : ["03:03", "03:01"], "DQB1" : ["03:02", "03:01"], "DRB1" : ["04:03", "04:03"]},
    "NA12878" : {"A" : ["01:01", "11:01"], "B" : ["08:01", "56:01"], "C" : ["01:02", "07:01"], "DQA1" : ["05:01", "01:01"], "DQB1" : ["02:01", "05:01"], "DRB1" : ["03:01", "01:01"]},
    "NA12879" : {"A" : ["01:01", "02:01"], "B" : ["08:01", "15:01"], "C" : ["03:04", "07:01"], "DQA1" : ["03:01", "05:01"], "DQB1" : ["03:02", "02:01"], "DRB1" : ["03:01", "04:03"]},
    "NA12880" : {"A" : ["02:01", "01:01"], "B" : ["15:01", "08:01"], "C" : ["03:04", "07:01"], "DQA1" : ["03:01", "05:01"], "DQB1" : ["03:02", "02:01"], "DRB1" : ["03:01", "04:03"]},
    "NA12881" : {"A" : ["03:01", "11:01"], "B" : ["44:02", "56:01"], "C" : ["05:01", "01:02"], "DQA1" : ["03:03", "01:01"], "DQB1" : ["03:01", "05:01"], "DRB1" : ["04:03", "01:01"]},
    "NA12882" : {"A" : ["02:01", "11:01"], "B" : ["15:01", "56:01"], "C" : ["01:02", "03:04"], "DQA1" : ["03:01", "01:01"], "DQB1" : ["03:02", "05:01"], "DRB1" : ["04:03", "01:01"]},
    "NA12883" : {"A" : ["03:01", "11:01"], "B" : ["44:02", "56:01"], "C" : ["01:02", "05:01"], "DQA1" : ["03:03", "01:01"], "DQB1" : ["03:01", "05:01"], "DRB1" : ["01:01", "04:03"]},
    "NA12884" : {"A" : ["02:01", "11:01"], "B" : ["15:01", "56:01"], "C" : ["01:02", "03:04"], "DQA1" : ["03:01", "01:01"], "DQB1" : ["03:02", "05:01"], "DRB1" : ["01:01", "04:03"]},
    "NA12885" : {"A" : ["03:01", "01:01"], "B" : ["44:02", "08:01"], "C" : ["05:01", "07:01"], "DQA1" : ["03:03", "05:01"], "DQB1" : ["03:01", "02:01"], "DRB1" : ["03:01", "04:03"]},
    "NA12886" : {"A" : ["03:01", "01:01"], "B" : ["44:02", "08:01"], "C" : ["07:01", "05:01"], "DQA1" : ["03:03", "05:01"], "DQB1" : ["02:01", "03:01"], "DRB1" : ["03:01", "04:03"]},
    "NA12887" : {"A" : ["02:01", "01:01"], "B" : ["15:01", "08:01"], "C" : ["03:04", "07:01"], "DQA1" : ["03:01", "05:01"], "DQB1" : ["03:02", "02:01"], "DRB1" : ["03:01", "04:03"]},
    "NA12888" : {"A" : ["01:01", "02:01"], "B" : ["08:01", "15:01"], "C" : ["07:01", "03:04"], "DQA1" : ["03:01", "05:01"], "DQB1" : ["03:02", "02:01"], "DRB1" : ["03:01", "04:03"]},
    "NA12889" : {"A" : ["03:01", "03:01"], "B" : ["07:02", "44:02"], "C" : ["05:01", "07:02"], "DQA1" : ["03:03", "01:02"], "DQB1" : ["03:01", "06:02"], "DRB1" : ["15:01", "04:03"]},
    "NA12890" : {"A" : ["03:01", "02:01"], "B" : ["44:03", "15:01"], "C" : ["16:01", "03:04"], "DQA1" : ["03:01", "02:01"], "DQB1" : ["03:02", "02:02"], "DRB1" : ["04:03", "07:01"]},
    "NA12891" : {"A" : ["24:02", "01:01"], "B" : ["08:01", "07:02"], "C" : ["07:02", "07:01"], "DQA1" : ["05:01", "01:02"], "DQB1" : ["06:02", "02:01"], "DRB1" : ["03:01", "15:01"]},
    "NA12892" : {"A" : ["02:01", "11:01"], "B" : ["15:01", "56:01"], "C" : ["01:02", "04:01"], "DQA1" : ["01:01", "01:01"], "DQB1" : ["05:01", "05:01"], "DRB1" : ["01:01", "01:02"]},
    "NA12893" : {"A" : ["02:01", "11:01"], "B" : ["15:01", "56:01"], "C" : ["01:02", "03:04"], "DQA1" : ["03:01", "01:01"], "DQB1" : ["03:02", "05:01"], "DRB1" : ["01:01", "04:03"]}
    }

# CEPH pedigree (17 family members)
pedigree = {
    "NA12889" : {"gender" : "M", "spouse" : "NA12890", "children" : ["NA12877"]},
    "NA12890" : {"gender" : "F", "spouse" : "NA12889", "children" : ["NA12877"]},
    "NA12877" : {"gender" : "M", "father" : "NA12889", "mother" : "NA12890", "spouse" : "NA12878", "children" : ["NA12879", "NA12880", "NA12881", "NA12882", "NA12883", "NA12884", "NA12885", "NA12886", "NA12887", "NA12888", "NA12893"]},

    "NA12891" : {"gender" : "M", "spouse" : "NA12892", "children" : ["NA12878"]},
    "NA12892" : {"gender" : "F", "spouse" : "NA12891", "children" : ["NA12878"]},
    "NA12878" : {"gender" : "F", "father" : "NA12892", "mother" : "NA12891", "spouse" : "NA12877", "children" : ["NA12879", "NA12880", "NA12881", "NA12882", "NA12883", "NA12884", "NA12885", "NA12886", "NA12887", "NA12888", "NA12893"]},

    "NA12879" : {"gender" : "F", "father" : "NA12877", "mother" : "NA12878"},
    "NA12880" : {"gender" : "F", "father" : "NA12877", "mother" : "NA12878"},
    "NA12881" : {"gender" : "F", "father" : "NA12877", "mother" : "NA12878"},
    "NA12882" : {"gender" : "M", "father" : "NA12877", "mother" : "NA12878"},
    "NA12883" : {"gender" : "M", "father" : "NA12877", "mother" : "NA12878"},
    "NA12884" : {"gender" : "M", "father" : "NA12877", "mother" : "NA12878"},
    "NA12885" : {"gender" : "F", "father" : "NA12877", "mother" : "NA12878"},
    "NA12886" : {"gender" : "M", "father" : "NA12877", "mother" : "NA12878"},
    "NA12887" : {"gender" : "F", "father" : "NA12877", "mother" : "NA12878"},
    "NA12888" : {"gender" : "M", "father" : "NA12877", "mother" : "NA12878"},
    "NA12893" : {"gender" : "M", "father" : "NA12877", "mother" : "NA12878"},
    }


"""
"""
def test_HLA_genotyping(reference_type,
                        hla_list,
                        aligners,
                        query_genomes,
                        exclude_allele_list,
                        num_mismatch,
                        verbose):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)

    if not os.path.exists("illumina/HLA"):
        print >> sys.stderr, "Error: illumina/HLA data is needed (please send an email to infphilo@gmail.com for getting the data)"
        sys.exit(1)

    num_test, num_success = 0, 0
    for genome in sorted(gold_allele_info.keys()):
        if not genome in query_genomes:
            continue
        genes = gold_allele_info[genome]
        read_fname_1, read_fname_2 = "illumina/HLA/%s.fished_1.fq" % genome, "illumina/HLA/%s.fished_2.fq" % genome
        if not os.path.exists(read_fname_1) or not os.path.exists(read_fname_2):
            continue
        print >> sys.stderr, genome        
        cmd_aligners = ['.'.join(aligners[i]) for i in range(len(aligners))]
        test_hla_script = os.path.join(ex_path, "hisat2_test_HLA_genotyping.py")
        for gene in sorted(genes.keys()):
            if not gene in hla_list:
                continue
            alleles = genes[gene]
            print >> sys.stderr, "\t%s - %s" % (gene, ' / '.join(alleles))            
            test_hla_cmd = [test_hla_script,
                            "--reference-type", reference_type,
                            "--hla-list", gene,
                            "--aligner-list", ','.join(cmd_aligners),
                            "--reads", "%s,%s" % (read_fname_1, read_fname_2),
                            "--best-alleles",
                            "--exclude-allele-list", ','.join(exclude_allele_list),
                            "--num-mismatch", str(num_mismatch)]

            if verbose:
                print >> sys.stderr, ' '.join(test_hla_cmd)
            
            proc = subprocess.Popen(test_hla_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
            num_test += 2
            test_alleles = set()
            for line in proc.stdout:
                print "\t\t", line,
                model, allele = line.split()[:2]
                if model != "SingleModel":
                    continue
                allele = allele.split('*')[1]
                allele = ':'.join(allele.split(':')[:2])
                test_alleles.add(allele)
            proc.communicate()
            for allele in alleles:
                if allele in test_alleles:
                    num_success += 1

    print >> sys.stderr, "%d/%d (%.2f%%)" % (num_success, num_test, num_success * 100.0 / num_test)


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping for Platinum Genomes')
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--hla-list",
                        dest="hla_list",
                        type=str,
                        default="A,B,C,DQA1,DQB1,DRB1",
                        help="A comma-separated list of HLA genes (default: A,B,C,DQA1,DQB1,DRB1)")
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="hisat2.graph",
                        help="A comma-separated list of aligners (default: hisat2.graph)")
    genomes_default = ','.join(gold_allele_info.keys())
    parser.add_argument("--genome-list",
                        dest="genome_list",
                        type=str,
                        default=genomes_default,
                        help="A comma-separated list of genomes (default: %s)" % genomes_default)
    parser.add_argument("--exclude-allele-list",
                        dest="exclude_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of allleles to be excluded")
    parser.add_argument("--num-mismatch",
                        dest="num_mismatch",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()

    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    args.hla_list = args.hla_list.split(',')
    if args.aligners == "":
        print >> sys.stderr, "Error: --aligners must be non-empty."
        sys.exit(1)    
    args.aligners = args.aligners.split(',')
    for i in range(len(args.aligners)):
        args.aligners[i] = args.aligners[i].split('.')
    args.genome_list = args.genome_list.split(',')
    args.exclude_allele_list = args.exclude_allele_list.split(',')

    test_HLA_genotyping(args.reference_type,
                        args.hla_list,
                        args.aligners,
                        args.genome_list,
                        args.exclude_allele_list,
                        args.num_mismatch,
                        args.verbose)
