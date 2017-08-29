#!/usr/bin/env python
#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT-genotype.
#
# HISAT-genotype is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT-genotype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT-genotype.  If not, see <http://www.gnu.org/licenses/>.
#


import sys, os, subprocess, re, resource
import inspect
import random
import glob
from argparse import ArgumentParser, FileType


"""
"""
if __name__ == '__main__':
    aligners = [["hisat2", "graph"], ["hisat2", "linear"], ["bowtie2", "linear"]]
    genes = ["A", "B", "C", "DQA1", "DQB1", "DRB1"]

    samples = ["NA12878", "LP6005041-DNA_A01", "LP6005045-DNA_D07"]
    for sample in samples:
        for aligner, type in aligners:
            sample_dir = "%s_%s" % (sample, aligner)
            if aligner == "hisat2":
                sample_dir += (".%s" % type)
            if not os.path.exists(sample_dir):
                continue

            fq_fnames = glob.glob("%s/*.fq.gz" % sample_dir)
            assert len(fq_fnames) == 2

            regions, region_loci, region_count, region_read1_count, region_read2_count = {}, {}, {}, {}, {}
            for line in open("%s.locus" % ("genotype_genome" if type == "graph" else "genotype_genome.linear")):
                family, allele_name, chr, left, right = line.strip().split()[:5]

                # DK - debugging purposes
                if family != "HLA":
                    continue
                
                region_name = "%s-%s" % (family, allele_name.split('*')[0])
                assert region_name not in regions
                regions[region_name] = allele_name
                left, right = int(left), int(right)
                if chr not in region_loci:
                    region_loci[chr] = {}
                region_loci[chr][region_name] = [allele_name, chr, left, right]
    

            aligner_cmd = [aligner]
            aligner_cmd += ["-x", "genotype_genome" if type == "graph" else "genotype_genome.linear"]
            if aligner == "hisat2":
                aligner_cmd += ["--no-spliced-alignment"]
            aligner_cmd += ["-X", "1000"]
            aligner_cmd += ["-1", fq_fnames[0],
                            "-2", fq_fnames[1]]
            # print >> sys.stderr, "Running:", ' '.join(aligner_cmd)
            print sample, aligner, type
            align_proc = subprocess.Popen(aligner_cmd,
                                          stdout=subprocess.PIPE,
                                          stderr=open("/dev/null", 'w'))

            prev_read_name, extract, read1_extract, read2_extract, read1_first, read2_first = "", set(), set(), set(), True, True
            for line in align_proc.stdout:
                if line.startswith('@'):
                    continue
                line = line.strip()
                cols = line.split()
                read_name, flag, chr, pos, mapQ, cigar, _, _, _, read, qual = cols[:11]
                flag, pos = int(flag), int(pos) - 1
                strand = '-' if flag & 0x10 else '+'                   
                AS, XS, NH = "", "", ""
                for i in range(11, len(cols)):
                    col = cols[i]
                    if col.startswith("AS"):
                        AS = int(col[5:])
                    elif col.startswith("XS"):
                        XS = int(col[5:])
                    elif col.startswith("NH"):
                        NH = int(col[5:])

                if read_name != prev_read_name:
                    for region in extract:
                        if region not in region_count:
                            region_count[region] = 0
                        region_count[region] += 1

                    for region in read1_extract:
                        if region not in region_read1_count:
                            region_read1_count[region] = 0
                        region_read1_count[region] += 1

                    for region in read2_extract:
                        if region not in region_read2_count:
                            region_read2_count[region] = 0
                        region_read2_count[region] += 1

                    prev_read_name, extract, read1_extract, read2_extract, read1_first, read2_first = "", set(), set(), set(), True, True

                if ((aligner == "hisat2" and NH == 1) or (aligner == "bowtie2" and AS > XS and read1_first if flag & 0x40 else read2_first)):
                    if chr in region_loci:
                        for region, loci in region_loci[chr].items():
                            _, _, loci_left, loci_right = loci
                            # there might be a different candidate region for each of left and right reads
                            if pos >= loci_left and pos < loci_right:
                                extract.add(region)
                                if flag & 0x40:
                                    read1_extract.add(region)
                                else:
                                    read2_extract.add(region)
                                break

                if flag & 0x40: # left read
                    read1_first = False
                else:
                    assert flag & 0x80 # right read
                    read2_first = False

                prev_read_name = read_name

            for gene in genes:
                gene = "HLA-" + gene
                if gene not in region_count:
                    continue
                print "\t%s pair: %d, left+right: %d" % (gene, region_count[gene], region_read1_count[gene] + region_read2_count[gene])
            
