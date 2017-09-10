#!/usr/bin/env python
#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com> and Bongsoo Park <genomicspark@gmail.com>
#
# This file is part of HISAT-bisulfite.
#
# HISAT-bisulfite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT-bisulfite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT-bisulfite.  If not, see <http://www.gnu.org/licenses/>.
#


import os, sys, subprocess, re
import heapq
import inspect
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common


"""
Download BigWig (bw) files from Smith Lab at USC, and convert them to a HISAT2's snp file
   http://smithlab.usc.edu/methbase/data/ENCODE-Project-Human-2011
"""
def create_snp_from_bigwig(base_fname,
                           chr_dic,
                           abundance_cutoff,
                           sample_cutoff):
    if not os.path.exists("%s.cpg" % base_fname):
        def get_html(url):
            download_cmd = ["wget",
                            "-O", "-",
                            url]
            proc = subprocess.Popen(download_cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=open("/dev/null", 'w'))
            output = ""
            for line in proc.stdout:
                output += line
            return output

        # Refer to Python's regular expression at https://docs.python.org/2/library/re.html
        sample_re = re.compile('^<td>(Human_\S+)</td>$')

        url = "http://smithlab.usc.edu/methbase/data/ENCODE-Project-Human-2011"
        content = get_html(url).split("\n")
        content = map(lambda x: x.strip(), content)
        samples = []
        for line in content:
            sample_match = sample_re.search(line)
            if not sample_match:
                continue
            samples.append(sample_match.group(1))

        print >> sys.stderr, "Number of samples:", len(samples)
        for s in range(len(samples)):
            sample = samples[s]
            print >> sys.stderr, "\tDownloading %s ... (%d/%d)" % (sample, s+1, len(samples))
            sample_cmd = ["wget",
                          "http://smithlab.usc.edu/methbase/data/ENCODE-Project-Human-2011/%s/tracks_hg38/%s.meth.bw" % (sample, sample)]
            proc = subprocess.Popen(sample_cmd,
                                    stderr=open("/dev/null", 'w'))
            proc.communicate()

            print >> sys.stderr, "\tConverting %s.meth.bw to %s.meth.bedgraph ..." % (sample, sample)
            convert_cmd = ["bigWigToBedGraph",
                           "%s.meth.bw" % sample,
                           "%s.meth.bedgraph" % sample]
            proc = subprocess.Popen(convert_cmd)
            proc.communicate()
            os.system("rm -f %s.meth.bw" % sample)

        def read_line(file):
            line = file.readline()
            line = line.strip()
            if len(line) == 0:
                return None, None, None
            chr, pos, _, abundance = line.split('\t')
            chr, pos, abundance = chr[3:], int(pos), float(abundance)
            if chr == "M":
                chr = "MT"
            return chr, pos, abundance

        # Use heap (a.k.a. priority queue) to handle merge in a memory lean way
        sample_files, heap = {}, []
        for sample in samples:
            sample_fname = "%s.meth.bedgraph" % sample
            if not os.path.exists(sample_fname):
                continue
            sample_file = open(sample_fname)
            sample_files[sample] = sample_file
            chr, pos, abundance = read_line(sample_file)
            if chr != None:
                heapq.heappush(heap, [chr, pos, abundance, sample])

        print >> sys.stderr, "Merging %s samples ..." % len(heap)
        CpG_file = open("%s.cpg" % base_fname, 'w')
        prev_chr, prev_pos, prev_sample, cur_samples = None, None, None, []
        while len(heap) > 0:
            chr, pos, abundance, sample = heapq.heappop(heap)
            if prev_chr != None and (chr != prev_chr or pos != prev_pos):
                print >> CpG_file, "%s\t%s\t%s" % (prev_chr, prev_pos, ','.join(cur_samples))
                cur_samples = ["%s:%.2f" % (sample, abundance)]
            else:
                cur_samples.append("%s:%.2f" % (sample, abundance))
            prev_chr, prev_pos, prev_sample = chr, pos, sample
            chr, pos, abundance = read_line(sample_files[sample])
            if chr != None:
                heapq.heappush(heap, [chr, pos, abundance, sample])

        if len(cur_samples) > 0:
            print >> CpG_file, "%s\t%s\t%s" % (prev_chr, prev_pos, ','.join(cur_samples))
        CpG_file.close()
        
        for sample in sample_files.keys():
            os.system("rm -f %s.meth.bedgraph" % sample)

    # Convert CpG sites to snps and haplotypes in file formats used for HISAT2
    print >> sys.stderr, "Converting to HISAT2's snp and haplotype files ..."
    snp_file = open("%s.snp" % base_fname, 'w')
    haplotype_file = open("%s.haplotype" % base_fname, 'w')
    CpG_file = open("%s.cpg" % base_fname)
    cpg_num = 0
    miscpg_dic = {}
    for line in CpG_file:
        chr, pos, samples = line.strip().split()
        pos = int(pos)
        abundances = [float(sample.split(':')[1]) for sample in samples.split(',')]
        abundances = sorted(abundances, reverse = True)
        if abundances[min(len(abundances), sample_cutoff) - 1] < abundance_cutoff:
            continue
        assert chr in chr_dic
        chr_seq = chr_dic[chr]
        if chr_seq[pos:pos+2] != "CG":
            if chr not in miscpg_dic:
                miscpg_dic[chr] = 1
            else:
                miscpg_dic[chr] += 1
            continue
        print >> snp_file, "cpg%d\tsingle\t%s\t%d\tT" % (cpg_num, chr, pos)
        print >> haplotype_file, "ht%d\t%s\t%d\t%d\tcpg%d" % (cpg_num, chr, pos, pos + 1, cpg_num)
        cpg_num += 1
        print >> snp_file, "cpg%d\tsingle\t%s\t%d\tA" % (cpg_num, chr, pos + 1)
        print >> haplotype_file, "ht%d\t%s\t%d\t%d\tcpg%d" % (cpg_num, chr, pos, pos + 1, cpg_num)
        cpg_num += 1

    for chr, num in miscpg_dic.items():
        print >> sys.stderr, "Warning: %d non-CpG sites on chromosome %s" % (num, chr)
        
    snp_file.close()
    haplotype_file.close()
        

"""
"""
def build_epigenome(base_fname,                          
                    abundance_cutoff,
                    sample_cutoff,
                    region,
                    threads,
                    aligner,
                    graph_index,
                    verbose):
    # Make sure the following programs installed
    programs = [aligner, "%s-build" % aligner, "samtools", "bigWigToBedGraph"]
    for program in programs:
        def which(file):
            for path in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(path, file)):
                    return os.path.join(path, file)
            return None        
        if not which(program):
            print >> sys.stderr, "Error: %s is required, please have it installed on your system." % program
            sys.exit(1)
    
    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not typing_common.check_files(HISAT2_fnames):
        typing_common.download_genome_and_index()

    # Load genomic sequences
    chr_dic, chr_names, chr_full_names = typing_common.read_genome(open("genome.fa"))

    # Extract variants from dbSNP database
    CpG_fnames = ["%s.snp" % base_fname,
                  "%s.haplotype" % base_fname]
    if not typing_common.check_files(CpG_fnames):
        create_snp_from_bigwig(base_fname,
                               chr_dic,
                               abundance_cutoff,
                               sample_cutoff)

    genome_fname = "genome.fa"
    # Because building an index for the whole genome takes a lot of memory (~200 GB) and a couple of hours,
    # we provide a way to build and test a small index using --region option.
    if len(region) > 0:
        extract_seq_cmd = ["samtools", "faidx", "genome.fa"]
        if len(region) == 1:
            region_chr, region_left, region_right = region[0], 1, len(chr_dic[region_chr])
            region_str = region_str2 = region_chr
            extract_seq_cmd.append(region_chr)
        else:
            region_chr, region_left, region_right = region
            region_str = "%s_%d_%d" % (region_chr, region_left, region_right)
            region_str2 = "%s:%d-%d" % (region_chr, region_left, region_right)
            extract_seq_cmd.append("%s:%d-%d" % (region_chr, region_left, region_right))
        extract_seq_proc = subprocess.Popen(extract_seq_cmd,
                                            stdout=open("genome_%s.fa" % region_str, 'w'),
                                            stderr=open("/dev/null", 'w'))

        region_left, region_right = region_left - 1, region_right - 1
        region_snp_file = open("%s_%s.snp" % (base_fname, region_str), 'w')
        for line in open("%s.snp" % base_fname):
            id, type, chr, pos, data = line.strip().split()
            pos = int(pos)
            if chr != region_chr:
                continue
            if pos < region_left or pos > region_right:
                continue
            assert chr_dic[region_chr][pos:pos+2] == "CG" or chr_dic[region_chr][pos-1:pos+1] == "CG"
            print >> region_snp_file, "%s\t%s\t%s\t%d\t%s" % (id, type, region_str2, pos - region_left, data)
        region_snp_file.close()
        
        region_haplotype_file = open("%s_%s.haplotype" % (base_fname, region_str), 'w')
        for line in open("%s.haplotype" % base_fname):
            id, chr, left, right, ids = line.strip().split()
            left, right = int(left), int(right)
            if chr != region_chr:
                continue
            if left < region_left or right > region_right:
                continue
            print >> region_haplotype_file, "%s\t%s\t%d\t%d\t%s" % (id, region_str2, left - region_left, right - region_left, ids)
        region_haplotype_file.close()
        base_fname = "%s_%s" % (base_fname, region_str)
        genome_fname = "genome_%s.fa" % region_str

    # Build indexes based on the above information
    if graph_index:
        assert aligner == "hisat2"
        build_cmd = ["hisat2-build",
                     "-p", str(threads),
                     "--snp", "%s.snp" % base_fname,
                     "--haplotype", "%s.haplotype" % base_fname,
                     genome_fname,
                     "%s" % base_fname]
    else:        
        assert aligner in ["hisat2", "bowtie2"]
        build_cmd = ["%s-build" % aligner,
                     "-p" if aligner == "hisat2" else "--threads", str(threads),
                     genome_fname,
                     "%s" % base_fname]
    print >> sys.stderr, "Building HISAT2 index:", ' '.join(build_cmd)        
    subprocess.call(build_cmd,
                    stdout=open("/dev/null", 'w'),
                    stderr=open("/dev/null", 'w'))

    if aligner == "hisat2":
        index_fnames = ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]
    else:
        index_fnames = ["%s.%d.bt2" % (base_fname, i+1) for i in range(4)]
        index_fnames += ["%s.rev.%d.bt2" % (base_fname, i+1) for i in range(2)]
    if not typing_common.check_files(index_fnames):
        print >> sys.stderr, "Error: indexing failed."
        sys.exit(1)

        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Build epigenome")
    parser.add_argument("--base", "--base-fname",
                        dest="base_fname",
                        type=str,
                        default="epigenome",
                        help="base filename for genotype genome (default: epigenome)")
    parser.add_argument("--abundance-cutoff",
                        dest="abundance_cutoff",
                        type=float,
                        default=0.0,
                        help="Abundance cutoff")
    parser.add_argument("--sample-cutoff",
                        dest="sample_cutoff",
                        type=int,
                        default=1,
                        help="Sample cutoff")
    parser.add_argument("--region",
                        dest="region",
                        type=str,
                        default="",
                        help="Genomic region (1-offset) instead of the whole genome (default: (empty), example: 7 and 10:1000000-2000000)")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--aligner",
                        dest="aligner",
                        type=str,
                        default="hisat2",
                        help="Aligner (default: hisat2)")
    parser.add_argument("--linear-index",
                        dest="graph_index",
                        action="store_false",
                        help="Build linear index")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.aligner not in ["hisat2", "bowtie2"]:
        print >> sys.stderr, "Error: --aligner should be either hisat2 or bowtie2."
        sys.exit(1)

    region = []
    if args.region != "":
        try:
            region_chr, region_range = args.region.split(':')
            if region_range != "":
                region_left, region_right = region_range.split('-')
                region_left, region_right = int(region_left), int(region_right)
                region = [region_chr, region_left, region_right]
        except ValueError:
            print >> sys.stderr, "Error: --region %s is ill-formatted." % args.region
            sys.exit(1)
        
    build_epigenome(args.base_fname,
                    args.abundance_cutoff,
                    args.sample_cutoff,
                    region,
                    args.threads,
                    args.aligner,
                    args.graph_index,
                    args.verbose)
    
