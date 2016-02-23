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


import sys, os, subprocess
from argparse import ArgumentParser, FileType


"""
"""
def read_genome(genome_file):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    return chr_dic


"""
"""
def generate_haplotypes(snp_file,
                        haplotype_file,
                        vars,
                        genotypes_list,
                        inter_gap,
                        intra_gap,
                        num_haplotypes):
    assert len(vars) == len(genotypes_list)
    assert len(vars) > 0
    num_chromosomes = len(genotypes_list[0])

    # Write SNPs into a file (.snp)
    for var in vars:
        varID, chr, pos, type, data = var
        if type == 'S':
            type = "single"
        elif type == 'D':
            type = "deletion"
        else:
            assert type == 'I'
            type = "insertion"
        print >> snp_file, "%s\t%s\t%s\t%s\t%s" % \
            (varID, type, chr, pos, data)

    # genotypes_list looks like
    #    Var0: 000001000
    #    Var1: 010000000
    #    Var2: 001100000
    # Get haplotypes from genotypes_list
    haplotypes = set()
    cnv_genotypes = ["" for i in range(num_chromosomes)]
    for genotypes in genotypes_list:
        for i in range(len(genotypes)):
            genotype = genotypes[i]
            cnv_genotypes[i] += genotype

    cnv_genotypes = set(cnv_genotypes)
    for raw_haplotype in cnv_genotypes:
        if '1' not in raw_haplotype:
            continue
        haplotype = ""
        for i in range(len(raw_haplotype)):
            if raw_haplotype[i] == '1':
                if haplotype == "":
                    haplotype = str(i)
                else:
                    haplotype += ("#%d" % i)                    
        assert haplotype != ""            
        haplotypes.add(haplotype)
    # haplotypes look like
    #    '8#10#12#23', '8#12#23', '5#8#12#23#30'

    # Split some haplotypes that include large gaps inside
    def split_haplotypes(haplotypes):
        split_haplotypes = set()
        for haplotype in haplotypes:
            haplotype = haplotype.split('#')
            assert len(haplotype) > 0
            if len(haplotype) == 1:
                split_haplotypes.add(haplotype[0])
                continue
            prev_s, s = 0, 1
            while s < len(haplotype):
                _, _, prev_locus, prev_type, prev_data = vars[int(haplotype[s-1])]
                _, _, locus, type, data = vars[int(haplotype[s])]
                prev_locus, locus = int(prev_locus), int(locus)
                if prev_type == 'D':
                    prev_locus += (int(prev_data) - 1)
                if prev_locus + intra_gap < locus:
                    split_haplotypes.add('#'.join(haplotype[prev_s:s]))
                    prev_s = s
                s += 1
                if s == len(haplotype):
                    split_haplotypes.add('#'.join(haplotype[prev_s:s]))
        return split_haplotypes

    haplotypes2 = split_haplotypes(haplotypes)

    def cmp_haplotype(a, b):
        a = a.split('#')
        _, _, a1_locus, _, _ = vars[int(a[0])]
        _, _, a2_locus, a2_type, a2_data = vars[int(a[-1])]
        a_begin, a_end = int(a1_locus), int(a2_locus)
        if a2_type == 'D':
            a_end += (int(a2_data) - 1)
        b = b.split('#')
        _, _, b1_locus, _, _ = vars[int(b[0])]
        _, _, b2_locus, b2_type, b2_data = vars[int(b[-1])]
        b_begin, b_end = int(b1_locus), int(b2_locus)
        if b2_type == 'D':
            b_end += (int(b2_data) - 1)
        if a_begin != b_begin:
            return a_begin - b_begin
        return a_end - b_end
    
    haplotypes = sorted(list(haplotypes2), cmp=cmp_haplotype)

    # daehwan - for debugging purposes
    """
    dis = prev_locus - locus
    print "\n[%d, %d]: %d haplotypes" % (i, j, len(haplotypes)), dis
    if len(cur_vars) in range(0, 1000):
        # print "vars:", sorted(list(cur_vars), cmp=cmp_varKey
        print "num:", len(haplotypes)
        for haplotype in haplotypes:
            print haplotype.split('#')
        print "\nnum:", len(haplotypes2)
        for haplotype in haplotypes2:
            print haplotype.split('#')
    """

    # Write haplotypes
    for h_i in range(len(haplotypes)):
        h = haplotypes[h_i].split('#')
        _, chr, h1_locus, _, _ = vars[int(h[0])]
        _, _, h2_locus, h2_type, h2_data = vars[int(h[-1])]
        h_begin, h_end = int(h1_locus), int(h2_locus)
        if h2_type == 'D':
            h_end += (int(h2_data) - 1)
        assert h_begin <= h_end
        h_new_begin = h_begin
        for h_j in reversed(range(0, h_i)):
            hc = haplotypes[h_j].split('#')
            _, _, hc_begin, hc_type, hc_data = vars[int(hc[-1])]
            hc_begin = int(hc_begin)
            hc_end = hc_begin
            if hc_type == 'D':
                hc_end += (int(hc_data) - 1)
            if hc_end + inter_gap < h_begin:
                break
            if h_new_begin > hc_end:
                h_new_begin = hc_end
        assert h_new_begin <= h_begin
        h_add = [vars[int(id)][0] for id in h]
        print >> haplotype_file, "ht%d\t%s\t%d\t%d\t%s" % \
            (num_haplotypes, chr, h_new_begin, h_end, ','.join(h_add))
        num_haplotypes += 1

    return num_haplotypes


"""
"""
digit2str = [str(i) for i in range(10)]
def main(genome_file,
         base_fname,
         species,
         assembly_version,
         genotype_vcf,
         genotype_gene_list,
         inter_gap,
         intra_gap,
         verbose = False):
    # load genomic sequences
    chr_dic = read_genome(genome_file)
    
    VCF_fnames = {}
    if species == "human":
        grch38_url_base = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions"
        grch37_url_base = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
        if assembly_version == "GRCh38":
            VCF_fnames["human"] = [["%s" % chr, "%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz" % (grch38_url_base, chr)] for chr in [str(i+1) for i in range(22)] + ["X", "Y"]]
        else:
            assert assebmly_version == "GRCh37"
            VCF_fnames["human"] = [["%s" % chr, "%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" % (grch17_url_base, chr)] for chr in [str(i+1) for i in range(22)] + ["X", "Y"]]
        # MT is common in both two assemblies.
        VCF_fnames["human"].append(["MT", "%s/ALL.chrMT.phase3_callmom.20130502.genotypes.vcf.gz" % grch37_url_base])
    
    if species not in VCF_fnames:
        print >> sys.stderr, "Error: %s is not supported." % species
        sys.exit(1)

    genotype_chr_list, genotype_var_list = {}, {}
    if genotype_vcf != "":
        assert len(genotype_gene_list) > 0
        gzip_cmd = ["gzip", "-cd", genotype_vcf]
        gzip_proc = subprocess.Popen(gzip_cmd,
                                     stdout=subprocess.PIPE,
                                     stderr=open("/dev/null", 'w'))
        for line in gzip_proc.stdout:
            if line.startswith("#"):
                continue

            chr, pos, varID, ref_allele, alt_alleles, qual, filter, info = line.strip().split()
            pos = int(pos) - 1
            include = False
            for gene in genotype_gene_list:
                if info.find(gene) != -1:
                    include = True
                    break
            if include:
                if chr not in genotype_chr_list:
                    genotype_chr_list[chr] = [pos, pos]
                else:
                    if pos < genotype_chr_list[chr][0]:
                        genotype_chr_list[chr][0] = pos
                    elif pos > genotype_chr_list[chr][1]:
                        genotype_chr_list[chr][1] = pos
                genotype_var_list[varID] = False

        print genotype_chr_list
        print len(genotype_var_list)

    SNP_file = open("%s.snp" % base_fname, 'w')
    haplotype_file = open("%s.haplotype" % base_fname, 'w')

    num_haplotypes = 0
    num_unassigned = 0
    for chr, url in VCF_fnames[species]:
        if chr not in chr_dic:
            continue
        if len(genotype_gene_list) > 0 and \
                chr not in genotype_chr_list:
            continue
        
        chr_seq = chr_dic[chr]
        fname = url.split('/')[-1]

        if not os.path.exists(fname):
            os.system("wget %s" % (url))
        assert os.path.exists(fname)
        gzip_cmd = ["gzip", "-cd", fname]
        gzip_proc = subprocess.Popen(gzip_cmd,
                                     stdout=subprocess.PIPE,
                                     stderr=open("/dev/null", 'w'))
        genomeIDs = []
        vars, genotypes_list = [], []
        curr_right = -1
        prev_varID, prev_pos = "", -1
        num_lines = 0
        for line in gzip_proc.stdout:
            num_lines += 1
            if line.startswith("##"):
                continue

            fields = line.strip().split()
            chr, pos, varID, ref_allele, alt_alleles, qual, filter, info, format = fields[:9]
            genotypes = fields[9:]

            if line.startswith("#"):
                genomeIDs = genotypes
                continue

            assert len(genotypes) == len(genomeIDs)

            if not varID.startswith("rs") or \
                    ';' in varID:
                continue

            if varID == prev_varID or \
                    pos == prev_pos:
                continue

            pos = int(pos) - 1
            if num_lines % 10000 == 1:
                print >> sys.stderr, "\t%s:%d\r" % (chr, pos),
            
            if len(genotype_gene_list) > 0:
                assert chr in genotype_chr_list
                min_pos, max_pos = genotype_chr_list[chr]
                assert min_pos < max_pos
                if pos < min_pos:
                    continue
                elif pos > max_pos:
                    break

            if len(genotype_gene_list) > 0:
                if varID not in genotype_var_list:
                    # daehwan - for debugging purposes
                    # continue
                    None
                else:
                    genotype_var_list[varID] = True
                
            if curr_right + inter_gap < pos and len(vars) > 0:
                assert len(vars) == len(genotypes_list)
                num_haplotypes = generate_haplotypes(SNP_file,
                                                     haplotype_file,
                                                     vars,
                                                     genotypes_list,
                                                     inter_gap,
                                                     intra_gap,
                                                     num_haplotypes)
                vars, genotypes_list = [], []

            assert ',' not in ref_allele
            alt_alleles = alt_alleles.split(',')
            assert len(alt_alleles) < 10
            for a in range(len(alt_alleles)):
                alt_allele = alt_alleles[a]
                ref_allele2, pos2 = ref_allele, pos
                min_len = min(len(ref_allele2), len(alt_allele))
                assert min_len >= 1
                if min_len > 1:
                    ref_allele2 = ref_allele2[min_len - 1:]
                    alt_allele = alt_allele[min_len - 1:]
                    pos2 += (min_len - 1)

                type, data = '', ''
                if len(ref_allele2) == 1 and len(alt_allele) == 1:
                    type = 'S'
                    data = alt_allele
                    assert ref_allele2 != alt_allele
                    if chr_seq[pos2] != ref_allele2:
                        continue
                elif len(ref_allele2) == 1:
                    assert len(alt_allele) > 1
                    type = 'I'
                    data = alt_allele[1:]
                    if len(data) > 32:
                        continue
                    if chr_seq[pos] != ref_allele2:
                        continue
                elif len(alt_allele) == 1:
                    assert len(ref_allele2) > 1
                    type = 'D'
                    data = len(ref_allele2) - 1
                    if chr_seq[pos2:pos2+data+1] != ref_allele2:
                        continue
                else:
                    assert False

                cnv_genotypes = []
                for genotype in genotypes:
                    P1, P2 = genotype[0], genotype[2]
                    if P1 == digit2str[a + 1]:
                        cnv_genotypes.append('1')
                    else:
                        cnv_genotypes.append('0')
                    if P2 == digit2str[a + 1]:
                        cnv_genotypes.append('1')
                    else:
                        cnv_genotypes.append('0')

                # Skip SNPs not present in a given population (e.g. 2,504 genomes in 1000 Genomes Project)
                if '1' not in cnv_genotypes:
                    continue

                tmp_varID = varID
                if len(alt_alleles) > 1:
                    tmp_varID += (".%d" % a)
                vars.append([tmp_varID, chr, pos2, type, data])
                genotypes_list.append(''.join(cnv_genotypes))
                right = pos2
                if type == 'D':
                    right += (int(data) - 1)
                if curr_right < right:
                    curr_right = right
                    
            prev_varID = varID
            prev_pos = pos

        if len(vars) > 0:
            assert len(vars) == len(genotypes_list)
            num_haplotypes = generate_haplotypes(SNP_file,
                                                 haplotype_file,
                                                 vars,
                                                 genotypes_list,
                                                 inter_gap,
                                                 intra_gap,
                                                 num_haplotypes)
            vars, genotypes_list = [], []

    SNP_file.close()
    haplotype_file.close()


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract haplotypes from a VCF file')
    parser.add_argument('genome_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input genome file')
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="",
                        help="base filename for SNPs and haplotypes (default: same as --species)")
    parser.add_argument('--species',
                        dest='species',
                        type=str,
                        default="human",
                        help='species (default: human)')
    parser.add_argument('--assembly-version',
                        dest='assembly_version',
                        type=str,
                        default="GRCh38",
                        help='species (default: GRCh38)')
    parser.add_argument('--genotype-vcf',
                        dest='genotype_vcf',
                        type=str,
                        default="",
                        help='VCF file name for genotyping (default: empty)')
    parser.add_argument('--genotype-gene-list',
                        dest='genotype_gene_list',
                        type=str,
                        default="",
                        help='A comma-separated list of genes to be genotyped (default: empty)')
    parser.add_argument("--inter-gap",
                        dest="inter_gap",
                        type=int,
                        default=30,
                        help="Maximum distance for variants to be in the same haplotype")
    parser.add_argument("--intra-gap",
                        dest="intra_gap",
                        type=int,
                        default=50,
                        help="Break a haplotype into several haplotypes")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    if args.base_fname == "":
        args.base_fname = args.species

    if args.genotype_vcf != "":
        if args.genotype_gene_list == "":
            print >> sys.stderr, "Error: please specify --genotype-gene-list."
            sys.exit(1)
        args.genotype_gene_list = args.genotype_gene_list.split(',')
    else:
        args.genotype_gene_list = []

    main(args.genome_file,
         args.base_fname,
         args.species,
         args.assembly_version,
         args.genotype_vcf,
         args.genotype_gene_list,
         args.inter_gap,
         args.intra_gap,
         args.verbose)
