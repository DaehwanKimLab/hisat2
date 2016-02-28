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

digit2str = [str(i) for i in range(10)]

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
Compare two variants [chr, pos, type, data, dic]
"""
def compare_vars(a, b):
    a_chr, a_pos, a_type, a_data = a[:4]
    b_chr, b_pos, b_type, b_data = b[:4]
    assert a_chr == b_chr
    if a_pos != b_pos:
        return a_pos - b_pos
    if a_type != b_type:
         if a_type == 'I':
             return -1
         elif b_type == 'I':
             return 1
         if a_type == 'S':
             return -1
         else:
             return 1
    if a_data < b_data:
        return -1
    elif a_data > b_data:
        return 1
    else:
        return 0


"""
"""
def compatible_vars(a, b):
    a_chr, a_pos, a_type, a_data = a[:4]
    b_chr, b_pos, b_type, b_data = b[:4]
    assert a_chr == b_chr
    assert a_pos <= b_pos
    if a_pos == b_pos:
        return False
    if a_type == 'D':
        if b_pos <= a_pos + a_data:
            return False
    return True


"""
Given a VCF line, the function reports a list of variants [pos, type, data]
type: 'S' for single nucleotide polymorphism, 'D' for deletion, and 'I' for insertion
"""
def extract_vars(chr_dic, chr, pos, ref_allele, alt_alleles, varID):
    chr_seq = chr_dic[chr]
    vars = []
    assert ',' not in ref_allele
    alt_alleles = alt_alleles.split(',')    
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
        varID2 = varID
        if len(alt_alleles) > 1:
            varID2 = "%s.%d" % (varID, a)
        vars.append([chr, pos2, type, data, {"id":varID, "id2":varID2}])
                    
    return vars


"""
"""
def generate_haplotypes(snp_file,
                        haplotype_file,
                        vars,
                        inter_gap,
                        intra_gap,
                        num_genomes,
                        num_haplotypes):
    assert len(vars) > 0

    # daehwan - for debugging purposes
    # Sort variants and remove redundant variants
    vars = sorted(vars, cmp=compare_vars)
    tmp_vars = []
    v = 0
    while v < len(vars):
        var = vars[v]
        for v2 in range(v + 1, len(vars)):
            var2 = vars[v2]
            if compare_vars(var, var2) == 0:
                v += 1
                if "CLNSIG" not in var[4]:
                    if "CLNSIG" in var2[4]:
                        var[4]["CLNSIG"] = var2[4]["CLNSIG"]
                if "genotype" not in var[4]:
                    if "genotype" in var2[4]:
                        var[4]["genotype"] = var2[4]["genotype"]
            else:
                assert compare_vars(var, var2) < 0
                break
        tmp_vars.append(var)
        v += 1
    vars = tmp_vars

    # Assign genotypes for those missing genotypes
    max_genotype_num = 1
    genotypes_list = []
    for v in range(len(vars)):
        var = vars[v]
        var_dic = var[4]
        if "genotype" not in var_dic:
            used = [True, True] + [False for i in range(8)]
            v2 = v - 1
            while v2 >= 0:
                var2 = vars[v2]
                if compatible_vars(var2, var):
                    break
                var2_dic = var2[4]
                assert "genotype" in var2_dic
                genotype_num = int(var2_dic["genotype"][0])
                used[genotype_num] = True
                v2 -= 1

            assert False in used
            for i in range(len(used)):
                if not used[i]:                
                    var_dic["genotype"] = ("%d" % i) * (num_genomes * 2)
                    if i > max_genotype_num:
                        max_genotype_num = i
                    break
        genotypes_list.append(var_dic["genotype"])

    # daehwan - for debugging purposes
    """
    for v in range(len(vars)):
        var = vars[v]
        var_chr, var_pos, var_type, var_data, var_dic = var
        print v, var_chr, var_pos, var_type, var_data, var_dic["id"], var_dic["id2"],
        if "CLNSIG" in var_dic:
            print "CLNSIG:", var_dic["CLNSIG"],
        if "genotype" in var_dic:
            print var_dic["genotype"][:50],
        print
    """

    num_chromosomes = len(genotypes_list[0])

    # Write SNPs into a file (.snp)
    for var in vars:
        chr, pos, type, data, var_dic = var
        varID = var_dic["id2"]
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
    #    Var3: 222222222
    # Get haplotypes from genotypes_list
    haplotypes = set()
    cnv_genotypes = ["" for i in range(num_chromosomes)]
    for genotypes in genotypes_list:
        for i in range(len(genotypes)):
            genotype = genotypes[i]
            cnv_genotypes[i] += genotype

    cnv_genotypes = set(cnv_genotypes)
    for raw_haplotype in cnv_genotypes:
        for num in range(1, max_genotype_num + 1):
            num_str = str(num)
            if num_str not in raw_haplotype:
                continue
            haplotype = ""
            for i in range(len(raw_haplotype)):
                if raw_haplotype[i] == num_str:
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
                _, prev_locus, prev_type, prev_data, _ = vars[int(haplotype[s-1])]
                _, locus, type, data, _ = vars[int(haplotype[s])]
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
        _, a1_locus, _, _, _ = vars[int(a[0])]
        _, a2_locus, a2_type, a2_data, _ = vars[int(a[-1])]
        a_begin, a_end = int(a1_locus), int(a2_locus)
        if a2_type == 'D':
            a_end += (int(a2_data) - 1)
        b = b.split('#')
        _, b1_locus, _, _, _ = vars[int(b[0])]
        _, b2_locus, b2_type, b2_data, _ = vars[int(b[-1])]
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
        chr, h1_locus, _, _, _ = vars[int(h[0])]
        _, h2_locus, h2_type, h2_data, _ = vars[int(h[-1])]
        h_begin, h_end = int(h1_locus), int(h2_locus)
        if h2_type == 'D':
            h_end += (int(h2_data) - 1)
        assert h_begin <= h_end
        h_new_begin = h_begin
        for h_j in reversed(range(0, h_i)):
            hc = haplotypes[h_j].split('#')
            _, hc_begin, hc_type, hc_data, _ = vars[int(hc[-1])]
            hc_begin = int(hc_begin)
            hc_end = hc_begin
            if hc_type == 'D':
                hc_end += (int(hc_data) - 1)
            if hc_end + inter_gap < h_begin:
                break
            if h_new_begin > hc_end:
                h_new_begin = hc_end
        assert h_new_begin <= h_begin
        h_add = []
        for id in h:
            var_dic = vars[int(id)][4]
            h_add.append(var_dic["id2"])
        print >> haplotype_file, "ht%d\t%s\t%d\t%d\t%s" % \
            (num_haplotypes, chr, h_new_begin, h_end, ','.join(h_add))
        num_haplotypes += 1

    return num_haplotypes


"""
"""
def main(genome_file,
         base_fname,
         species,
         assembly_version,
         reference_type,
         genotype_vcf,
         genotype_gene_list,
         extra_files,
         inter_gap,
         intra_gap,
         verbose = False):
    # Load genomic sequences
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
            VCF_fnames["human"].append(["MT", "%s/ALL.chrMT.phase3_callmom.20130502.genotypes.vcf.gz" % grch37_url_base])
    
    if species not in VCF_fnames:
        print >> sys.stderr, "Error: %s is not supported." % species
        sys.exit(1)

    # List of variants (e.g. ClinVar database)
    genotype_var_list = {}
    # List of genomic regions to be processed
    genotype_ranges = {}
    if genotype_vcf != "":
        var_set = set()
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

            gene = None
            for g in genotype_gene_list:
                if info.find(g) != -1:
                    gene = g
                    break
            if not gene:
                continue

            for item in info.split(';'):
                if not item.startswith("CLNSIG"):
                    continue
                try:
                    key, value = item.split('=')
                    CLNSIG = int(value)
                except ValueError:
                    continue
            
            vars = extract_vars(chr_dic, chr, pos, ref_allele, alt_alleles, varID)
            if len(vars) == 0:
                continue

            if chr not in genotype_var_list:
                genotype_var_list[chr] = []
                genotype_ranges[chr] = {}
            if gene not in genotype_ranges[chr]:
                genotype_ranges[chr][gene] = [len(chr_dic[chr]), -1]

            for var in vars:
                var_chr, var_pos, var_ref_allele, var_alt_allele = var[:4]
                var_str = "%s-%d-%s-%s" % (var_chr, var_pos, var_ref_allele, var_alt_allele)
                if var_str in var_set:
                    continue
                var[4]["CLNSIG"] = CLNSIG
                genotype_var_list[chr].append(var)
                if var_pos < genotype_ranges[chr][gene][0]:
                    genotype_ranges[chr][gene][0] = var_pos
                if var_pos > genotype_ranges[chr][gene][1]:
                    genotype_ranges[chr][gene][1] = var_pos
                    
                var_set.add(var_str)

        print >> sys.stderr, "Number of variants in %s is:" % (genotype_vcf)
        for chr, vars in genotype_var_list.items():
            vars = sorted(vars, cmp=compare_vars)
            print >> sys.stderr, "\tChromosome %s: %d variants" % (chr, len(vars))

        for chr, gene_ranges in genotype_ranges.items():
            for gene, value in gene_ranges.items():
                gene_ranges[gene] = [value[0] - 100, value[1] + 100]
                value = genotype_ranges[chr][gene]
                print >> sys.stderr, "%s\t%s\t%d-%d" % (chr, gene, value[0], value[1])
            
        clnsig_file = open("%s.clnsig" % base_fname, 'w')
        for chr, vars in genotype_var_list.items():
            for var in vars:
                varID = var[4]["id2"]
                CLNSIG = var[4]["CLNSIG"]
                print >> clnsig_file, "%s\t%d" % (varID, CLNSIG)
        clnsig_file.close()

    SNP_file = open("%s.snp" % base_fname, 'w')
    haplotype_file = open("%s.haplotype" % base_fname, 'w')
    if extra_files:
        backbone_file = open("%s_backbone.fa" % base_fname, 'w')
        ref_file = open("%s.ref" % base_fname, 'w')

    num_haplotypes = 0
    num_unassigned = 0
    for chr, url in VCF_fnames[species]:
        if chr not in chr_dic:
            continue
        if reference_type != "genome" and \
                len(genotype_gene_list) > 0 and \
                chr not in genotype_var_list:
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

        chr_genotype_vars = []
        chr_genotype_ranges = {}
        if len(genotype_gene_list) > 0:
            assert chr in genotype_var_list
            chr_genotype_vars = genotype_var_list[chr]
            assert chr in genotype_ranges
            chr_genotype_ranges = genotype_ranges[chr]
        
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
                num_genomes = len(genomeIDs)
                continue

            assert len(genotypes) == len(genomeIDs)

            if not varID.startswith("rs") or \
                    ';' in varID:
                continue

            if varID == prev_varID:
                continue

            pos = int(pos) - 1

            # daehwan - for debugging purposes
            if pos >= 44000000:
                break

            if num_lines % 10000 == 1:
                print >> sys.stderr, "\t%s:%d\r" % (chr, pos),
            
            if pos == prev_pos:
                continue

            if chr_genotype_ranges:
                skip = True
                for gene, range_ in chr_genotype_ranges.items():
                    if pos > range_[0] and pos < range_[1]:
                        skip = False
                        break
                if skip:
                    continue
                if len(vars) == 0:
                    for var in chr_genotype_vars:
                        var_chr, var_pos, var_type, var_data, var_dic = var
                        if var_pos < range_[0]:
                            continue
                        if var_pos > range_[1]:
                            break
                        vars.append([var_chr, var_pos, var_type, var_data, var_dic])
                    curr_right = range_[1]

            def add_vars(pos,
                         varID,
                         ref_allele,
                         alt_alleles,
                         genotypes):
                tmp_vars = extract_vars(chr_dic, chr, pos, ref_allele, alt_alleles, varID)
                max_right = -1
                for v in range(len(tmp_vars)):
                    var = tmp_vars[v]
                    _, pos2, type, data = var[:4]
                    cnv_genotypes = []
                    for genotype in genotypes:
                        P1, P2 = genotype[0], genotype[2]
                        if P1 == digit2str[v + 1]:
                            cnv_genotypes.append('1')
                        else:
                            cnv_genotypes.append('0')
                        if P2 == digit2str[v + 1]:
                            cnv_genotypes.append('1')
                        else:
                            cnv_genotypes.append('0')

                    # Skip SNPs not present in a given population (e.g. 2,504 genomes in 1000 Genomes Project)
                    if '1' not in cnv_genotypes:
                        continue

                    tmp_varID = var[4]["id2"]
                    var_dic = {"id":varID, "id2":tmp_varID, "genotype":''.join(cnv_genotypes)}
                    vars.append([chr, pos2, type, data, var_dic])
                    right = pos2
                    if type == 'D':
                        right += (int(data) - 1)
                    if max_right < right:
                        max_right = right
                return max_right

            if curr_right + inter_gap < pos and len(vars) > 0:
                num_haplotypes = generate_haplotypes(SNP_file,
                                                     haplotype_file,
                                                     vars,
                                                     inter_gap,
                                                     intra_gap,
                                                     num_genomes,
                                                     num_haplotypes)
                vars = []
            right = add_vars(pos,
                             varID,
                             ref_allele,
                             alt_alleles,
                             genotypes)
            if curr_right < right:
                curr_right = right

            prev_varID = varID
            prev_pos = pos

        if len(vars) > 0:
            num_haplotypes = generate_haplotypes(SNP_file,
                                                 haplotype_file,
                                                 vars,
                                                 inter_gap,
                                                 intra_gap,
                                                 num_genomes,
                                                 num_haplotypes)
            vars = []

    SNP_file.close()
    haplotype_file.close()

    if extra_files:
        backbone_file.close()
        ref_file.close()

    if genotype_vcf != "":
        clnsig_file.close()
        


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
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="genome",
                        help="Reference type: gene, chromosome, and genome (default: genome)")
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
    parser.add_argument('--extra-files',
                        dest='extra_files',
                        action='store_true',
                        help='Output extra files such as _backbone.fa and .ref')
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
         args.reference_type,
         args.genotype_vcf,
         args.genotype_gene_list,
         args.extra_files,
         args.inter_gap,
         args.intra_gap,
         args.verbose)
