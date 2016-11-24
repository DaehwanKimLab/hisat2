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


import os, sys, subprocess, re
import inspect
from argparse import ArgumentParser, FileType


"""
"""
def read_genome(genome_file):
    chr_dic, chr_names, chr_full_names = {}, [], []
    chr_name, chr_full_name, sequence = "", "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
                chr_names.append(chr_name)
            chr_full_name = line.strip()[1:]
            chr_name = line.strip().split()[0][1:]
            chr_full_names.append(chr_full_name)
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
        chr_names.append(chr_name)
        chr_full_names.append(chr_full_name)
    return chr_dic, chr_names, chr_full_names


"""
"""
def read_sequences(fname):
    allele_seqs = {}
    allele_name, sequence = "", ""
    for line in open(fname):
        if line.startswith(">"):
            if allele_name != "" and allele_name not in allele_seqs:
                allele_seqs[allele_name] = sequence
            allele_name = line.strip()[1:]
            sequence = ""
        else:
            sequence += line.strip()
    if allele_name != "" and allele_name not in allele_seqs:
        allele_seqs[allele_name] = sequence
    return allele_seqs


"""
"""
def read_variants(fname):
    allele_vars = {}
    for line in open(fname):
        var_id, type, allele_name, left, data = line.strip().split()
        left = int(left)
        if type == "deletion":
            data = int(data)
        if allele_name not in allele_vars:
            allele_vars[allele_name] = []
        allele_vars[allele_name].append([left, type, data, var_id])
    return allele_vars


"""
"""
def read_haplotypes(fname):
    allele_haplotypes = {}
    for line in open(fname):
        haplotype_id, allele_name, left, right, vars = line.strip().split()
        vars = vars.split(',')
        left, right = int(left), int(right)
        if allele_name not in allele_haplotypes:
            allele_haplotypes[allele_name] = []
        allele_haplotypes[allele_name].append([left, right, vars])
    return allele_haplotypes


"""
"""
def read_links(fname):
    links = []
    for line in open(fname):
        var_id, allele_names = line.strip().split('\t')
        links.append([var_id, allele_names])
    return links


"""
"""
def read_clnsig(fname):
    clnsig_dic = {}
    for line in open(fname):
        var_id, gene, clnsig = line.strip().split('\t')
        clnsig_dic[var_id] = [gene, clnsig]
    return clnsig_dic


"""
Compare two variants
"""
def compare_vars(a, b):
    a_pos, a_type, a_data = a[:3]
    b_pos, b_type, b_data = b[:3]

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
def build_genotype_genome(reference,
                          base_fname,                          
                          partial,
                          inter_gap,
                          intra_gap,
                          threads,
                          use_clinvar,
                          verbose):    
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(build_genotype_genome))
    ex_path = os.path.dirname(curr_script)

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")

    # Load genomic sequences
    chr_dic, chr_names, chr_full_names = read_genome(open(reference))

    if use_clinvar:
        # Extract variants from the ClinVar database
        CLINVAR_fnames = ["clinvar.vcf.gz",
                          "clinvar.snp",
                          "clinvar.haplotype",
                          "clinvar.clnsig"]

        if not check_files(CLINVAR_fnames):
            if not os.path.exists("clinvar.vcf.gz"):
                os.system("wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz")
            assert os.path.exists("clinvar.vcf.gz")

            extract_clinvar_script = os.path.join(ex_path, "hisat2_extract_snps_haplotypes_VCF.py")
            extract_cmd = [extract_clinvar_script]
            extract_cmd += ["--inter-gap", str(inter_gap),
                            "--intra-gap", str(intra_gap),
                            "--genotype-vcf", "clinvar.vcf.gz",
                            reference, "/dev/null", "clinvar"]
            if verbose:
                print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
            proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()
            if not check_files(CLINVAR_fnames):
                print >> sys.stderr, "Error: extract variants from clinvar failed!"
                sys.exit(1)

        # Read variants to be genotyped
        genotype_vars = read_variants("clinvar.snp")

        # Read haplotypes
        genotype_haplotypes = read_haplotypes("clinvar.haplotype")

        # Read information about clinical significance
        genotype_clnsig = read_clnsig("clinvar.clnsig")
    else:
        genotype_vars, genotype_haplotypes, genotype_clnsig = {}, {}, {}

    # Genes to be genotyped
    genotype_genes = {}

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        os.system("git clone https://github.com/jrob119/IMGTHLA.git")

    # Extract HLA variants, backbone sequence, and other sequeces
    HLA_fnames = ["hla_backbone.fa",
                  "hla.ref",
                  "hla.snp",
                  "hla.haplotype",
                  "hla.link"]

    if not check_files(HLA_fnames):
        extract_hla_script = os.path.join(ex_path, "hisat2_extract_HLA_vars.py")
        extract_cmd = [extract_hla_script]
        if partial:
            extract_cmd += ["--partial"]
        extract_cmd += ["--inter-gap", str(inter_gap),
                        "--intra-gap", str(intra_gap)]
        if verbose:
            print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: extract_HLA_vars failed!"
            sys.exit(1)

    # Read HLA genes
    if os.path.exists("hla.ref"):
        for line in open("hla.ref"):
            HLA_name, chr, left, right, length, exon_str = line.strip().split()
            left, right = int(left), int(right)
            length = int(length)
            if chr not in chr_names:
                continue
            if chr not in genotype_genes:
                genotype_genes[chr] = []
            genotype_genes[chr].append([left, right, length, HLA_name, "hla"])

    # Write genotype genome
    var_num, haplotype_num = 0, 0
    genome_out_file = open("%s.fa" % base_fname, 'w')
    gene_out_file = open("%s.gene" % base_fname, 'w')
    var_out_file = open("%s.snp" % base_fname, 'w')
    haplotype_out_file = open("%s.haplotype" % base_fname, 'w')
    link_out_file = open("%s.link" % base_fname, 'w')
    coord_out_file = open("%s.coord" % base_fname, 'w')
    clnsig_out_file = open("%s.clnsig" % base_fname, 'w')    
    for c in range(len(chr_names)):
        chr = chr_names[c]
        chr_full_name = chr_full_names[c]
        assert chr in chr_dic
        chr_seq = chr_dic[chr]
        chr_len = len(chr_seq)
        if chr in genotype_genes:
            chr_genes = genotype_genes[chr]
            def gene_cmp(a, b):
                a_left, a_right, a_length = a[:3]
                b_left, b_right, b_length = b[:3]
                if a_left != b_left:
                    return a_left - b_left
                if a_right != b_right:
                    return a_right - b_right
                return a_lenght - b_length
            chr_genes = sorted(chr_genes, cmp=gene_cmp)
        else:
            chr_genes = []

        chr_genotype_vars, chr_genotype_vari = [], 0
        if chr in genotype_vars:
            chr_genotype_vars = genotype_vars[chr]
        chr_genotype_haplotypes, chr_genotype_hti = [], 0
        if chr in genotype_haplotypes:
            chr_genotype_haplotypes = genotype_haplotypes[chr]

        def add_vars(left, right, chr_genotype_vari, chr_genotype_hti, haplotype_num):
            # Output variants with clinical significance
            while chr_genotype_vari < len(chr_genotype_vars):
                var_left, var_type, var_data, var_id =  chr_genotype_vars[chr_genotype_vari]
                var_right = var_left
                if var_type == "deletion":
                    var_right += var_data
                if var_right > right:
                    break
                if var_right >= left:
                    continue

                print >> var_out_file, "%s\t%s\t%s\t%d\t%s" % \
                    (var_id, var_type, chr, var_left + off, var_data)

                assert var_id in genotype_clnsig
                var_gene, clnsig = genotype_clnsig[var_id]
                print >> clnsig_out_file, "%s\t%s\t%s" % \
                    (var_id, var_gene, clnsig)
                
                chr_genotype_vari += 1

            # Output haplotypes
            while chr_genotype_hti < len(chr_genotype_haplotypes):
                ht_left, ht_right, ht_vars =  chr_genotype_haplotypes[chr_genotype_hti]
                if ht_right > right:
                    break
                if ht_right >= left:
                    continue

                print >> haplotype_out_file, "ht%d\t%s\t%d\t%d\t%s" % \
                    (haplotype_num, chr, ht_left + off, ht_right + off, ','.join(ht_vars))
                chr_genotype_hti += 1
                haplotype_num += 1

            return chr_genotype_vari, chr_genotype_hti, haplotype_num

        out_chr_seq = ""
        
        off = 0
        prev_right = 0
        for gene in chr_genes:
            left, right, length, name, family = gene

            chr_genotype_vari, chr_genotype_hti, haplotype_num = add_vars(left, right, chr_genotype_vari, chr_genotype_hti, haplotype_num)

            # Read HLA backbone sequences
            allele_seqs = read_sequences("%s_backbone.fa" % family)

            # Read HLA variants
            allele_vars = read_variants("%s.snp" % family)

            # Read HLA haplotypes
            allele_haplotypes = read_haplotypes("%s.haplotype" % family)

            # Read HLA link information between haplotypes and variants
            links = read_links("%s.link" % family)

            if name not in allele_seqs or \
                    name not in allele_vars or \
                    name not in allele_haplotypes:
                continue
            allele_seq = allele_seqs[name]
            vars = allele_vars[name]
            haplotypes = allele_haplotypes[name]
            assert length == len(allele_seq)
            assert left < chr_len and right < chr_len
            # Skipping overlapping genes
            if left < prev_right:
                print >> sys.stderr, "Warning: skipping %s ..." % (name)
                continue

            varID2htID = {}

            assert left < right
            prev_length = right - left + 1
            assert prev_length <= length

            if prev_right < left:
                # print >> coord_out_file, "%d\t%d\t%d" % \
                #    (len(out_chr_seq), prev_right, left - prev_right)
                out_chr_seq += chr_seq[prev_right:left]

            # Output gene (genotype_genome.gene)
            print >> gene_out_file, "%s\t%s\t%s\t%d\t%d" % \
                (family.upper(), name, chr, len(out_chr_seq), len(out_chr_seq) + length - 1)

            # Output coord (genotype_genome.coord)
            print >> coord_out_file, "%s\t%d\t%d\t%d" % \
                (chr, len(out_chr_seq), left, right - left + 1)
            out_chr_seq += allele_seq

            # Output variants (genotype_genome.snp)
            for var in vars:
                var_left, var_type, var_data, var_id = var
                new_var_id = "hv%d" % var_num
                varID2htID[var_id] = new_var_id
                new_var_left = var_left + left + off
                assert var_type in ["single", "deletion"]
                assert new_var_left < len(out_chr_seq)
                if var_type == "single":                    
                    assert out_chr_seq[new_var_left] != var_data
                else:
                    assert var_type == "deletion"
                    assert new_var_left + var_data <= len(out_chr_seq)
                    
                print >> var_out_file, "%s\t%s\t%s\t%d\t%s" % \
                    (new_var_id, var_type, chr, new_var_left, var_data)
                var_num += 1
                
            # Output haplotypes (genotype_genome.haplotype)
            for haplotype in haplotypes:
                ht_left, ht_right, ht_vars = haplotype
                new_ht_left = ht_left + left + off
                assert new_ht_left < len(out_chr_seq)
                new_ht_right = ht_right + left + off
                assert new_ht_left <= new_ht_right
                assert new_ht_right <= len(out_chr_seq)
                new_ht_vars = []
                for var_id in ht_vars:
                    assert var_id in varID2htID
                    new_ht_vars.append(varID2htID[var_id])
                print >> haplotype_out_file, "ht%d\t%s\t%d\t%d\t%s" % \
                    (haplotype_num, chr, new_ht_left, new_ht_right, ','.join(new_ht_vars))
                haplotype_num += 1

            # Output link information between alleles and variants (genotype_genome.link)
            for link in links:
                var_id, allele_names = link
                if var_id not in varID2htID:
                    continue
                new_var_id = varID2htID[var_id]
                print >> link_out_file, "%s\t%s" % (new_var_id, allele_names)
                
            off += (length - prev_length)

            prev_right = right + 1

        # Write the rest of the Vars
        chr_genotype_vari, chr_genotype_hti, haplotype_num = add_vars(10000000000, 10000000000, chr_genotype_vari, chr_genotype_hti, haplotype_num)            
            
        print >> coord_out_file, "%s\t%d\t%d\t%d" % \
            (chr, len(out_chr_seq), prev_right, len(chr_seq) - prev_right)
        out_chr_seq += chr_seq[prev_right:]

        assert len(out_chr_seq) == len(chr_seq) + off

        # Output chromosome sequence
        print >> genome_out_file, ">%s" % (chr_full_name)
        line_width = 60
        for s in range(0, len(out_chr_seq), line_width):
            print >> genome_out_file, out_chr_seq[s:s+line_width]

    genome_out_file.close()
    gene_out_file.close()
    var_out_file.close()
    haplotype_out_file.close()
    link_out_file.close()
    coord_out_file.close()
    clnsig_out_file.close()

    # Build HISAT2 graph indexes based on the above information
    hisat2_index_fnames = ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]
    hisat2_build = os.path.join(ex_path, "hisat2-build")
    build_cmd = [hisat2_build,
                 "-p", str(threads),
                 "--snp", "%s.snp" % base_fname,
                 "--haplotype", "%s.haplotype" % base_fname,
                 "%s.fa" % base_fname,
                 "%s" % base_fname]
    if verbose:
        print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
    proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
    proc.communicate()        
    if not check_files(hisat2_index_fnames):
        print >> sys.stderr, "Error: indexing failed!  Perhaps, you may have forgotten to build hisat2 executables?"
        sys.exit(1)


        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Extract HLA variants from HLA multiple sequence alignments")
    parser.add_argument("reference",
                        nargs='?',
                        type=str,
                        help="Reference genome")
    parser.add_argument("base_fname",
                        nargs='?',
                        type=str,
                        help="base filename for genotype genome")
    parser.add_argument('--partial',
                        dest='partial',
                        action='store_true',
                        help='Include partial alleles (e.g. A_nuc.fasta)')
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
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads") 
    parser.add_argument("--no-clinvar",
                        dest="use_clinvar",
                        action="store_false",
                        help="")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if not args.reference or not args.base_fname:
        parser.print_help()
        sys.exit(1)
    if args.inter_gap > args.intra_gap:
        print >> sys.stderr, "Error: --inter-gap (%d) must be smaller than --intra-gap (%d)" % (args.inter_gap, args.intra_gap)
        sys.exit(1)
    build_genotype_genome(args.reference,
                          args.base_fname,
                          args.partial,
                          args.inter_gap,
                          args.intra_gap,
                          args.threads,
                          args.use_clinvar,
                          args.verbose)
