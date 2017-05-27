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
import hisatgenotype_typing_common as typing_common, hisatgenotype_gene_typing as gene_typing


"""
"""
def read_clnsig(fname):
    clnsig_dic = {}
    for line in open(fname):
        var_id, gene, clnsig = line.strip().split('\t')
        clnsig_dic[var_id] = [gene, clnsig]
    return clnsig_dic


"""
"""
def build_genotype_genome(base_fname,                          
                          inter_gap,
                          intra_gap,
                          threads,
                          database_list,
                          use_clinvar,
                          use_commonvar,
                          verbose):    
    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not typing_common.check_files(HISAT2_fnames):
        typing_common.download_genome_and_index()

    # Load genomic sequences
    chr_dic, chr_names, chr_full_names = typing_common.read_genome(open("genome.fa"))

    genotype_vars, genotype_haplotypes, genotype_clnsig = {}, {}, {}
    if use_clinvar:
        # Extract variants from the ClinVar database
        CLINVAR_fnames = ["clinvar.vcf.gz",
                          "clinvar.snp",
                          "clinvar.haplotype",
                          "clinvar.clnsig"]

        if not typing_common.check_files(CLINVAR_fnames):
            if not os.path.exists("clinvar.vcf.gz"):
                os.system("wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive/2017/clinvar_20170404.vcf.gz")
            assert os.path.exists("clinvar.vcf.gz")

            extract_cmd = ["hisat2_extract_snps_haplotypes_VCF.py"]
            extract_cmd += ["--inter-gap", str(inter_gap),
                            "--intra-gap", str(intra_gap),
                            "--genotype-vcf", "clinvar.vcf.gz",
                            "genome.fa", "/dev/null", "clinvar"]
            if verbose:
                print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
            proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()
            if not typing_common.check_files(CLINVAR_fnames):
                print >> sys.stderr, "Error: extract variants from clinvar failed!"
                sys.exit(1)

        # Read variants to be genotyped
        genotype_vars = typing_common.read_variants("clinvar.snp")

        # Read haplotypes
        genotype_haplotypes = typing_common.read_haplotypes("clinvar.haplotype")

        # Read information about clinical significance
        genotype_clnsig = typing_common.read_clnsig("clinvar.clnsig")

    if use_commonvar:
        # Extract variants from dbSNP database
        commonvar_fbase = "snp144Common"
        commonvar_fnames = ["%s.snp" % commonvar_fbase,
                            "%s.haplotype" % commonvar_fbase]
        if not typing_common.check_files(commonvar_fnames):
            if not os.path.exists("%s.txt.gz" % commonvar_fbase):
                os.system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/%s.txt.gz" % commonvar_fbase)
            assert os.path.exists("%s.txt.gz" % commonvar_fbase)
            os.system("gzip -cd %s.txt.gz | awk 'BEGIN{OFS=\"\t\"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == \"M\") {$2 = \"MT\"} print}' > %s.txt" % (commonvar_fbase, commonvar_fbase))
            extract_cmd = ["hisat2_extract_snps_haplotypes_UCSC.py",
                           "--inter-gap", str(inter_gap),
                           "--intra-gap", str(intra_gap),
                           "genome.fa", "%s.txt" % commonvar_fbase, commonvar_fbase]
            if verbose:
                print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
            proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()
            if not typing_common.check_files(commonvar_fnames):
                print >> sys.stderr, "Error: extract variants from clinvar failed!"
                sys.exit(1)

        # Read variants to be genotyped
        genotype_vars = typing_common.read_variants("%s.snp" % commonvar_fbase)

        # Read haplotypes
        genotype_haplotypes = typing_common.read_haplotypes("%s.haplotype" % commonvar_fbase)

    # Genes to be genotyped
    genotype_genes = {}

    # Read genes or genomics regions
    for database_name in database_list:
        # Extract HLA variants, backbone sequence, and other sequeces
        typing_common.extract_database_if_not_exists(database_name,
                                                     [],            # locus_list
                                                     inter_gap,
                                                     intra_gap,
                                                     True,          # partial?
                                                     verbose)
        locus_fname = "%s.locus" % database_name
        assert os.path.exists(locus_fname)
        for line in open(locus_fname):
            HLA_name, chr, left, right, length, exon_str, strand = line.strip().split()
            left, right = int(left), int(right)
            length = int(length)
            if chr not in chr_names:
                continue
            if chr not in genotype_genes:
                genotype_genes[chr] = []
            genotype_genes[chr].append([left, right, length, HLA_name, database_name, exon_str, strand])

    # Write genotype genome
    var_num, haplotype_num = 0, 0
    genome_out_file = open("%s.fa" % base_fname, 'w')
    locus_out_file = open("%s.locus" % base_fname, 'w')
    var_out_file = open("%s.snp" % base_fname, 'w')
    index_var_out_file = open("%s.index.snp" % base_fname, 'w')
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
                    chr_genotype_vari += 1
                    continue

                out_str = "%s\t%s\t%s\t%d\t%s" % (var_id, var_type, chr, var_left + off, var_data)
                print >> var_out_file, out_str
                print >> index_var_out_file, out_str

                if var_id in genotype_clnsig:
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
                    chr_genotype_hti += 1
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
            left, right, length, name, family, exon_str, strand = gene

            chr_genotype_vari, chr_genotype_hti, haplotype_num = add_vars(left, right, chr_genotype_vari, chr_genotype_hti, haplotype_num)

            # Read HLA backbone sequences
            allele_seqs = typing_common.read_allele_sequences("%s_backbone.fa" % family)

            # Read HLA variants
            allele_vars = typing_common.read_variants("%s.snp" % family)
            allele_index_vars = typing_common.read_variants("%s.index.snp" % family)
                
            # Read HLA haplotypes
            allele_haplotypes = typing_common.read_haplotypes("%s.haplotype" % family)

            # Read HLA link information between haplotypes and variants
            links = typing_common.read_links("%s.link" % family)

            if name not in allele_seqs or \
                    name not in allele_vars or \
                    name not in allele_haplotypes:
                continue
            allele_seq = allele_seqs[name]
            vars, index_vars = allele_vars[name], allele_index_vars[name]
            index_var_ids = set()
            for _, _, _, var_id in index_vars:
                index_var_ids.add(var_id)

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
                out_chr_seq += chr_seq[prev_right:left]

            # Output gene (genotype_genome.gene)
            print >> locus_out_file, "%s\t%s\t%s\t%d\t%d\t%s\t%s" % \
                (family.upper(), name, chr, len(out_chr_seq), len(out_chr_seq) + length - 1, exon_str, strand)

            # Output coord (genotype_genome.coord)
            print >> coord_out_file, "%s\t%d\t%d\t%d" % \
                (chr, len(out_chr_seq), left, right - left + 1)
            out_chr_seq += allele_seq

            # Output variants (genotype_genome.snp and genotype_genome.index.snp)
            for var in vars:
                var_left, var_type, var_data, var_id = var
                new_var_id = "hv%d" % var_num
                varID2htID[var_id] = new_var_id
                new_var_left = var_left + left + off
                assert var_type in ["single", "deletion", "insertion"]
                assert new_var_left < len(out_chr_seq)
                if var_type == "single":                    
                    assert out_chr_seq[new_var_left] != var_data
                elif var_type == "deletion":
                    assert new_var_left + var_data <= len(out_chr_seq)
                else:
                    assert var_type == "insertion"

                out_str = "%s\t%s\t%s\t%d\t%s" % (new_var_id, var_type, chr, new_var_left, var_data)
                print >> var_out_file, out_str
                if var_id in index_var_ids:
                    print >> index_var_out_file, out_str
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
        chr_genotype_vari, chr_genotype_hti, haplotype_num = add_vars(sys.maxint, sys.maxint, chr_genotype_vari, chr_genotype_hti, haplotype_num)            
            
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
    locus_out_file.close()
    var_out_file.close()
    index_var_out_file.close()
    haplotype_out_file.close()
    link_out_file.close()
    coord_out_file.close()
    clnsig_out_file.close()

    partial_out_file = open("%s.partial" % base_fname, 'w')
    for database in database_list:
        for line in open("%s.partial" % database):
            allele_name = line.strip()
            print >> partial_out_file, "%s\t%s" % (database.upper(), allele_name)
    partial_out_file.close()

    # Index genotype_genome.fa
    index_cmd = ["samtools", "faidx", "%s.fa" % base_fname]
    subprocess.call(index_cmd)

    # Build HISAT-genotype graph indexes based on the above information
    hisat2_index_fnames = ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]
    build_cmd = ["hisat2-build",
                 "-p", str(threads),
                 "--snp", "%s.index.snp" % base_fname,
                 "--haplotype", "%s.haplotype" % base_fname,
                 "%s.fa" % base_fname,
                 "%s" % base_fname]
    if verbose:
        print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
        
    subprocess.call(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
    if not typing_common.check_files(hisat2_index_fnames):
        print >> sys.stderr, "Error: indexing failed!  Perhaps, you may have forgotten to build hisat2 executables?"
        sys.exit(1)

        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Build genotype genome")
    parser.add_argument("--base", "--base-fname",
                        dest="base_fname",
                        type=str,
                        default="genotype_genome",
                        help="base filename for genotype genome (default: genotype_genome)")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--database-list",
                        dest="database_list",
                        type=str,
                        default="",
                        help="A comma-separated list of databases (default: hla,codis,cyp)")
    parser.add_argument("--commonvar",
                        dest="use_commonvar",
                        action="store_true",
                        help="Include common variants from dbSNP")
    parser.add_argument("--clinvar",
                        dest="use_clinvar",
                        action="store_true",
                        help="Include variants from ClinVar database")
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
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.inter_gap > args.intra_gap:
        print >> sys.stderr, "Error: --inter-gap (%d) must be smaller than --intra-gap (%d)" % (args.inter_gap, args.intra_gap)
        sys.exit(1)
        
    if args.database_list == "":
        database_list = []
    else:
        database_list = args.database_list.split(',')

    if args.use_clinvar and args.use_commonvar:
        print >> sys.stderr, "Error: both --clinvar and --commonvar cannot be used together."
        sys.exit(1)
        
        
    build_genotype_genome(args.base_fname,
                          args.inter_gap,
                          args.intra_gap,
                          args.threads,
                          database_list,
                          args.use_clinvar,
                          args.use_commonvar,
                          args.verbose)
    
