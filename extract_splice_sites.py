#!/usr/bin/env python

#
# Copyright 2014, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT.
#
# HISAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
#

import os, sys

"""
"""
def extract_transcripts(gtf_filename, verbose = True):
    genes = {}
    trans = {}
    tran_ids = []

    gtf_file = open(gtf_filename)
    for line in gtf_file:
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        chr, protein, type, left, right, comma1, strand, comma2, values = line[:-1].split('\t')
            
        values = values.strip().split(';')
        left, right = int(left), int(right)

        if type != "exon":
            continue

        if left >= right:
            continue

        value_dic = {}
        for value in values[:-1]:
            name, value = value.strip().split()
            value = value[1:-1]
            value_dic[name] = value

        if 'gene_id' not in value_dic or \
                'transcript_id' not in value_dic:
            continue

        gene_id = value_dic['gene_id']
        if gene_id not in genes:
            genes[gene_id] = []

        transcript_id = value_dic['transcript_id']
        if transcript_id not in trans:
            trans[transcript_id] = [chr, strand, [[left, right]]]
            tran_ids.append(transcript_id)
            genes[gene_id].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    gtf_file.close()

    for tran_id in tran_ids:
        chr, strand, exons = trans[tran_id]
        assert len(exons) >= 1
        def exon_cmp(a, b):
            if a[0] != b[0]:
                if a[0] < b[0]:
                    return -1
                else:
                    return 1

            if a[1] != b[1]:
                if a[1] < b[1]:
                    return -1
                else:
                    return 1

            assert False
            return 0

        exons = sorted(exons, cmp=exon_cmp)
        temp_exons = [exons[0]]
        for i in range(1, len(exons)):
            prev_left, prev_right = temp_exons[-1]
            left, right = exons[i]

            assert prev_left < left
            if left - prev_right <= 5:
                temp_exons[-1][1] = right
            else:
                temp_exons.append([left, right])
        trans[tran_id] = [chr, strand, temp_exons]

    if verbose:
        transcript_lengths, transcript_counts, transcript_length_sum = [0 for i in range(100000)], 0, 0
        exon_lengths, exon_counts, exon_length_sum = [0 for i in range(100000)], 0, 0
        intron_lengths, intron_counts, intron_length_sum = [0 for i in range(10000000)], 0, 0

        for tran_id in tran_ids:
            chr, strand, exons = trans[tran_id]
            transcript_length = 0
            for i in range(len(exons)):
                exon_start, exon_end = exons[i]
                exon_length = exon_end - exon_start + 1
                if exon_length >= len(exon_lengths):
                    continue

                exon_lengths[exon_length] += 1
                exon_counts += 1
                exon_length_sum += exon_length

                transcript_length += exon_length

                if i == 0:
                    continue

                intron_length = exon_start - exons[i-1][1]
                if intron_length >= len(intron_lengths):
                    continue

                intron_lengths[intron_length] += 1
                intron_counts += 1
                intron_length_sum += intron_length

            transcript_counts += 1
            transcript_length_sum += transcript_length

        gene_counts, gene_multi_isoforms_counts = 0, 0
        for g_id, t_ids in genes.items():
            gene_counts += 1
            if len(t_ids) > 1:
                gene_multi_isoforms_counts += 1

        transcript_avg_length = float(transcript_length_sum) / max(transcript_counts, 1)
        exon_avg_length = float(exon_length_sum) / max(exon_counts, 1)
        intron_avg_length = float(intron_length_sum) / max(intron_counts, 1)


        print >> sys.stderr, "gene counts: %d, gene counts (multiple isoforms): %d (%.2f)" % (gene_counts, gene_multi_isoforms_counts, float(gene_multi_isoforms_counts) / gene_counts)
        print >> sys.stderr, "transcript counts: %d, transcript avg. length: %.2f" % (transcript_counts, transcript_avg_length)
        print >> sys.stderr, "exon counts: %d, exon avg. length: %.2f" % (exon_counts, exon_avg_length)
        print >> sys.stderr, "intron counts: %d, intron avg. length: %.2f" % (intron_counts, intron_avg_length)
        print >> sys.stderr, "average number of exons per transcript: %.2f" % (float(exon_counts) / transcript_counts)
        
        for i in range(20):
            sum = 0
            for j in range(i * 10 + 1, (i+1) * 10 + 1):
                sum += exon_lengths[j]
                
            print >> sys.stderr, "[%d, %d] : %.2f" % (i * 10 + 1, (i+1) * 10, float(sum) * 100 / exon_counts)
            
    return trans, tran_ids


"""
"""
def extract_junctions(trans_dic):
    junctions_dic = {}
    for tran_id, values in trans_dic.items():
        chr, strand, exons = values
        for i in range(1, len(exons)):
            left, right = exons[i-1][1], exons[i][0]
            junction = "%s:%d:%d:%s" % (chr, left, right, strand)

            if junction not in junctions_dic:
                junctions_dic[junction] = []

            junctions_dic[junction].append(tran_id)
            

    return junctions_dic


"""
"""
def extract_splice_sites(gtf_fname):
    if not os.path.exists(gtf_fname):
        print >> sys.stderr, ""
        sys.exit(1)

    trans_dic, trans_ids = extract_transcripts(gtf_fname, verbose = False)
    junctions_dic = extract_junctions(trans_dic)
    junctions = []

    for junction, value in junctions_dic.items():
        chr, left, right, strand = junction.split(":")
        left, right = int(left) - 1, int(right) - 1
        assert left >= 0 and right >= 0 and right > left
        junctions.append([chr, left, right, strand])

    def junction_cmp(a, b):
        if a[0] != b[0]:
            if a[0] < b[0]:
                return -1
            else:
                return 1

        if a[1] != b[1]:
            if a[1] < b[1]:
                return -1
            else:
                return 1

        if a[2] != b[2]:
            if a[2] < b[2]:
                return -1
            else:
                return 1

        if a[3] != b[3]:
            if a[3] == "+":
                return -1
            else:
                return 1

        assert False
        return 0


    junctions = sorted(junctions, cmp=junction_cmp)
    for junction in junctions:
        chr, left, right, strand = junction
        print "%s\t%d\t%d\t%s" % (chr, left, right, strand)

    
"""
"""
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print >> sys.stderr, "Usage:"
        usage = "\tpython extract_splice_sites.py in_gtf_filename > out_splice_site_filename"
        print >> sys.stderr, usage
        sys.exit(1)
        
    extract_splice_sites(sys.argv[1])
