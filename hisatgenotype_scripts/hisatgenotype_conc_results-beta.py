#!/usr/bin/env python
#
# Copyright 2019, Christopher Bennett <christopher@bennett-tech.dev>
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

import os, sys, glob
from argparse import ArgumentParser


def build_tree(vlist, tree, leaf):
    if len(vlist) == 0:
        return {'score' : leaf, 'children' : None}

    field = vlist[0]
    if field not in tree['children']:
        tree['children'][field] = build_tree(vlist[1:], {'score' : 0, 'children' : {}}, leaf)
    else:
        tree['children'][field] = build_tree(vlist[1:], tree['children'][field], leaf)

    tree['score'] += leaf
    return tree

def call_nuance_results(nfile):
    datatree = { 'EM' : {} , 'Allele splitting' : {}, 'Assembly' : {} }

    viterbi = False
    with open(nfile, "r") as ifi:
        for line in ifi:
            line = line.strip()
            if line.startswith("Assembly"):
                viterbi = True
                continue

            if 'abundance' not in line and not viterbi:
                continue

            if viterbi:
                ix = line.find(':')
                datatree['Assembly'][line[:ix]] = line[ix+2:]
                continue

            if '***' in line:
                split_by = 3
            else:
                split_by = 2

            gene = line.split()[split_by].split('*')[0]
            ix = line.find(gene)
            line = line[ix:]
            if gene not in datatree['EM']:
                datatree['EM'][gene] = []
                datatree['Allele splitting'][gene] = {'score' : 0, 'children' : {}}

            datatree['EM'][gene].append(line)

            replacement = ['(', ')']
            for sym in replacement:
                line = line.replace(sym, '')

            allele, _, percent = line.split()

            allele = allele.split('*')[-1].split(':')
            datatree['Allele splitting'][gene] = build_tree(allele, datatree['Allele splitting'][gene], round(float(percent[:-1])/100,4))

    return datatree


def flatten(tree, prev_key = '', sep = '*'):
    items = []
    for key, value in tree.items():
        new_key = prev_key + sep + key if prev_key else key
        try:
            items.extend(flatten(value['children'], new_key, ':').items())
            items.append((new_key + ' - Partial', value['score']))
        except:
            items.append((new_key, value['score']))

    if sep == ":":
        return dict(items)
    else:
        return sorted(items, key=lambda tup : (tup[1], len(tup[0].split()[0])), reverse = True)

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Script for simplifying HISAT-genotype results')

    parser.add_argument('--in-dir',
                        dest="read_dir",
                        type=str,
                        default=".",
                        help='Input directory (e.g. read_input) (default: (empty))')
    parser.add_argument("--csv",
                        dest="csv",
                        action = "store_true",
                        help='Save Results as CSV dataframe')

    args = parser.parse_args()

    if args.read_dir:
        indir = args.read_dir
    else:
        indir = '.'

    reports = glob.glob('%s/*.report' % indir)

    report_results = {}
    for report in reports:
        report_results[report] = call_nuance_results(report)

    scores = []
    genecol = ["File"]
    for file_ in report_results:
        print('File: %s' % file_)
        scores.append([file_])

        for type_ in report_results[file_]:
            print("\tAnalysis - %s" % type_)
            if type_ == 'Allele splitting':
                tree = report_results[file_]['Allele splitting']
                for gene in tree:
                    collab = "%s: %s" % (type_, gene)
                    if collab not in genecol:
                        genecol.append(collab)                    

                    print('\t\tGene: %s (score: %.2f)' % (gene, tree[gene]['score']))
                    flattened_tree = flatten(tree[gene]['children'], gene)

                    pastv, pastn = 0, ''
                    report_line = ''
                    for tup in flattened_tree:
                        # This is to filter out any missilaneous partial alleles with similar score to parent (ex remove A*01 - Partial 50% if A*01:01:01:01 50%)
                        if tup[1] < 0.2 or (pastv == tup[1] and "Partial" in tup[0]): 
                            continue

                        result_str = '%s (score: %.4f)' % (tup[0], tup[1])
                        report_line += result_str + ","

                        print("\t\t\t" + result_str)
                        pastn, pastv = tup
                    
                    scores[-1].append(report_line[:-1])
            
            else:
                for gene, data in report_results[file_][type_].items():
                    collab = "%s: %s" % (type_, gene)
                    if collab not in genecol:
                        genecol.append(collab)

                    print('\t\tGene: %s' % gene)

                    if isinstance(data, list):
                        scores[-1].append(",".join(data))
                        for line in data:
                            print('\t\t\t%s' % line)
                    else:
                        scores[-1].append(data)
                        print('\t\t\t%s' % data)

    if args.csv:
        scores.insert(0, genecol)
        with open("HG_report_results.csv", "w") as ofo:
            for line in scores:
                line = "\t".join(line)
                ofo.write(line + "\n")
