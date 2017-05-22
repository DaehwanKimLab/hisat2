#!/usr/bin/env python

#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>
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
import inspect, operator
from argparse import ArgumentParser, FileType

# sequences for DNA fingerprinting loci are available at http://www.cstl.nist.gov/biotech/strbase/seq_ref.htm

CODIS_loci = ["CSF1PO", "FGA", "TH01", "TPOX", "VWA", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11"]


"""
## Download variant information from website
"""
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


"""
Download CODIS data
"""
def download_codis(base_dname,
                   base_fname,
                   locus_list,
                   verbose):    
    # CODIS database base URL
    base_url = "http://www.cstl.nist.gov/biotech/strbase"
    
    # Refer to Python's regular expression at https://docs.python.org/2/library/re.html
    #   <td width="16%" align="center"><font size="4">47.2 </font> </td>
    allele_re = re.compile('>(\d+\.?\d?\"?\'*\(?\d*\.?\d?\"?\'*\)?\*?)</')
    #   <td width="35%"><font size="2">[TTTC]<sub>4</sub>TTTT TT<span style="mso-spacerun: yes"> </span>[CTTT]<sub>14</sub>[CTGT]<sub>3</sub>[CTTT]<sub>14 </sub>[CTTC]<sub>4</sub>[CTTT]<sub>3</sub>CTCC[TTCC]<sub>4</sub></font> </td>
    # repeat_re = re.compile('^(\[[ACGT]+\]\d+|[ACGT]+)+$')
    repeat_re = re.compile('^(\[[ACGT]+\]\d+|[ACGT]+|\s)+$')
    # Remove extra tags
    tag_re = re.compile('(<[^>]*>)')
    nbsp_re = re.compile('&nbsp;')
    quot_re = re.compile('&quot;')
    codis_data_file = open(base_fname + ".dat", 'w')
    for locus_name in CODIS_loci:
        if len(locus_list) > 0 and locus_name not in locus_list:
            continue
        url = "%s/str_%s.htm" % (base_url, locus_name)
        content = get_html(url).split("\r\n")
        content = map(lambda x: x.strip(), content)
        content2 = []
        for line in content:
            if line.startswith("<t") or \
               line.startswith("</tr") or \
               len(content2) == 0:
                content2.append(line)
            else:
                content2[-1] += line
                
        content = content2
        alleles = []
        l = 0
        while l < len(content):
            line = content[l]
            if line.startswith("<tr"):
                l += 1
                if l < len(content):
                    line = content[l]
                    line = re.sub(nbsp_re, '', line)
                    line = re.sub(quot_re, "''", line)
                    line = line.replace(' ', '')
                    allele_match = allele_re.search(line)
                    if not allele_match:
                        continue
                    allele_id = allele_match.group(1)                        
                    l += 1
                    repeat_match = None
                    while l < len(content):
                        line = content[l]                        
                        if not line.startswith("<td"):
                            break
                        line = re.sub(tag_re, '', line)
                        line = re.sub(nbsp_re, '', line)
                        repeat_match = repeat_re.search(line)
                        if repeat_match:
                            break
                        l += 1
                        
                    if not repeat_match:
                        continue

                    repeat_st = line
                    alleles.append([allele_id, repeat_st])
            else:
                l += 1

        for allele_id, repeat_st in alleles:
            print >> codis_data_file, "%s\t%s\t%s" % (locus_name, allele_id, repeat_st)

    codis_data_file.close()


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Extract multiple sequence alignments for DNA Fingerprinting loci")
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="codis",
                        help="base filename (default: codis)")
    parser.add_argument("--locus-list",
                        dest="locus_list",
                        type=str,
                        default="",
                        help="base filename (default: empty)")    
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.base_fname.find('/') != -1:
        elems = args.base_fname.split('/')
        base_fname = elems[-1]
        base_dname = '/'.join(elems[:-1])
    else:
        base_fname = args.base_fname
        base_dname = ""
    if args.locus_list != "":
        locus_list = args.locus_list.split(',')
    else:
        locus_list = []
        
    download_codis(base_dname,
                   base_fname,
                   locus_list,
                   args.verbose)

