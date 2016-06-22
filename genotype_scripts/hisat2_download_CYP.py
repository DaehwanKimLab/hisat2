#!/usr/bin/env python

#
# Copyright 2016, Raymon Cao <rcao5@jhu.edu> and Daehwan Kim <infphilo@gmail.com>
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


def download_CYP(verbose):
    # CYP database base URL
    base_url = "http://www.cypalleles.ki.se"
    
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(download_CYP))
    ex_path = os.path.dirname(curr_script)

    # Refer to Python's regular expression at https://docs.python.org/2/library/re.html
    cyp_re = re.compile('http://www.cypalleles.ki.se/cyp\w+.htm')
    output = get_html(base_url)
    cyp_urls = cyp_re.findall(output)
    # Original list had duplicate urls, removes duplicates
    cyp_urls = set(cyp_urls)

    cyp_file = open("cyp.web.output", 'w')    
    for cyp_url in cyp_urls:
        cyp_gene_name = cyp_url.split('/')[-1]
        cyp_gene_name = cyp_gene_name.split('.')[0]
        
        # Hardcoded for cyp21 database (has inconsistant url naming) 
        if cyp_gene_name.lower() == "cyp21".lower():
            cyp_gene_name = cyp_gene_name + "a2" 

        # Changed to match all instances of "cyp"
        if not re.compile("cyp[\d\w]+", re.IGNORECASE).search(cyp_gene_name):
            continue

        # Raymon - future work - CYP2A6        
        print >> sys.stderr, "\n\n", cyp_url, cyp_gene_name
        print >> cyp_file, "\n\n", cyp_url, cyp_gene_name

        cyp_output = get_html(cyp_url)
        if cyp_output == "":
            continue

        listA = cyp_output.split("<tr style=")

        indStart = -1
        foundStart = False
        while not foundStart:
            indStart += 1
            foundStart = (cyp_gene_name + '*').upper() in listA[indStart].upper()
            
        # Look for first occurance of "[cyp_gene_name]*"
        listA = listA[indStart:]

        # Look for last occurance of "[cyp_gene_name]*"
        indEnd = 0
        foundEnd = False
        while not foundEnd:
            indEnd -= 1
            foundEnd = (cyp_gene_name + '*').upper()  in listA[indEnd].upper()

        listA = listA[:(indEnd + 1)]
        
        for itemA in listA:
            tabRow = itemA.split("</td>")
            for ind in range(len(tabRow)):
                tabRow[ind] = tabRow[ind].replace("\r\n","")

            allele_name_re = re.compile(cyp_gene_name.upper() + '\*[\w\d]+')
            varInfo_re = re.compile('-?\d+[ACGT]\&gt;[ACGT]|-?\d+_?-?\d+?del[ACGT]+|-?\d+_?-?\d+?ins[ACGT]+')

            alleleName = allele_name_re.findall(tabRow[0])
            if len(alleleName) > 0:
                alleleName = alleleName[0]

            # @RaymonFix - some databases have extra table, ignores headers (CYP2A6)
            # @Daehwan - some databases (e.g. http://www.cypalleles.ki.se/cyp3a4.htm)
            #            have 2 rows of Nucleotide changes (cDNA and Gene), might need
            #            to look at all rows for snps
            try:
                varInfo = varInfo_re.findall(tabRow[2])
            except IndexError:
                continue

            for varInd in range(len(varInfo)):
                varInfo[varInd] = varInfo[varInd].replace('&gt;','>')
                
        
            if isinstance(alleleName, basestring):
                print >> cyp_file, (str(alleleName) + "\t" + ','.join(varInfo))
            
    cyp_file.close()
            
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Download CYP gene information")
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="hla",
                        help="base filename")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    download_CYP(args.verbose)
