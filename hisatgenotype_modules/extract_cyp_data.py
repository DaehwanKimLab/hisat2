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
import inspect, operator
import glob
from argparse import ArgumentParser, FileType


global gene_names
gene_names = ['cyp1a1','cyp1a2','cyp1b1','cyp2a6',
              'cyp2a13','cyp2b6','cyp2c8','cyp2c9',
              'cyp2c19','cyp2d6','cyp2e1','cyp2f1',
              'cyp2j2','cyp2r1','cyp2S1','cyp2w1',
              'cyp3a4','cyp3a5','cyp3a7','cyp3a43',
              'cyp4a11','cyp4a22','cyp4b1','cyp4f2',
              'cyp5a1','cyp8a1','cyp19a1','cyp21a2',
              'cyp26a1']

"""
Download variant information from website database
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
    print("Downloading data from:")
    
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

    os.system('mkdir cyp_var_files')
    for cyp_url in cyp_urls:
        cyp_gene_name = cyp_url.split('/')[-1]
        cyp_gene_name = cyp_gene_name.split('.')[0]
        
        # Hardcoded for cyp21 database (has inconsistant url naming) 
        if cyp_gene_name.lower() == "cyp21".lower():
            cyp_gene_name = cyp_gene_name + "a2" 

        # Changed to match all instances of "cyp"
        if not re.compile("cyp[\d\w]+", re.IGNORECASE).search(cyp_gene_name):
            continue

        # Open file to write on
        cyp_file = open("cyp_var_files/%s.var" % (cyp_gene_name), 'w')
        
        print >> sys.stderr, cyp_url, cyp_gene_name
        print >> cyp_file, cyp_url, cyp_gene_name

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
            varInfo_re = re.compile('-?\d+[ACGT]\&gt;[ACGT]|-?\d+_?-?\d+?del[ACGT]+|-?\d+_?-?\d+?ins[ACGT]+|None')

            alleleName = allele_name_re.findall(tabRow[0])
            if len(alleleName) > 0:
                alleleName = alleleName[0]

            # @RaymonFix - some databases have extra table, ignores headers (CYP2A6)
            # @Daehwan - some databases (e.g. http://www.cypalleles.ki.se/cyp3a4.htm)
            #            have 2 rows of Nucleotide changes (cDNA and Gene), might need
            #            to look at all rows for snps
            #
            # @RaymonFix - look in 4th column for "Gene" nt changes first, then consider cDNA if applicable; updated re to remove "<>" formating expressions 

            if cyp_url == 'http://www.cypalleles.ki.se/cyp21.htm': # Hardcoded for special format for cyp21a2
                try:
                    varInfo = varInfo_re.findall(re.sub('<[^>]+>', '',tabRow[1]))
                except IndexError:
                    continue
                
            else:
                try:
                    varInfo = varInfo_re.findall(re.sub('<[^>]+>', '',tabRow[3]))
                    if len(varInfo) == 0:
                        varInfo = varInfo_re.findall(re.sub('<[^>]+>', '',tabRow[2]))
                except IndexError:
                    continue

            for varInd in range(len(varInfo)):
                varInfo[varInd] = varInfo[varInd].replace('&gt;','>')

            if 'None' in varInfo:
                try:
                    assert len(varInfo) == 1
                except:
                    varInfo = filter(lambda a: a != 'None', varInfo)
                
        
            if isinstance(alleleName, basestring):
                print >> cyp_file, (str(alleleName) + "\t" + ','.join(varInfo))
            
        cyp_file.close()

         
"""
Make MSF files from variants
"""

def checkNTloc(fasta_fileName,var_fileName,gene_name):
    print "\nGene: %s" % gene_name
    seq = ""
    for line in open(fasta_fileName,'r'):
        if line[0] == '>':
            continue
        seq += line.strip()

    cyp_var_file = open(var_fileName,'r')
    cyp_var_dict = makeVarDict(cyp_var_file)
    cyp_var_file.close()

    print "len:", len(seq)
    varsPos = set()
    varsNeg = set()

    for varList in cyp_var_dict.values():
        for var in varList:
            if ">" in var: # is SNP
                posNt = int(var[:-3])
                ntChange = var[-3:].replace('>','')
                assert len(ntChange) == 2
                for nt in ntChange:
                    assert nt in "ACGT"

                if posNt > 0:
                    varsPos.add(str(posNt) + '->' + ntChange[0])
                else:
                    assert posNt < 0
                    varsNeg.add(str(posNt) + '->' + ntChange[0])
                    
            elif "del" in var: # is deletion
                posNt = var.split('del')[0].split('_')
                posNt = [int(p) for p in posNt]
                ntDel = var.split('del')[1]
                for nt in ntDel:
                    assert nt in "ACGT"

                if len(posNt) == 1: # single nt deletion
                    assert len(ntDel) == 1
                    if posNt[0] > 0:
                        varsPos.add(str(posNt[0]) + '->' + ntDel)
                    else:
                        assert posNt[0] < 0
                        varsNeg.add(str(posNt[0]) + '->' + ntDel)

                else: # mutliple nt deletion
                    assert len(posNt) == 2
                    try:
                        assert posNt[1] - posNt[0] + 1 == len(ntDel)
                    except AssertionError:
                        print "Incorrect deletion format: %s , skipping variation" % (var)
                        '''sys.exit(1)'''
                        continue
                    ntDelList = list(ntDel)
                    for i in range(posNt[0],posNt[1] + 1):
                        if i > 0:
                            varsPos.add(str(i) + '->' + ntDelList.pop(0))
                        else:
                            assert i < 0
                            varsNeg.add(str(i) + '->' + ntDelList.pop(0))
                    assert len(ntDelList) == 0
                    
            else:
                assert ("ins" in var) or ("None" in var)
                continue
    
    scorePos = {} # { position offset : number of alignments } for positive positions
    for i in range(-len(seq), len(seq)):
        align_score = 0
        for var in varsPos:
            pos, base = var.split('->')
            pos = int(pos)
            
            try:
                seq[pos+i]
            except IndexError:
                continue
            
            if seq[pos+i] == base:
                align_score += 1

        scorePos[i] = align_score
    oSetPos = max(scorePos.iteritems(), key=operator.itemgetter(1))[0]
    print "Positive postitions offset: %d" % oSetPos
    print "Score: %d out of %d\n" % (scorePos[oSetPos], len(varsPos))
    

    print "Checking negative position offset: %d" % (oSetPos + 1)
    align_score = 0
    oSetNeg = oSetPos + 1
    for var in varsNeg:
        pos, base = var.split('->')
        pos = int(pos)
        
        try:
            seq[pos + oSetNeg]
        except IndexError:
            continue
        
        if seq[pos + oSetNeg] == base:
            align_score += 1
    print "Score: %d out of %d\n\n" % (align_score, len(varsNeg))

    if len(varsNeg) == 0 and len(varsPos) != 0:
        return oSetPos, oSetNeg, float(scorePos[oSetPos])/float(len(varsPos)), 1.0, float(scorePos[oSetPos] + align_score)/float(len(varsPos) + len(varsNeg))
    elif len(varsNeg) != 0 and len(varsPos) == 0:
        return oSetPos, oSetNeg, 1.0, float(align_score)/float(len(varsNeg)), float(scorePos[oSetPos] + align_score)/float(len(varsPos) + len(varsNeg))
    elif len(varsNeg) == 0 and len(varsPos) == 0:
        return oSetPos, oSetNeg, 1.0, 1.0, 1.0
    else:
        assert len(varsNeg) != 0 and len(varsPos) != 0
        return oSetPos, oSetNeg, float(scorePos[oSetPos])/float(len(varsPos)), float(align_score)/float(len(varsNeg)), float(scorePos[oSetPos] + align_score)/float(len(varsPos) + len(varsNeg))
        

def create_map(seq):
    seq_map = {}
    count = 0
    for i in range(len(seq)):
        bp = seq[i]
        if bp == '.':
            continue
        assert bp.upper() in "ACGT"
        seq_map[count] = i
        count += 1
    return seq_map

def splitString(someStr,posList):
    posList.insert(0,-1)
    posList.append(len(someStr) - 1)
    splitStr = []
    for i in range(len(posList) - 1):
        left = posList[i] + 1
        right = posList[i+1] + 1
        splitStr.append(someStr[left:right])

    return splitStr

def extractSeq(faFile):
    seq = ""
    for line in faFile:
        if line.startswith(">"):
            continue

        seq += line.strip()

    return seq

def makeVarDict(fname):
    alleleVarDict = {}

    allLines = [line.strip() for line in fname]

    ref_al_id_present = False
    for line in allLines[1:]:
        if 'None' in line:
            ref_al_id_present = True

    line_num = 0
    for line in allLines[1:]:
        line_num += 1
        assert line.upper().startswith("CYP")
        alleleName = line.split("\t")[0].upper()

        if (not ref_al_id_present) and line_num == 1:
            varList = ['None']            
        else:
            try:
                varList = line.split("\t")[1].split(',')
            except IndexError:
                continue
        
        try:
            assert not alleleName in alleleVarDict
            alleleVarDict[alleleName] = set(varList)
        except:
            print >> sys.stdout, ("Warning, %s allele is already represented" % alleleName)
            alleleVarDict[alleleName] = alleleVarDict[alleleName] | set(varList)

    return alleleVarDict

def makeSnp(oldSeq, pos, oldNt, newNt):
    assert oldSeq[pos] == oldNt
    newSeq = oldSeq[:pos] + newNt + oldSeq[pos+1:]
    assert len(newSeq) == len(oldSeq)
    return newSeq

def makeDel(oldSeq, left, right, toDel):
    assert right - left + 1 == len(toDel)
    assert oldSeq[left:right + 1] == toDel
    newSeq = oldSeq[:left] + '.'*len(toDel) + oldSeq[right + 1:]
    assert len(newSeq) == len(oldSeq)
    return newSeq
    
def makeIns(oldSeq,left,right,toIns):
    assert right - left - 1 >= len(toIns)
    for nt in oldSeq[left + 1:right]:
      assert nt == '.'
    remDots = right - left - 1 - len(toIns)
    newSeq = oldSeq[:left + 1] + toIns + '.'*remDots + oldSeq[right:]
    assert len(newSeq) == len(oldSeq)
    return newSeq
    

def makeMSF(gene_name, oSetPos, oSetNeg):
    cyp_var_file = open("cyp_var_files/%s.var" % gene_name,'r')
    cyp_var_dict = makeVarDict(cyp_var_file)
    cyp_var_file.close()

    if len(cyp_var_dict) < 2:
        print('\tOnly reference allele included, skipping gene')
        return

    try:
        blast_allele_var = extract_var_from_blast('cyp_blast_alignment/%s_blast.align' % gene_name)
        if len(blast_allele_var) > 0:
            cyp_var_dict[gene_name.upper() + '*REFGRCH38P7'] = set(blast_allele_var)
    except IOError:
        print('\t%s blast file was skipped.' % gene_name)

    cyp_faFile = open("cyp_fasta/%s.fasta" % gene_name,'r')
    cyp_seq = extractSeq(cyp_faFile)
    cyp_faFile.close()
    preBackbone_seq = ''
    

    msfTable = {}

    # Building backbone structure (augment length with insertions)
    longestIns = {} # { key = position : value = length }
    for allele,varList in cyp_var_dict.items():
        for var in varList:
            if not "ins" in var:
                continue
            pos = var.split('ins')[0].split('_')
            pos = [int(p) for p in pos]
            ntIns = var.split('ins')[1]
            correctFormat = len(pos) == 2 and pos[1] - pos[0] == 1
            if not correctFormat:
                correctFormat = len(pos) == 1
            try:
                assert correctFormat
            except:
                print >> sys.stdout, "\tIncorrect format for insertion: variation %s on allele %s" % (var, allele)
                continue

            # convert to position in string
            if not 'GRCH38' in allele:
                if pos[0] > 0:
                    pos = pos[0] + oSetPos
                else:
                    pos = pos[0] + oSetNeg
            else:
                pos = pos[0]
                
            # Make dictionary of longest insertions
            if not pos in longestIns:
                longestIns[pos] = len(ntIns)
            else:
                if len(ntIns) > longestIns[pos]:
                    longestIns[pos] = len(ntIns)
    
    posInsList = sorted(longestIns.keys())
    
    splitSeq = splitString(cyp_seq,posInsList)
    posInsList = posInsList[1:-1]

    for i in range(len(posInsList)):
        splitSeq[i] += '.' * longestIns[posInsList[i]]

    for subseq in splitSeq:
        try:
            assert len(subseq) > 0 and not subseq.startswith('.')
            preBackbone_seq += subseq
        except:
            continue
    # pre-backbone built


    map_cyp = create_map(preBackbone_seq) # { Index of bp in original seq : Actual index in string }
    

    for allele,varList in cyp_var_dict.items():
        for var in varList:
            isSnp = False
            isDel = False
            isIns = False
        
            if ">" in var:
                isSnp = True
            elif "del" in var:
                isDel = True
            elif "ins" in var:
                isIns = True
            else:
                assert("None" in var)
                isRef = True

            if isSnp:
                pos = int(var[:-3])
                dbPos = pos
                ntChange = var[-3:].replace('>','')
                assert len(ntChange) == 2
                for nt in ntChange:
                    assert nt in "ACGT"

                if not 'GRCH38' in allele:
                    if pos > 0:
                        pos = pos + oSetPos
                    else:
                        pos = pos + oSetNeg

                if pos < 0 or pos > len(cyp_seq) - 1:
                    print >> sys.stdout, "\tWarning: position %d out of bounds" % (dbPos)
                    print >> sys.stdout, "\t\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                    continue
                    
                try:
                    assert(preBackbone_seq[map_cyp[pos]] == ntChange[0]) # nt at pos in seq must match database
                except:
                    print >> sys.stdout, "\tWarning: position %d in sequence contains %s, but expected %s from database" % (dbPos, preBackbone_seq[map_cyp[pos]], ntChange[0])
                    print >> sys.stdout, "\t\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                    continue
                
                # Adding to msf table
                if not allele in msfTable:
                    msfTable[allele] = makeSnp(preBackbone_seq, map_cyp[pos], ntChange[0], ntChange[1])
                else:
                    msfTable[allele] = makeSnp(msfTable[allele], map_cyp[pos], ntChange[0], ntChange[1])
                    
            elif isDel:
                pos = var.split('del')[0].split('_')
                pos = [int(p) for p in pos]
                if len(pos) == 1: # Handle single deletion with format for multi deletion with one location (e.g. [1707] -> [1707,1707])  
                    pos.append(pos[0])
                assert len(pos) == 2
                dbPos = pos
                ntDel = var.split('del')[1]
                for nt in ntDel:
                    assert nt in "ACGT"

                if not 'GRCH38' in allele:
                    for i in range(len(pos)):
                        if pos[i] > 0:
                            pos[i] = pos[i] + oSetPos
                        else:
                            pos[i] = pos[i] + oSetNeg

                skipDel = False
                for i in range(len(pos)):
                    if pos[i] < 0 or pos[i] > len(cyp_seq) - 1:
                        print >> sys.stdout, "\tWarning: position %d out of bounds" % (dbPos[i])
                        print >> sys.stdout, "\t\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                        skipDel = True

                if skipDel:
                    continue
                        
            
                try:
                    assert pos[1] - pos[0] + 1 == len(ntDel)
                except:
                    print >> sys.stdout, "\tIncorrect deletion data with %s on allele %s. Skipping variation." % (var, allele)
                    continue
                            
                try:
                    assert preBackbone_seq[ map_cyp[pos[0]] : map_cyp[pos[1]] + 1 ] == ntDel
                except:
                    print >> sys.stdout, "\tWarning, positions %d to %d in sequence contains %s, but expected %s from database" % \
                          (dbPos[0], dbPos[1], preBackbone_seq[ map_cyp[pos[0]] : map_cyp[pos[1]] + 1 ], ntDel)
                    print >> sys.stdout, "\t\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                    continue


                # Adding to msf table
                if not allele in msfTable:
                    msfTable[allele] = makeDel(preBackbone_seq, map_cyp[pos[0]], map_cyp[pos[1]], ntDel)
                else:
                    msfTable[allele] = makeDel(msfTable[allele], map_cyp[pos[0]], map_cyp[pos[1]], ntDel)

                        
            elif isIns:
                pos = var.split('ins')[0].split('_')
                pos = [int(p) for p in pos]
                if len(pos) == 1:
                    pos.append(pos[0] + 1)
                assert len(pos) == 2
                dbPos = pos
                try:
                    assert pos[1] - pos[0] == 1
                except AssertionError:
                    print >> sys.stdout, "\tIncorrect insertion data with %s on allele %s. Skipping variation." % (var, allele)
                    continue 
                ntIns = var.split('ins')[1]
                for nt in ntIns:
                    assert nt in "ACGT"

                if not 'GRCH38' in allele:
                    for i in range(len(pos)):
                        if pos[i] > 0:
                            pos[i] = pos[i] + oSetPos
                        else:
                            pos[i] = pos[i] + oSetNeg

                skipIns = False
                for i in range(len(pos)):
                    if pos[i] < 0 or pos[i] > len(cyp_seq) - 1:
                        print >> sys.stdout, "Warning: position %d out of bounds" % (dbPos[i])
                        print >> sys.stdout, "\tError occured on variation %s on allele %s. Skipping variation." % (var, allele)
                        skipIns = True

                if skipIns:
                    continue


                # Adding to msf table
                if not allele in msfTable:
                    msfTable[allele] = makeIns(preBackbone_seq, map_cyp[pos[0]], map_cyp[pos[1]], ntIns)
                else:
                    msfTable[allele] = makeIns(msfTable[allele], map_cyp[pos[0]], map_cyp[pos[1]], ntIns)


            else:
                assert isRef
                assert not allele in msfTable
                msfTable[allele] = preBackbone_seq

    # Sanity checking
    seq_len = 0
    for allele, msf_seq in msfTable.items():
        if seq_len == 0:
            seq_len = len(msf_seq)
        else:
            assert seq_len == len(msf_seq)
    assert seq_len > 0

    # Follow MSF style of IMGT/HLA database
    msfFile = open('cyp_msf/%s_gen.msf' % gene_name[3:].upper(),'w')
    for i in range(0, seq_len, 50):
        for allele, msf_seq in msfTable.items():
            output = "%12s" % allele[3:].upper()
            for j in range(i, i+50, 10):
                if j >= seq_len:
                    break
                if j == i:
                    output += "\t"
                else:
                    output += " "
                output += msf_seq[j:j+10]
            print >> msfFile, output
        print >> msfFile

    msfFile.close()


def build_msf_files():
    os.system('mkdir cyp_msf')

    oSetPos = 0
    oSetNeg = 0
    oSetScorePos = 0.0
    oSetScoreNeg = 0.0
    tot_score = 0.0
        
    print('\nBuilding MSF files:')
    for gene_name in gene_names:
        oSetPos, oSetNeg, oSetScorePos, oSetScoreNeg, tot_score = checkNTloc("cyp_fasta/%s.fasta" % gene_name,"cyp_var_files/%s.var" % gene_name,gene_name)
        if not (tot_score >= 0.95):
            print "\tLess than 95% match, skipping gene."
            continue
        
        makeMSF(gene_name, oSetPos, oSetNeg)


'''
Check MSF files against variants files
'''

global incorrect_msf_entries
incorrect_msf_entries = []

def create_inv_map(seq):
    seq_map = {}
    count = 0
    for i in range(len(seq)):
        bp = seq[i]
        if bp == '.':
            continue
        assert bp.upper() in "ACGT"
        seq_map[i] = count
        count += 1
    return seq_map

def readMSF(msf_fname): # { Allele name : MSF sequence }
    msf_dict = {}
    all_lines = [line for line in msf_fname]
    for line in all_lines:
        line = line.strip().replace(' ','')
        if len(line) == 0 : continue
        allele_name = 'CYP' + line.split('\t')[0]
        msf_seq = line.split('\t')[1]
        if not allele_name in msf_dict:
            msf_dict[allele_name] = msf_seq
        else:
            msf_dict[allele_name] = msf_dict[allele_name] + msf_seq

    return msf_dict

def msf_removeIns(ref_seq, al_seq):
    assert len(ref_seq) == len(al_seq)
    ins_ind_list = []
    for i in range(len(ref_seq)):
        if ref_seq[i] == '.':
            ins_ind_list.append(i)

    ori_ref_seq = ref_seq.replace('.','')
    ori_al_seq = list(al_seq)

    for i in ins_ind_list:
        ori_al_seq[i] = '-'

    ori_al_seq = ''.join(ori_al_seq).replace('-','')

    assert len(ori_ref_seq) == len(ori_al_seq)
    return ori_ref_seq, ori_al_seq

def msfToVarList(ref_seq, al_seq):
    var_list = []
    
    assert len(ref_seq) == len(al_seq)
    for bp in ref_seq: assert bp in "ACGT."
    for bp in al_seq: assert bp in "ACGT."
    inv_map = create_inv_map(ref_seq)
    
    ins_re = re.compile('[ACGT]\.+')
    ins_subStrPos = [(m.start(0), m.end(0)) for m in re.finditer(ins_re, ref_seq)] # list of duples of indicies of insertions in ref_seq
    ins_pos_length = [(tup[0], tup[1] - tup[0] - 1) for tup in ins_subStrPos]

    for tup in ins_pos_length:
        ins_pos, ins_length = tup[0], tup[1]
        ins_seq = al_seq[ins_pos + 1: ins_pos + ins_length  + 1]
        ins_seq = ins_seq.replace('.','')
        if len(ins_seq) == 0:
            continue
        ins_str_data = str(inv_map[tup[0]]) + '_' + str(inv_map[tup[0]] + 1) + 'ins' + ins_seq
        var_list.append(ins_str_data)

    # insertions finished
    
    ori_ref_seq, ori_al_seq = msf_removeIns(ref_seq, al_seq)

    for i in range(len(ori_ref_seq)):
        if ori_al_seq[i] == '.':
            continue 
        elif ori_al_seq[i] != ori_ref_seq[i]: # snp
            var_list.append(str(i) + ori_ref_seq[i] + '>' + ori_al_seq[i])

    del_subStrPos = [(m.start(0), m.end(0)) for m in re.finditer(ins_re, ori_al_seq)] # list of duples of indicies of deletions in ori_al_seq
    del_pos_length = [(tup[0], tup[1] - tup[0] - 1) for tup in del_subStrPos]

    for tup in del_pos_length:
        del_pos, del_length = tup[0], tup[1]
        del_seq = ori_ref_seq[del_pos + 1 : del_pos + del_length + 1]
        if del_length == 1:
            assert len(del_seq) == 1
            del_str_data = str(tup[0] + 1) + 'del' + del_seq
        else:
            del_str_data = str(tup[0] + 1) + '_' + str(tup[0] + tup[1]) + 'del' + del_seq
        var_list.append(del_str_data)

    # deletions finished

    return var_list

def checkMSFfile(gene_name, msf_fname, var_fname, fasta_filename):
    oSetPos, oSetNeg, oSet_pos_score, oSet_neg_score, tot_score = checkNTloc(fasta_filename, var_fname, gene_name)
    
    try:
        msf_file = open(msf_fname,'r')
        msf_dict = readMSF(msf_file) # { Allele name : MSF sequence }
        msf_file.close()
    except IOError:
        print("\t%s msf file was skipped.\n" % (gene_name))
        return

    var_file = open(var_fname,'r')
    var_dict = makeVarDict(var_file)
    var_file.close()

    try:
        blast_allele_var = extract_var_from_blast('cyp_blast_alignment/%s_blast.align' % gene_name)
        if len(blast_allele_var) > 0:
            var_dict[gene_name.upper() + '*REFGRCH38P7'] = set(blast_allele_var)
    except IOError:
        print('\t%s blast file was skipped.' % gene_name)
    
    fa_file = open(fasta_filename,'r')
    oriSeq = extractSeq(fa_file)
    fa_file.close()


    # Find reference allele
    ref_allele = ''
    for allele_name in var_dict.keys():
        if len(var_dict[allele_name]) == 1 and list(var_dict[allele_name])[0] == "None":
            assert ref_allele == ''
            ref_allele = allele_name
    assert not ref_allele == ''


    # Check if ref allele seq in msf matches fasta
    assert ref_allele in msf_dict

    try:
        assert msf_dict[ref_allele].replace('.','') == oriSeq
        print("Sequences match for reference allele %s" % ref_allele)
    except AssertionError:
        print("Warning: sequences do not match for reference allele %s" % ref_allele)
        sys.exit(1)


    # Check all alleles are included
    try:
        assert set([k.upper() for k in msf_dict.keys()]).issubset(set([k.upper() for k in var_dict.keys()]))
    except AssertionError:
        print("Extra alleles in MSF!\n")
        print(sorted(msf_dict.keys()))
        print("\n\n")
        print(sorted(var_dict.keys()))
        sys.exit(1)


    # Convert from database positions to sequence positions (using offset)
    for allele, var_list in var_dict.items():
        oSet_var_list = []
        for var in var_list:
            if '>' in var: # snp
                pos = int(var.split('>')[0][:-1])
                ntSnp = [var.split('>')[0][-1]]
                ntSnp.append(var.split('>')[1])
                assert len(ntSnp) == 2
                if not 'GRCH38' in allele:
                    if pos > 0:
                        pos = pos + oSetPos
                    else:
                        pos = pos + oSetNeg

                if pos < 0 or pos > len(oriSeq) - 1: # out of bounds
                    continue
                if oriSeq[pos] != ntSnp[0]: # mismatch
                    print('\tMismatch on variation %s' % var)
                    continue

                oSet_var = str(pos) + ntSnp[0] + '>' + ntSnp[1]
                oSet_var_list.append(oSet_var)

            elif 'del' in var: # deletion
                pos = var.split('del')[0].split('_')
                pos = [int(p) for p in pos]
                if len(pos) == 1: # Handle single deletion with format for multi deletion with one location (e.g. [1707] -> [1707,1707])  
                    pos.append(pos[0])
                assert len(pos) == 2
                ntDel = var.split('del')[1]
                for nt in ntDel:
                    assert nt in "ACGT"

                skipDel = False
                if not 'GRCH38' in allele:
                    for i in range(len(pos)):
                        if pos[i] > 0:
                            pos[i] = pos[i] + oSetPos
                        else:
                            pos[i] = pos[i] + oSetNeg
                        if pos[i] < 0 or pos[i] > len(oriSeq) - 1: # out of bounds
                            skipDel = True
                if (oriSeq[ pos[0] : pos[1] + 1 ] != ntDel): # mismatch
                    print('\tMismatch on variation %s' % var)
                    continue

                if skipDel:
                    continue

                assert pos[1] - pos[0] + 1 == len(ntDel)

                oSet_var = 'del' + ntDel
                if pos[0] == pos[1]:
                    oSet_var = str(pos[0]) + oSet_var
                else:
                    oSet_var = str(pos[0]) + '_' + str(pos[1]) + oSet_var

                oSet_var_list.append(oSet_var)                        

            elif 'ins' in var: # insertion
                pos = var.split('ins')[0].split('_')
                pos = [int(p) for p in pos]
                if len(pos) == 1:
                    pos.append(pos[0] + 1)
                assert len(pos) == 2
                try:
                    assert pos[1] - pos[0] == 1
                except AssertionError:
                    print('\tIncorrect insertion format on variation %s' % var)
                    continue
                ntIns = var.split('ins')[1]
                for nt in ntIns:
                    assert nt in "ACGT"

                skipIns = False
                if not 'GRCH38' in allele:
                    for i in range(len(pos)):
                        if pos[i] > 0:
                            pos[i] = pos[i] + oSetPos
                        else:
                            pos[i] = pos[i] + oSetNeg
                        if pos[i] < 0 or pos[i] > len(oriSeq) - 1: # out of bounds
                            skipIns = True

                if skipIns:
                    continue

                oSet_var = str(pos[0]) + '_' + str(pos[1]) + 'ins' + ntIns
                oSet_var_list.append(oSet_var)

            else:
                assert allele == ref_allele
                assert var == 'None'
                assert len(oSet_var_list) == 0
                oSet_var_list.append('None')

        var_dict[allele] = set(oSet_var_list)

    # Check variants created from MSF file against variants list
    num_correct_alleles = 0
    for allele, msf_seq in msf_dict.items():
        if allele == ref_allele:
            num_correct_alleles += 1
            continue
        msf_var_list = msfToVarList(msf_dict[ref_allele], msf_seq)
        '''print('\t' + str(var_dict[allele] == set(msf_var_list)) + '\t' + str(allele) + '\t' + str(msf_var_list))'''

        try:
            assert var_dict[allele] == set(msf_var_list)
            num_correct_alleles += 1
        except AssertionError:
            incorrect_msf_entries.append(allele)
            print('\n')
            print('\t\tVar File:\t' + str(var_dict[allele]))
            print('\t\tMSF File:\t' + str(set(msf_var_list)))
            print('\t\tDifference:\t' + str(var_dict[allele] - set(msf_var_list)) + '\n')
            '''sys.exit(1)'''

    print("\t%d out of %d alleles have correct msf sequences\n" % (num_correct_alleles, len(msf_dict)))

def check_msf_files():
    print("\nChecking MSF files:")

    for gene_name in gene_names:
        checkMSFfile(gene_name, 'cyp_msf/%s_gen.msf' % gene_name[3:].upper(), 'cyp_var_files/%s.var' % gene_name, 'cyp_fasta/%s.fasta' % gene_name)

    print('\n\n%d incorrect msf entries on alleles %s\n' % (len(incorrect_msf_entries), str(incorrect_msf_entries)))


"""
Write allele sequences to fasta for each gene
"""

def writeGenFasta(gene_name, msf_fname, line_length):
    try:
        msf_file = open(msf_fname,'r')
        msf_seq_dict = readMSF(msf_file)
        msf_file.close()
    except IOError:
        print("\t%s msf file was skipped." % (gene_name))
        return

    gen_fasta_file = open('gen_fasta/%s_gen.fasta' % gene_name[3:].upper(), 'w')
    
    for allele, seq in msf_seq_dict.items():
        seq = seq.replace('.','')
        print >> gen_fasta_file, ('>' + allele[3:].upper() + ' ' + str(len(seq)) + ' bp')
        seq_lines = [seq[i:i+line_length] for i in range(0, len(seq), line_length)]
        print >> gen_fasta_file, ('\n'.join(seq_lines))

    gen_fasta_file.close()
    print('%s_gen.fasta completed' % gene_name)

def build_gen_fasta_files():
    os.system('mkdir gen_fasta')

    print("\nBuilding alleles sequence fasta files:")
    for gene_name in gene_names:
        writeGenFasta(gene_name, 'cyp_msf/%s_gen.msf' % gene_name[3:].upper(), 60)


"""
Run script
"""

def extract_cyp_data():
    download_CYP(True)
    build_msf_files()
    check_msf_files()
    build_gen_fasta_files()

####################################################################################################
## Debuging BLASTN alignment ref alleles

def adjust_blast_vars(blast_vars_list,qry_pos):
    if len(blast_vars_list) == 0:
        return []

    qry_pos = qry_pos - 1
    adj_blst_var_list = []

    for var in blast_vars_list:
        if '>' in var: # SNP
            old_pos = int(var[:-3])
            adj_var = str(old_pos + qry_pos) + var[-3:]
            adj_blst_var_list.append(adj_var)
        elif 'del' in var: # deletion
            old_pos = var.split('del')[0].split('_')
            old_pos = [int(i) for i in old_pos]
            old_pos = [i + qry_pos for i in old_pos]
            if len(old_pos) == 1:
                adj_var = str(old_pos[0]) + 'del' + var.split('del')[1]
            else:
                assert len(old_pos) == 2
                adj_var = str(old_pos[0]) + '_' + str(old_pos[1]) + 'del' + var.split('del')[1]
            adj_blst_var_list.append(adj_var)
        else: # insertion
            assert 'ins' in var
            old_pos = var.split('ins')[0].split('_')
            old_pos = [int(i) for i in old_pos]
            old_pos = [i + qry_pos for i in old_pos]
            assert len(old_pos) == 2 and (old_pos[1] - old_pos[0] == 1)
            adj_var = str(old_pos[0]) + '_' + str(old_pos[1]) + 'ins' + var.split('ins')[1]
            adj_blst_var_list.append(adj_var)

    return adj_blst_var_list

def extract_var_from_blast(cyp_blast_fname):
    blastn_file = open(cyp_blast_fname,'r')
    all_lines = [line.strip() for line in blastn_file if not (len(line.strip()) == 0 or line.strip().startswith('|'))]
    blastn_file.close()

    id_match = [m.group(0) for l in all_lines[0:25] for m in [re.compile('.*(Identities.*).*').search(l)] if m][0]
    id_match = id_match.split('%')[0].split(' (')[0].split('= ')[1].split('/')
    id_match = [int(i) for i in id_match]

    # print(id_match)    
    assert len(id_match) == 2 and id_match[1] - id_match[0] >= 0
    if id_match[1] - id_match[0] == 0:
        print('\tPerfect match using blastn')
        return []
    
    
    start = -1
    end = -1
    for i in range(len(all_lines)): # Get rid of headers and footers
        if all_lines[i].startswith('Score ='):
            assert start == -1
            start = i

        if all_lines[i].startswith('Lambda'):
            assert start != -1 and end == -1
            end = i
            break

    all_lines = all_lines[start + 3 : end]
    # print('\n'.join(all_lines))

    blastn_var_list = []
    for i in range(0,len(all_lines),2):
        qry_seq = '\t'.join(all_lines[i].split())
        qry_seq_pos = int(qry_seq.split('\t')[1])
        sbj_seq = '\t'.join(all_lines[i + 1].split())
        qry_seq = qry_seq.split('\t')[2].replace('-','.').upper()
        sbj_seq = sbj_seq.split('\t')[2].replace('-','.').upper()
        #print(qry_seq)
        #print(sbj_seq)

        temp_var_list = msfToVarList(qry_seq, sbj_seq)
        #print(str(qry_seq_pos) + '\t' + str(temp_var_list) +  '\t' + str(adjust_blast_vars(temp_var_list,qry_seq_pos)))
        temp_var_list = adjust_blast_vars(temp_var_list,qry_seq_pos)
        blastn_var_list = blastn_var_list + temp_var_list
        
    return blastn_var_list

# extract_var_from_blast('cyp_blast_alignment/cyp2d6_blast.align')

extract_cyp_data()
