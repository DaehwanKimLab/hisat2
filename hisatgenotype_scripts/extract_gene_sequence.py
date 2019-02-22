import sys, os, subprocess, re
import urllib2
import xml.etree.ElementTree as etree
from multiprocessing import Pool
from Bio import Entrez, SeqIO
import time
from argparse import ArgumentParser, FileType
import shutil
import copy

def get_xml(url):
    file = urllib2.urlopen(url)
    data = file.read()
    file.close()

    tree = etree.fromstring(data)
    return tree

def get_RefSeqID(locus = {}, verbose = False):
    refseq = {}
    for gene in locus:
        refseq.update({ gene : [] })

    webaddress = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene'
    try:
        website = get_website(webaddress, verbose)
    except ValueError:
        if verbose: print >> sys.stderr, 'Cannot access refseq database at this time'
        raise ValueError()

    for line in website:
        line = line.strip().split('\t')
        if line[2] in refseq:
            refseq[line[2]].append(line[3])
      
    return refseq

def get_website(url, verbose = False):
    tries = 0
    while tries < 4:
        tries += 1
        try:
            if verbose: print >> sys.stderr, 'Attempt %d of 3 to connect to URL' % tries
            resp = urllib2.urlopen(url)
        except IOError:
            if verbose: print >> sys.stderr, 'No Connection to URL'
            if tries == 3:
                raise ValueError()
        else:
            if verbose: print >> sys.stderr, 'Connection successful!'
            file = resp.read()
            return file.splitlines(True)

def build_symbol_uid_dic(verbose = False):
    webaddress = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene'
    try:
        website = get_website(webaddress, verbose)
    except ValueError:
        if verbose: print >> sys.stderr, 'Cannot access refseq database at this time'
        raise ValueError()
    symbol_uid_dic = {}

    for line in website:
        line = line.strip().split('\t')
        uid, symbol = line[1:3]
        symbol_uid_dic[symbol] = uid
      
    return (symbol_uid_dic)

def symbol2accession(gene_symbol_list, dic):

    def get_ncbi_tree(uid):
        page_address = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s&rettype=xml" % str(uid)

        page = urllib2.Request(page_address)
        xml  = urllib2.urlopen(page)
        tree = etree.parse(xml)
        xml.close()

        return(tree)

    def get_offical_symbol(mother_tree):
        for commentary in mother_tree.iter('Gene-commentary'):
            for commentary_label in commentary.iter('Gene-commentary_label'):
                if commentary_label.text == 'Official Symbol' and commentary[2].tag == 'Gene-commentary_text':
                    symbol = commentary[2].text
                    
                    return(symbol)    

    def get_accession_tree(mother_tree):
        for commentary in mother_tree.iter('Gene-commentary'):
            for commentary_heading in commentary.iter('Gene-commentary_heading'):
                if commentary_heading.text == 'Related Sequences':
                    accession_tree = commentary[2]
                    
                    return(accession_tree)

    def get_accession_list(accession_tree):
        accession_list = []

        for i in range(len(accession_tree)):
            for child in accession_tree[i]:
                if child.tag == 'Gene-commentary_accession':
                    accession = child.text

                elif child.tag == 'Gene-commentary_version':
                    version = child.text

            full_accession = accession + '.' + version
            accession_list.append(full_accession)
        return(accession_list)
    
    accession_dict = {}
    
    for gene_symbol in gene_symbol_list:
        
        uid = dic[gene_symbol]
        tree = get_ncbi_tree(uid)

        #offical_symbol = get_offical_symbol(tree)
        accession_tree = get_accession_tree(tree)

        accession_list = get_accession_list(accession_tree)
        accession_list = list(set(accession_list))
        
        accession_dict[gene_symbol] = accession_list
        time.sleep(1)
    return(accession_dict)

def run_clustalo(in_file):
    out_file = in_file.replace('.fasta', '.msf')
    cmd = 'clustalo --in=%s --full --outfmt=msf --out=%s --threads=1 --force' % (in_file, out_file)
    os.system(cmd)

def run_clustalo_pile(in_file):
    out_file = in_file.replace('.fasta', '.msf')
    cmd = 'clustalo --in=%s --pileup --outfmt=msf --out=%s --threads=1 --force' % (in_file, out_file)
    # cmd = 'clustalw2 -infile=%s -output=GDE -outfile=%s -quiet -case=upper > /dev/null' % (in_file, out_file)
    os.system(cmd)

def get_range(tree):
    output_range = []
    for seq_interval in tree.getiterator('Seq-interval'):
        interval_from = seq_interval.getiterator('Seq-interval_from')[0].text
        interval_to   = seq_interval.getiterator('Seq-interval_to')[0].text
        output_range.append((int(interval_from), int(interval_to)))
    if not output_range:
        point = tree.getiterator('Seq-point_point')[0].text
        output_range.append((int(point), int(point) + 1))
    return(output_range)

def range_2_seq(whole_sequence, seq_ranges):
    seqs = []
    for seq_range in seq_ranges:
        seq_start, seq_end = seq_range
        seq = whole_sequence[seq_start: seq_end + 1]
        seqs.append(seq)
    return(seqs)

def range_2_exon_intron(whole_sequence, exon_range):
    exons_seq   = []
    introns_seq = []
    output_dict = {}
    range_label = []

    for i in range(len(exon_range)):
        exon_start, exon_end = exon_range[i]
        exon_seq = whole_sequence[exon_start: exon_end + 1]
        exons_seq.append(exon_seq)
        range_label.append([exon_start, exon_end, 'exon'])
        if i > 0:
            intron_start = exon_range[i - 1][-1] + 1
            intron_end   = exon_start - 1
            intron_seq   = whole_sequence[intron_start: intron_end + 1]
            introns_seq.append(intron_seq)
            range_label.append([intron_start, intron_end, 'intron'])

    output_dict['exons']   = exons_seq
    output_dict['introns'] = introns_seq
    range_label.sort()
    return(output_dict, range_label)

def range_2_UTR(whole_sequence, exon_range, boundry_range):
    output_dict = {}
    boundry_start, boundry_end = boundry_range
    exon_start   , exon_end    = exon_range[0][0], exon_range[-1][-1]
    UTR_5 = whole_sequence[boundry_start: exon_start]
    UTR_3 = whole_sequence[exon_end + 1: boundry_end + 1]
    output_dict['UTR_5'] = UTR_5
    output_dict['UTR_3'] = UTR_3
    return(output_dict)

def get_sequence(symbol_list, symbol_uid_dic, verbose = False):

    step = 32
    accession_lists = symbol2accession(symbol_list, symbol_uid_dic)
    ref_access_list = get_RefSeqID(symbol_list)
    master_output_dict = {}
    ref_range_label    = {}
    ref_sequence       = {}

    for symbol in symbol_list:
        ref_access     = ref_access_list[symbol][0]
        accession_list = accession_lists[symbol] 
        accession_list.insert(0, ref_access)
        accession_list = list(set(accession_list))
        master_output_dict[symbol] = {}
        i = 0
        
        while i <= len(accession_list):
            time.sleep(3) 
            accession_sublist = accession_list[i:i + step]

            webaddress = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&dopt=gb&sendto=on&id=%s&rettype=xml' % ','.join(accession_sublist)
            try:
                page = urllib2.Request(webaddress)
            except ValueError:
                if verbose: print >> sys.stderr, 'No data available for accession: %s' % access
                raise ValueError()

            xml = urllib2.urlopen(page)
            master_tree = etree.parse(xml)
            xml.close()
            root = master_tree.getroot()

            handle = Entrez.efetch(db="nuccore", id=accession_sublist, rettype="fasta", retmode="text")
            records = handle.read().strip().split('>')
            handle.close()

            whole_sequence_dict = {}

            for record in records:
                if record:
                    whole_sequence = ''.join(record.split('\n')[1:])
                    access_id = record.split('\n')[0]
                    access_id = access_id.split(' ')[0]
                    whole_sequence_dict[access_id] = whole_sequence

            output_dict = {}
  
            seq_entrys = root.findall('./Bioseq-set_seq-set/Seq-entry')

            for access in accession_sublist:

                access_seq_dict = {}

                whole_sequence = whole_sequence_dict[access]

                if '.' in access:  
                    access_no_version = access.split('.')[0]
                else:
                    access_no_version = access


                for seq_entry in seq_entrys:
                    for accession_id in seq_entry.getiterator('Textseq-id_accession'):
                        if accession_id.text == access_no_version:
                            target_seq_entry = seq_entry

                cds_range        = []
                rna_range        = []
                gene_range       = []
                pseudo_cds_range = []
                exon_range       = []
                intron_range     = []
                UTR_5_range      = []
                UTR_3_range      = []
                access_seq_dict  = {}

                for bioseq in target_seq_entry.getiterator('Bioseq'):
                    check = bioseq.getiterator('Bioseq')
                    if len(check) != 1:
                        continue
                    if bioseq.getiterator('Textseq-id_accession')[0].text == access_no_version:
                        info = bioseq

                set_annot = target_seq_entry.find('./*/*/Bioseq-set_annot')
                if set_annot is not None:
                    for seq_feat in set_annot.getiterator('Seq-feat'):
                        for category in seq_feat.iter():
                            if category.tag == 'SeqFeatData_cdregion' or category.tag == 'SeqFeatData_rna':
                                cds_range = get_range(seq_feat.find('./Seq-feat_location')) 
                else:
                    cds_range = []

                for Seq_feat in info.getiterator('Seq-feat'):
                    Seq_feat_partial = Seq_feat.find('./Seq-feat_partial')
                    if Seq_feat_partial is not None:
                        continue

                    Seq_feat_data = Seq_feat.find('Seq-feat_data')
                    for category in Seq_feat_data.iter():
                        if category.tag == 'SeqFeatData_gene':
                            for ref_locus in Seq_feat.iter('Gene-ref_locus'):
                                if ref_locus.text == symbol:
                                    gene_range = get_range(Seq_feat)
                                    continue
                        elif category.tag == 'SeqFeatData_rna':
                            rna_range = get_range(Seq_feat)
                            continue
                        elif category.tag == 'SeqFeatData_cdregion':
                            cds_range = get_range(Seq_feat) 
                            continue
                        elif category.tag == 'Seq-feat_pseudo':
                            check = Seq_feat.getiterator('SeqFeatData_gene')
                            if not check:
                                pseudo_cds_range = get_range(Seq_feat)
                                continue
                        elif category.tag == 'SeqFeatData_imp':
                            for feat_key in Seq_feat.getiterator('Imp-feat_key'):
                                if feat_key.text == 'exon':
                                    exon_range.append(get_range(Seq_feat)[0])
                                    continue
                                elif feat_key.text == 'intron':
                                    intron_range.append(get_range(Seq_feat)[0])
                                    continue
                                elif feat_key.text == '3\'UTR':
                                    UTR_3_range.append(get_range(Seq_feat)[0])
                                    continue
                                elif feat_key.text == '5\'UTR':
                                    UTR_5_range.append(get_range(Seq_feat)[0])
                                    continue

                if gene_range:
                    boundry_range = gene_range[0]
                elif rna_range:
                    rna_range.sort()
                    boundry_range = (rna_range[0][0], rna_range[-1][-1])
                elif cds_range:
                    cds_range.sort()
                    boundry_range = (cds_range[0][0], cds_range[-1][-1])
                elif pseudo_cds_range:
                    pseudo_cds_range.sort()
                    boundry_range = (pseudo_cds_range[0][0], pseudo_cds_range[-1][-1])

                if cds_range:
                    cds_range.sort()
                    access_seq_dict, range_label = range_2_exon_intron(whole_sequence, cds_range)
                    UTR_dict = range_2_UTR(whole_sequence, cds_range, boundry_range)
                    label_range = cds_range
                elif rna_range:
                    rna_range.sort()
                    access_seq_dict, range_label = range_2_exon_intron(whole_sequence, rna_range)
                    UTR_dict = range_2_UTR(whole_sequence, rna_range, boundry_range)
                    label_range = rna_range
                elif pseudo_cds_range:
                    pseudo_cds_range.sort()
                    access_seq_dict, range_label = range_2_exon_intron(whole_sequence, pseudo_cds_range)
                    UTR_dict = range_2_UTR(whole_sequence, pseudo_cds_range, boundry_range)
                    label_range = pseudo_cds_range
                else:
                    UTR_dict = {}
                    access_seq_dict['exons']   = []
                    access_seq_dict['introns'] = []
                    access_seq_dict['UTR_5']   = []
                    access_seq_dict['UTR_3']   = []

                access_seq_dict.update(UTR_dict)

                if len(access_seq_dict['exons']) == 0 and exon_range:
                    exon_range.sort()
                    access_seq_dict['exons'] = range_2_seq(whole_sequence, exon_range)
                if len(access_seq_dict['introns']) == 0 and intron_range:
                    intron_range.sort()
                    access_seq_dict['introns'] = range_2_seq(whole_sequence, intron_range) 
                if len(access_seq_dict['UTR_5']) == 0 and UTR_5_range:
                    if len(label_range) > 0:
                        if UTR_5_range[-1] < label_range[0][0]:
                            access_seq_dict['UTR_5'] = range_2_seq(whole_sequence, UTR_5_range)[0]
                if len(access_seq_dict['UTR_3']) == 0 and UTR_3_range:
                    if len(label_range) > 0:
                        if UTR_3_range[0] > label_range[-1][-1]:
                            access_seq_dict['UTR_3'] = range_2_seq(whole_sequence, UTR_3_range)[0]
                
                if access == ref_access:
                    range_label.append([boundry_range[0],  label_range[0][0] - 1,   'UTR_5'])
                    range_label.append([label_range[-1][-1] + 1, boundry_range[-1], 'UTR_3'])
                    range_label.sort()
                    
                    ref_boundary = (range_label[0][0], range_label[-1][1])
                    range_label_dict = {}
                    for label_i in range(len(range_label)):
                        range_label_dict[label_i] = [(range_label[label_i][0] - ref_boundary[0],
                                                range_label[label_i][1] - ref_boundary[0]),
                                                range_label[label_i][2]]  
                    ref_range_label[symbol] = range_label_dict
                    
                    ref_sequence[symbol] = whole_sequence[ref_boundary[0]: ref_boundary[1] + 1]
                else:
                    output_dict[access] = access_seq_dict
                

            master_output_dict[symbol].update(output_dict)
            i = i + step
            
    return(master_output_dict, ref_range_label, ref_sequence)

def load_local_sequence(symbol_list, local_filename):
    whole_sequence_dict = {}
    for symbol in symbol_list:
        whole_sequence_dict[symbol] = {}
    
    f = open(local_filename)
    for line in f:
        symbol, access, sequence = line.strip().split('\t')
        if symbol == 'REM1':
            continue
        whole_sequence_dict[symbol][access] = sequence
        
    return(whole_sequence_dict)

def get_sequence_local(symbol_list, local_filename, verbose = False):
    
    whole_sequence_dict = load_local_sequence(symbol_list, local_filename)
    access_info_path = '/work/bioinformatics/s429377/genbank_test/'
    accession_lists = {}
    for symbol, access_list in whole_sequence_dict.items():
        accession_lists[symbol] = access_list.keys()
    
    ref_access_list = get_RefSeqID(symbol_list)
    master_output_dict = {}
    ref_range_label    = {}
    ref_sequence       = {}

    for symbol, access_list in accession_lists.items():
        
        ref_access = ref_access_list[symbol][0]
        master_output_dict[symbol] = {}
        
        for access in access_list:

            whole_sequence = whole_sequence_dict[symbol][access]

            if '.' in access:  
                access_no_version = access.split('.')[0]
            else:
                access_no_version = access

            xml_filename = access_info_path + access +'.xml'
            xml = open(xml_filename)
            tree = etree.parse(xml)
            xml.close()
            root = tree.getroot()
            
            cds_range        = []
            rna_range        = []
            gene_range       = []
            pseudo_cds_range = []
            exon_range       = []
            intron_range     = []
            UTR_5_range      = []
            UTR_3_range      = []
            access_seq_dict  = {}

            for bioseq in root.getiterator('Bioseq'):
                check = bioseq.getiterator('Bioseq')
                if len(check) != 1:
                    continue
                if bioseq.getiterator('Textseq-id_accession')[0].text == access_no_version:
                    info = bioseq

            set_annot = root.find('./*/*/Bioseq-set_annot')
            if set_annot is not None:
                for seq_feat in set_annot.getiterator('Seq-feat'):
                    for category in seq_feat.iter():
                        if category.tag == 'SeqFeatData_cdregion' or category.tag == 'SeqFeatData_rna':
                            cds_range = get_range(seq_feat.find('./Seq-feat_location')) 
            else:
                cds_range = []

            for Seq_feat in info.getiterator('Seq-feat'):
                Seq_feat_partial = Seq_feat.find('./Seq-feat_partial')
                if Seq_feat_partial is not None:
                    continue
                        
                Seq_feat_data = Seq_feat.find('Seq-feat_data')
                for category in Seq_feat_data.iter():
                    if category.tag == 'SeqFeatData_gene':
                        for ref_locus in Seq_feat.iter('Gene-ref_locus'):
                            if ref_locus.text == symbol:
                                gene_range = get_range(Seq_feat.find('./Seq-feat_location'))
                                continue
                    elif category.tag == 'SeqFeatData_rna':
                        rna_range = get_range(Seq_feat.find('./Seq-feat_location'))
                        continue
                    elif category.tag == 'SeqFeatData_cdregion':
                        cds_range = get_range(Seq_feat.find('./Seq-feat_location')) 
                        continue
                    elif category.tag == 'Seq-feat_pseudo':
                        check = Seq_feat.getiterator('SeqFeatData_gene')
                        if not check:
                            pseudo_cds_range = get_range(Seq_feat.find('./Seq-feat_location'))
                            continue
                    elif category.tag == 'SeqFeatData_imp':
                        for feat_key in Seq_feat.getiterator('Imp-feat_key'):
                            if feat_key.text == 'exon':
                                exon_range.append(get_range(Seq_feat.find('./Seq-feat_location'))[0])
                                continue
                            elif feat_key.text == 'intron':
                                intron_range.append(get_range(Seq_feat.find('./Seq-feat_location'))[0])
                                continue
                            elif feat_key.text == '3\'UTR':
                                UTR_3_range.append(get_range(Seq_feat.find('./Seq-feat_location'))[0])
                                continue
                            elif feat_key.text == '5\'UTR':
                                UTR_5_range.append(get_range(Seq_feat.find('./Seq-feat_location'))[0])
                                continue

            if gene_range:
                boundry_range = gene_range[0]
            elif rna_range:
                rna_range.sort()
                boundry_range = (rna_range[0][0], rna_range[-1][-1])
            elif cds_range:
                cds_range.sort()
                boundry_range = (cds_range[0][0], cds_range[-1][-1])
            elif pseudo_cds_range:
                pseudo_cds_range.sort()
                boundry_range = (pseudo_cds_range[0][0], pseudo_cds_range[-1][-1])

            if cds_range:
                cds_range.sort()
                access_seq_dict, range_label = range_2_exon_intron(whole_sequence, cds_range)
                UTR_dict = range_2_UTR(whole_sequence, cds_range, boundry_range)
                label_range = cds_range
            elif rna_range:
                rna_range.sort()
                access_seq_dict, range_label = range_2_exon_intron(whole_sequence, rna_range)
                UTR_dict = range_2_UTR(whole_sequence, rna_range, boundry_range)
                label_range = rna_range
            elif pseudo_cds_range:
                pseudo_cds_range.sort()
                access_seq_dict, range_label = range_2_exon_intron(whole_sequence, pseudo_cds_range)
                UTR_dict = range_2_UTR(whole_sequence, pseudo_cds_range, boundry_range)
                label_range = pseudo_cds_range
            else:
                UTR_dict = {}
                access_seq_dict['exons']   = []
                access_seq_dict['introns'] = []
                access_seq_dict['UTR_5']   = []
                access_seq_dict['UTR_3']   = []

            access_seq_dict.update(UTR_dict)

            if len(access_seq_dict['exons']) == 0 and exon_range:
                exon_range.sort()
                access_seq_dict['exons'] = range_2_seq(whole_sequence, exon_range)
            if len(access_seq_dict['introns']) == 0 and intron_range:
                intron_range.sort()
                access_seq_dict['introns'] = range_2_seq(whole_sequence, intron_range) 
            if len(access_seq_dict['UTR_5']) == 0 and UTR_5_range:
                if len(label_range) > 0:
                    if UTR_5_range[-1] < label_range[0][0]:
                        access_seq_dict['UTR_5'] = range_2_seq(whole_sequence, UTR_5_range)[0]
            if len(access_seq_dict['UTR_3']) == 0 and UTR_3_range:
                if len(label_range) > 0:
                    if UTR_3_range[0] > label_range[-1][-1]:
                        access_seq_dict['UTR_3'] = range_2_seq(whole_sequence, UTR_3_range)[0]

            if access == ref_access:
                range_label.append([boundry_range[0],  label_range[0][0] - 1,   'UTR_5'])
                range_label.append([label_range[-1][-1] + 1, boundry_range[-1], 'UTR_3'])
                range_label.sort()

                ref_boundary = (range_label[0][0], range_label[-1][1])
                range_label_dict = {}
                for label_i in range(len(range_label)):
                    range_label_dict[label_i] = [(range_label[label_i][0] - ref_boundary[0],
                                            range_label[label_i][1] - ref_boundary[0]),
                                            range_label[label_i][2]]  
                ref_range_label[symbol] = range_label_dict

                ref_sequence[symbol] = whole_sequence[ref_boundary[0]: ref_boundary[1] + 1]
            else:
                master_output_dict[symbol][access] = access_seq_dict
            
    return(master_output_dict, ref_range_label, ref_sequence)

def build_seqs_fasta(access_seqs, out_dir):
    seq_filenames = {}
    gene_list = access_seqs.keys()
    for symbol in gene_list:
        file_path = out_dir + '/%s/' % symbol
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        filename = file_path + '/%s_seqs.fasta' % symbol
        seq_filenames[symbol] = filename
        ofile = open(filename, 'w')
        for access, access_seq in access_seqs[symbol].items():
            for idx, exon_seq in enumerate(access_seq['exons']):
                if len(exon_seq) >= 25 and (exon_seq != len(exon_seq) * 'N'):
                    ofile.write('>' + access + '_exon_' + str(idx) + '\n')
                    ofile.write(exon_seq + '\n')
            for idx, intron_seq in enumerate(access_seq['introns']):
                if len(intron_seq) >= 25 and intron_seq != len(intron_seq) * 'N':
                    ofile.write('>' + access + '_intron_' + str(idx) + '\n')
                    ofile.write(intron_seq + '\n')
            if len(access_seq['UTR_5']) > 2 and access_seq['UTR_5'] != len(access_seq['UTR_5']) * 'N':
                ofile.write('>' + access + '_UTR_5' + '\n') 
                ofile.write(access_seq['UTR_5'] + '\n')
            if len(access_seq['UTR_3']) > 2 and access_seq['UTR_3'] != len(access_seq['UTR_3']) * 'N':
                ofile.write('>' + access + '_UTR_3' + '\n') 
                ofile.write(access_seq['UTR_3'] + '\n')
        ofile.close
        
    return(seq_filenames)

def build_ref_fasta(ref_sequence, out_dir):
    ref_filenames = {}
    for symbol, ref_seq in ref_sequence.items():
        file_path = out_dir + '/%s/' % symbol
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        filename = file_path + '/%s_ref.fasta' % symbol
        ref_filenames[symbol] = filename
        ofile = open(filename, 'w')
        ofile.write('>' + symbol + '_ref' + '\n')
        ofile.write(ref_seq + '\n')
        ofile.close
    return(ref_filenames)

def build_hisat2_index(ref_filenames):
    indexed_filnames = {}
    for symbol, ref_filename in ref_filenames.items():
        indexed_name = ref_filename.replace('.fasta', '_indexed')
        indexed_filnames[symbol] = indexed_name
        cmd_hisat2_build = ['hisat2-build', '-q', ref_filename, indexed_name]
        os.system(' '.join(cmd_hisat2_build))
    return(indexed_filnames)

def hisat2_alignment(seq_filenames, indexed_filenames):
    alignment_filenames = {}
    for symbol, indexed_filename in indexed_filenames.items():
        seq_filename = seq_filenames[symbol]
        alignment_filename = indexed_filename.replace('_ref_indexed', '_mapping_result')
        alignment_filenames[symbol] = alignment_filename
        cmd_hisat2_alignment = ['hisat2',
                                '-x',indexed_filename,
                                '-f',seq_filename,
                                '--no-spliced-alignment', '-S',
                                alignment_filename]
        os.system(' '.join(cmd_hisat2_alignment))
    return(alignment_filenames)

def assign_segment_correct_cigar(ref_ranges, symbol_exon_dict, mapping_filenames):

    symbol_seg_access_dict  = {}
    symbol_exon_access_dict  = {}
    seg_seq_access_dict  = {}
    back_mapped_path = {}

    
    for symbol in ref_ranges.keys():
        
        present_seg = []
        present_exon = []
        
        current_ref_range = copy.deepcopy(ref_ranges[symbol])
        file_name = mapping_filenames[symbol]
        f = open(file_name)

        seg_seq_access_dict[symbol]  = {}

        
        back_mapped_result_path = file_name.replace(symbol + '_mapping_result', 'umapped_mapping_result/')
        back_mapped_path[symbol] = back_mapped_result_path
        os.makedirs(back_mapped_result_path)
        
        ref_segments_file = file_name.replace('_mapping_result', '_ref_segments.fasta')
        
        #filename = '/home2/s429377/ABO/ABO_ref_segments.fasta'
        ref_segments_f = open(ref_segments_file, 'w')
        for seg_i, info in ref_ranges[symbol].items():
            seq = ref_seq[symbol][info[0][0]: info[0][1] + 1]
            ref_segments_f.write('>' + symbol + '_' + str(seg_i) + '\n')
            ref_segments_f.write(seq + '\n')
        ref_segments_f.close
        
        exons = []
        for seg, info in ref_ranges[symbol].items():
            seg_seq_access_dict[symbol][seg] = {}
            if info[1] == 'exon':
                exons.append(seg)

        for line in f:
            if line.startswith('@'):
                continue
        
            line = line.split('\t')

            if line[2] == '*':
                folder_name = file_name.replace('mapping_result', line[0])
                os.makedirs(folder_name)
                indexed_filename = folder_name + '/' + line[0] + '_indexed'
                cmd_hisat2_build = ['hisat2-build', '-q', '-c', line[9], indexed_filename]
                os.system(' '.join(cmd_hisat2_build))
                
                mapping_result_filename = back_mapped_result_path + line[0]
                cmd_hisat2_alignment = ['hisat2',
                                        '-x',indexed_filename,
                                        '-f',ref_segments_file,
                                        '--no-spliced-alignment', '-S',
                                        mapping_result_filename]
                os.system(' '.join(cmd_hisat2_alignment))
                
                shutil.rmtree(folder_name)
                #unmapped_f.write('>' + line[0] + '\n')
                #unmapped_f.write(line[9] + '\n')
                continue

            if len(line[0].split('_')) == 3:
                access = line[0].split('_')[0]
            else:
                access = '_'.join(line[0].split('_')[:-2])

            current_ref_range = copy.deepcopy(ref_ranges[symbol])
            align_start  = int(line[3]) - 1
            cigar_string = line[5]
            seq = line[9]
    
            cigar_segments = re.findall('\d+\w', cigar_string)


            assigned_seq = ''

            edited_seq = ''
            #previous_length = 0
            for cigar_segment in cigar_segments:

                current_cigar_seg_len   = int(cigar_segment[:-1])
                current_cigar_seg_label = cigar_segment[-1]

                if current_cigar_seg_label == 'M':
                    edited_seq += seq[:current_cigar_seg_len]
                    seq = seq[current_cigar_seg_len:]
                elif current_cigar_seg_label == 'D':
                    edited_seq += '*' * current_cigar_seg_len
                elif current_cigar_seg_label == 'S':
                    if (current_cigar_seg_len > 1) and (len(edited_seq) == 0):
                        seq = seq[current_cigar_seg_len:]
                        align_start -=current_cigar_seg_len
                    else:
                        edited_seq += seq[:current_cigar_seg_len]
                        seq = seq[current_cigar_seg_len:]
                elif current_cigar_seg_label == 'I':
                    I_position = align_start + len(edited_seq)
                    for seg_i, info in current_ref_range.items():
                        seg_start, seg_end = info[0]
                        if seg_start > I_position:
                            seg_start += current_cigar_seg_len
                        if seg_end > I_position:
                            seg_end += current_cigar_seg_len
                        current_ref_range[seg_i] = [(seg_start, seg_end), info[1]]

            last_end = 0
            align_end = align_start + len(edited_seq)
            align_range = set(range(align_start, align_end))
            indentation = align_start

            for seg_i, info in current_ref_range.items():
                seg_start, seg_end = info[0]
                seg_range = range(seg_start, seg_end + 1)
                if seg_end < align_start or seg_start > align_end:
                    continue
                overlap = align_range.intersection(seg_range)
                if len(overlap) >= 20 or len(overlap) == len(seg_range):

                    if last_end == 0:
                        start = 0
                    else:
                        start = last_end
                    end = seg_end - indentation + 1
                    seg_seq = edited_seq[int(start): int(end)]
                    #seg_seq = seg_seq.replace('~', '')
                    if seg_seq in seg_seq_access_dict[symbol][seg_i].keys():   
                        seg_seq_access_dict[symbol][seg_i][seg_seq].append(access)
                    else:
                        seg_seq_access_dict[symbol][seg_i][seg_seq] = [access]
                    present_seg.append(access)
                    if seg_i in exons:
                            present_exon.append(access)
                            
                    last_end = end
                    
                elif seg_start <= align_start and align_end >= seg_end:
                    seg_seq = edited_seq[int(align_start): int(align_end)]
                    if len(seg_seq) >= 10:
                        #seg_seq = seg_seq.replace('~', '')
                        if seg_seq in seg_seq_access_dict[symbol][seg_i].keys():   
                            seg_seq_access_dict[symbol][seg_i][seg_seq].append(access)
                        else:
                            seg_seq_access_dict[symbol][seg_i][seg_seq] = [access]
                        present_seg.append(access)
                        if seg_i in exons:
                            present_exon.append(access)
                        
        symbol_seg_access_dict[symbol] = present_seg
        symbol_exon_access_dict[symbol] = present_exon

                
    return(seg_seq_access_dict, symbol_seg_access_dict, symbol_exon_access_dict, back_mapped_path)

def find_exon_seg(ref_ranges):
    symbol_exon_dict = {}
    
    for symbol, segment in ref_ranges.items():
        symbol_exon_dict[symbol] = []
        for seg_i, info in segment.items():
            range_, label = info
            if label == 'exon':
                symbol_exon_dict[symbol].append(seg_i)
                
    return(symbol_exon_dict)

def make_msf(threads, symbol_seg_filenames, symbol_exon_filenames):
    pool = Pool(int(threads))
    
    for symbol, seg_names in symbol_seg_filenames.items():
        pool.map(run_clustalo, seg_names)

    for symbol, seg_names in symbol_exon_filenames.items():
        pool.map(run_clustalo_pile, seg_names)

    pool.close()

def write_fasta_by_seg(symbol_seg_dict, remapped_access_seq_dict, symbol_exon_dict, ref_seq, ref_ranges, out_dir):
    symbol_seg_filenames  = {}
    symbol_exon_filenames = {}
    ref_access_dict = get_RefSeqID(symbol_seg_dict.keys())
    
    datas = [symbol_seg_dict, remapped_access_seq_dict]
        
    for data in datas:
        for symbol, symbol_seg_dict_ in data.items():
            if data == symbol_seg_dict:
                symbol_seg_filenames[symbol] = []
                symbol_exon_filenames[symbol] = []
                
            file_path = out_dir + '/%s/segments/' % symbol
            if not os.path.exists(file_path):
                os.makedirs(file_path)
            file_path = out_dir + '/%s/exons/' % symbol
            if not os.path.exists(file_path):
                os.makedirs(file_path)

            ref_access = ref_access_dict[symbol][0]

            for seg_i, seg_seqs in symbol_seg_dict_.items():

                file_name = out_dir + '/%s/segments/%s_%s.fasta' % (symbol, symbol, seg_i)
                f = open(file_name, 'a')
                if data == symbol_seg_dict:
                    symbol_seg_filenames[symbol].append(file_name)

                    ref_start, ref_end = ref_ranges[symbol][seg_i][0]
                    ref_sequence = ref_seq[symbol][ref_start: ref_end + 1]

                    f.write('>' + ref_access + '\n')
                    f.write(ref_sequence + '\n')

                if seg_i in symbol_exon_dict[symbol]:
                    file_name_exon = out_dir + '/%s/exons/%s_%s.fasta' % (symbol, symbol, seg_i)
                    f_exon = open(file_name_exon, 'a')
                    if data == symbol_seg_dict:
                        symbol_exon_filenames[symbol].append(file_name_exon)
                        f_exon.write('>' + ref_access + '\n')
                        f_exon.write(ref_sequence + '\n')

                for seg_seq, access_list in seg_seqs.items():
                    for access in access_list:

                        f.write('>' + access + '\n')
                        f.write(seg_seq + '\n')
                        if seg_i in symbol_exon_dict[symbol]:  
                            f_exon.write('>' + access + '\n')
                            f_exon.write(seg_seq + '\n')

    return(symbol_seg_filenames, symbol_exon_filenames)

def combine_msf(ref_ranges, symbol_exon_dict, symbol_seg_access_dict, symbol_exon_access_dict, out_dir):
    
    for category in ['segments', 'exons']:
        if category == 'exons':
            symbol_target_access_dict = symbol_exon_access_dict
        else:
            symbol_target_access_dict = symbol_seg_access_dict
    
        symbol_access_whole_seq = {}
        ref_access_list = get_RefSeqID(ref_ranges.keys())

        for symbol in symbol_target_access_dict.keys():
            symbol_target_access_dict[symbol] = symbol_target_access_dict[symbol] + ref_access_list[symbol]

        for symbol, segment_info in ref_ranges.items():
            symbol_access_whole_seq[symbol] = {}
            ref_access = ref_access_list[symbol][0]

            for seg_i, info in segment_info.items():
                if category == 'exons' and seg_i not in symbol_exon_dict[symbol]:
                    continue
                
                filename = out_dir + '/%s/%s/%s_%s.msf' % (symbol, category, symbol, seg_i)
                check_exist = os.path.isfile(filename)
                access_seq  = {}
                if check_exist:
                    with open(filename) as f:  
                        for line in f:
                            line = line.strip()
                            line = line.split(' ')

                            line[0] = line[0].split('_')
                            if len(line[0]) < 1:
                                continue
                            line[0] = '_'.join(line[0][:-1]) + '.' + line[0][-1]
                            if line[0] in symbol_target_access_dict[symbol]:
                                tmp_sequence = ''
                                access = line[0]
                                for sequence_piece in line[1:]:
                                    if len(sequence_piece) > 0:
                                        tmp_sequence = tmp_sequence + sequence_piece

                                if access in access_seq.keys():
                                    access_seq[access] = access_seq[access] + tmp_sequence
                                else:
                                    access_seq[access] = tmp_sequence
                else:
                    filename = filename.replace('.msf', '.fasta')
                    with open(filename) as f:
                        for line in f:
                            if line.startswith('>'):
                                access = line[1:].strip()
                            else:
                                sequence = line.strip()
                    access_seq[access] = sequence

                seq_len = len(access_seq[ref_access])    
                for access in symbol_target_access_dict[symbol]:
                    if access in access_seq.keys():
                        continue
                    access_seq[access] = '~' * seq_len

                for access in symbol_target_access_dict[symbol]:
                    if access in symbol_access_whole_seq[symbol].keys():
                        symbol_access_whole_seq[symbol][access] = symbol_access_whole_seq[symbol][access] + access_seq[access]
                    else:
                        symbol_access_whole_seq[symbol][access] = access_seq[access]

        for symbol, access_whole_seq in symbol_access_whole_seq.items():
            ref_access = ref_access_list[symbol][0]
            seq_len  = len(symbol_access_whole_seq[symbol][ref_access])
            position = 0
            filename = out_dir + '/%s/%s_all_%s.msf' % (symbol, symbol, category)
            ofile    = open(filename, 'w')
            while position <= seq_len:
                for access, whole_sequence in symbol_access_whole_seq[symbol].items():
                    seq_list = []
                    tmp_seq  = whole_sequence[position: position + 50].replace('.', '~')
                    sub_position = 0
                    while sub_position != 50:
                        seq_list.append(tmp_seq[sub_position: sub_position + 10])
                        sub_position = sub_position + 10
                    ofile.write(access + '\t' + ' '.join(seq_list) + '\n')
                ofile.write('\n')
                position = position + 50

def clean_up(out_dir, gene_list):
    
    for symbol in gene_list:
        directory = out_dir + '/%s/' % symbol
        name_list = os.listdir(directory)
        for name in name_list:
            if name.endswith('_all_exons.msf') or name.endswith('_all_segments.msf'):
                continue
            try: 
                os.remove(directory + name)
            except:
                shutil.rmtree(directory + name)

def back_mapping(back_mapped_path, out_dir):
    access_seg_dict = {}
    access_seq_dict = {}
    #back_umapped_filenames = {}
    
    for symbol, path in back_mapped_path.items():
        path = back_mapped_path[symbol]

        filenames = os.listdir(path)


        access_seg_dict[symbol] = {}
        access_seq_dict[symbol] = {}

        back_umapped_seq_filename = out_dir + symbol + '/' + symbol + '_back_umapped_seq.fasta'
        back_umapped_f = open(back_umapped_seq_filename, 'w')
        start_positions = []
        for filename in filenames:

            full_filename = path + filename
            access, seg_type, seg = filename.split('_')

            f = open(full_filename)
            mapped_ref_seg_list = []
            mapped_ref_seg_seq  = {}
            
            
            for line in f:
                if line.startswith('@'):
                    continue

                line = line.split('\t')

                if line[2] == '*':
                    continue

                if seg_type == 'exon':
                    access_seq = access_seqs[symbol][access]['exons'][int(seg)]
                elif seg_type == 'intron':
                    access_seq = access_seqs[symbol][access]['introns'][int(seg)]
                else:
                    access_seq = access_seqs[symbol][access][seg_type + '_' + seg]

                ref_seg = int(line[0].split('_')[-1])

                align_start  = int(line[3]) - 1
                current_position = align_start
                cigar_string = line[5]
                ref_seq_len = len(line[9])

                cigar_segments = re.findall('\d+\w', cigar_string)

                assigned_seq = ''

                edited_seq = ''
                
                
                for cigar_segment in cigar_segments:
    
                    current_cigar_seg_len   = int(cigar_segment[:-1])
                    current_cigar_seg_label = cigar_segment[-1]

                    if current_cigar_seg_label == 'M':
                        current_position += current_cigar_seg_len

                    elif current_cigar_seg_label == 'S':
                        if current_cigar_seg_len > 1:
                            if current_position == align_start:
                                ref_seq_len -= current_cigar_seg_len
                                align_start += current_cigar_seg_len
                            else:
                                ref_seq_len -= current_cigar_seg_len
                        else:
                            current_position += current_cigar_seg_len

                    elif current_cigar_seg_label == 'I':
                        access_seq = access_seq[:current_position] + '*' * current_cigar_seg_len + access_seq[current_position:]
                        current_position += current_cigar_seg_len

                    seq = access_seq[align_start: align_start + ref_seq_len]
                    start_positions.append(align_start)

                    if (len(mapped_ref_seg_list) == 0) and align_start > 11:
                        back_umapped_f.write('>' + filename + '\n')
                        back_umapped_f.write(access_seq[:align_start - 1] + '\n')
                
                if 'N' in seq:
                    continue
                        
                mapped_ref_seg_list.append(ref_seg)
                mapped_ref_seg_seq[seq] = {ref_seg : access}
                if ref_seg in access_seq_dict[symbol]:
                    if seq in access_seq_dict[symbol][ref_seg]:
                        access_seq_dict[symbol][ref_seg][seq].append(access)
                    else:
                        access_seq_dict[symbol][ref_seg][seq] = [access]
                else:
                    access_seq_dict[symbol][ref_seg] = {seq: [access]}
                    #access_seq_dict[symbol][ref_seg][seq] = access
            
            if sorted(start_positions) != start_positions:
                continue
            
            if (len(mapped_ref_seg_list) != 0) and (len(access_seq) - (align_start + ref_seq_len) >= 10):
                back_umapped_f.write('>' + filename + '\n')
                back_umapped_f.write(access_seq[align_start + ref_seq_len:] + '\n')

            if len(mapped_ref_seg_list) != 0:
                if access in access_seg_dict[symbol]:
                    access_seg_dict[symbol][access] = access_seg_dict[symbol][access] + mapped_ref_seg_list
                else:
                    access_seg_dict[symbol][access] = mapped_ref_seg_list
                    
                
    return(access_seg_dict, access_seq_dict)


def re_mapping(symbol_list, out_dir, back_access_seg_dict, back_access_seq_dict):
    remapping_filenames = {}
    for symbol in symbol_list:
        unmapped_filename = out_dir + symbol + '/' + symbol + '_back_umapped_seq.fasta'

        if os.path.isfile(unmapped_filename):

            indexed_filename = indexed_filenames[symbol]
            alignment_filename = unmapped_filename.replace('_back_umapped_seq.fasta', '_remapping_result')
            cmd_hisat2_alignment = ['hisat2',
                                    '-x',indexed_filename,
                                    '-f',unmapped_filename,
                                    '--no-spliced-alignment', '-S',
                                    alignment_filename]
            
            os.system(' '.join(cmd_hisat2_alignment))
            #remapping_filenames[symbol] = alignment_filename
        
            remapping_f = open(alignment_filename)

            for line in remapping_f:

                if line.startswith('@'):
                    continue

                line = line.split('\t')
                if line[2] == '*':
                    continue
                #print(line)
                if len(line[0].split('_')) == 3:
                    access = line[0].split('_')[0]
                else:
                    access = '_'.join(line[0].split('_')[:-2])
                
                if access not in back_access_seq_dict[symbol]:
                    continue
                current_ref_range = copy.deepcopy(ref_ranges[symbol])
                align_start  = int(line[3]) - 1
                cigar_string = line[5]
                seq = line[9]

                cigar_segments = re.findall('\d+\w', cigar_string)

                assigned_seq = ''

                edited_seq = ''

                for cigar_segment in cigar_segments:

                    current_cigar_seg_len   = int(cigar_segment[:-1])
                    current_cigar_seg_label = cigar_segment[-1]

                    if current_cigar_seg_label == 'M':
                        edited_seq += seq[:current_cigar_seg_len]
                        seq = seq[current_cigar_seg_len:]
                    elif current_cigar_seg_label == 'D':
                        edited_seq += '*' * current_cigar_seg_len
                    elif current_cigar_seg_label == 'S':
                        if (current_cigar_seg_len > 1) and (len(edited_seq) == 0):
                            seq = seq[current_cigar_seg_len:]
                            align_start -=current_cigar_seg_len
                        else:
                            edited_seq += seq[:current_cigar_seg_len]
                            seq = seq[current_cigar_seg_len:]
                    elif current_cigar_seg_label == 'I':
                        I_position = align_start + len(edited_seq)
                        for seg_i, info in current_ref_range.items():
                            seg_start, seg_end = info[0]
                            if seg_start > I_position:
                                seg_start += current_cigar_seg_len
                            if seg_end > I_position:
                                seg_end += current_cigar_seg_len
                            current_ref_range[seg_i] = [(seg_start, seg_end), info[1]]

                last_end = 0
                align_end = align_start + len(edited_seq)
                align_range = set(range(align_start, align_end))
                indentation = align_start

                for seg_i, info in current_ref_range.items():
                    seg_start, seg_end = info[0]
                    seg_range = range(seg_start, seg_end + 1)
                    if seg_end < align_start or seg_start > align_end:
                        continue
                    overlap = align_range.intersection(seg_range)
                    if len(overlap) >= 10 or len(overlap) == len(seg_range):

                        if last_end == 0:
                            start = 0
                        else:
                            start = last_end
                        end = seg_end - indentation + 1
                        seg_seq = edited_seq[int(start): int(end)]

                        seg_list = back_access_seg_dict[symbol][access]
                        if (seg_i - max(seg_list) <= 2) or (mix(seg_list) - seg_i <= 2):
                            if seg_i in back_access_seq_dict[symbol]:
                                if seg_seq in back_access_seq_dict[symbol][seg_i]:
                                    back_access_seq_dict[symbol][seg_i][seg_seq].append(access)
                                else:
                                    back_access_seq_dict[symbol][seg_i][seg_seq] = [access]
                            else:
                                back_access_seq_dict[symbol][seg_i] = {seg_seq: [access]}


                        last_end = end

                    elif seg_start <= align_start and align_end >= seg_end:
                        seg_seq = edited_seq[int(align_start): int(align_end)]
                        if len(seg_seq) >= 10:
                            if (seg_i - max(seg_list) <= 2) or (mix(seg_list) - seg_i <= 2):
                                if seg_i in back_access_seq_dict[symbol]:
                                    if seg_seq in back_access_seq_dict[symbol][seg_i]:
                                        back_access_seq_dict[symbol][seg_i][seg_seq].append(access)
                                    else:
                                        back_access_seq_dict[symbol][seg_i][seg_seq] = [access]
                                else:
                                    back_access_seq_dict[symbol][seg_i] = {seg_seq: [access]}
    return(back_access_seq_dict)

def extract_gene_sequence(gene_list, out_dir, threads, verbose):
    
    symbol_uid_dic = build_symbol_uid_dic()
    for symbol in gene_list:
        symbol = [symbol]
   
        access_seqs, ref_ranges, ref_seq = get_sequence(symbol, symbol_uid_dic, verbose)

        ref_filenames = build_ref_fasta(ref_seq, out_dir)
        seq_filenames = build_seqs_fasta(access_seqs, out_dir)

        indexed_filenames = build_hisat2_index(ref_filenames)
        alignment_filenames = hisat2_alignment(seq_filenames, indexed_filenames)
        symbol_exon_dict = find_exon_seg(ref_ranges)
        symbol_seg_dict, symbol_seg_access_dict, symbol_exon_access_dict = assign_segment(ref_ranges, symbol_exon_dict, alignment_filenames)
        symbol_seg_filenames, symbol_exon_filenames = write_fasta_by_seg(symbol_seg_dict, symbol_exon_dict, ref_seq, ref_ranges, out_dir)

        make_msf(threads, symbol_seg_filenames, symbol_exon_filenames)

        combine_msf(ref_ranges, symbol_exon_dict, symbol_seg_access_dict, symbol_exon_access_dict, out_dir)
    clean_up(out_dir, gene_list)

if __name__ == '__main__':
    
    Entrez.email = 'Yun.Zhang@UTSouthWestern.edu'
    
    parser = ArgumentParser(
        description='Extract Gene and Related Senquence')
    parser.add_argument("--out-dir",
                        dest="out_dir",
                        type=str,
                        default=".",
                        help="Directory name for extracted read files")
    parser.add_argument("--gene-list",
                        dest="gene_list",
                        type=str,
                        default="",
                        help="A comma-separated list of genes (default: empty)")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="Verbose Output")

    args = parser.parse_args()

    gene_list = []
    if args.gene_list != '':
        gene_list = args.gene_list.upper().split(',')

    extract_gene_sequence(gene_list, 
                          args.out_dir, 
                          args.threads, 
                          args.verbose)
