#!/usr/bin/env python3

import os, sys, math, random, re
from collections import defaultdict, Counter
from argparse import ArgumentParser, FileType
import pprint
import bisect

cigar_re = re.compile('\d+\w')

def read_len_cigar(cigars):
    read_len = 0
    for cigar in cigars:
        length, cigar_op = cigar
        if cigar_op in "MISH":
            read_len += int(length)

    return read_len


class Transinfo:
    def __init__(self):
        # tbl['tid,offset'] = (gene_id, tlen, ..., list of exons)
        self.tbl = {}
        # offset_lookup_table[tid] = [list of offsets, sorted]
        self.offset_lookup_table = {}

        # transcript-exon mapping table (per gene)
        self.tmap = {}

    def loadfromfile(self, info_fp):
        self.tbl = {}
        self.offset_lookup_table = {}

        current_trans = None
        # parse
        for line in info_fp:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('>'):
                line = line[1:].split('\t')

                # tid, chr_name, gene_id, offset, tran_gene_len
                current_tid = line[0]
                chr_name = line[1]

                strand = ''
                if len(line) < 3:
                    tran_gene_len = 100000000
                    gene_id = ''
                    offset = 0
                else:
                    gene_id = line[2]
                    offset = int(line[3])
                    tran_gene_len = int(line[4])

                tbl_key = '{},{}'.format(current_tid, offset)


                # add to offset_lookup_table
                if not current_tid in self.offset_lookup_table:
                    self.offset_lookup_table[current_tid] = [offset]
                else:
                    self.offset_lookup_table[current_tid].append(offset)

                if not tbl_key in self.tbl:
                    #                   chr_name, strand, len, exons, gene_id
                    self.tbl[tbl_key] = [chr_name, strand, tran_gene_len, list(), gene_id]
                    current_trans = self.tbl[tbl_key]

                else:
                    print('Duplicated tid', tbl_key)
                    current_trans = self.tbl[tbl_key]

            else:
                field = line.split('\t')

                for item in field:
                    chr_name, genomic_position, exon_len = item.split(':')[0:3]
                    current_trans[3].append([int(genomic_position), int(exon_len)])

        # sort offset_lookup_table
        for tid in self.offset_lookup_table:
            self.offset_lookup_table[tid].sort()


    def load_tmapfile(self, tmap_fp):
        self.tmap = {}

        current_tmap = None
        # parse file
        for line in tmap_fp:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('>'):
                # header
                line = line[1:].split('\t')

                # stome, chrname, stxname, bit length, number of tx
                stxname = line[2]
                bitlength = int(line[3])
                ntx = int(line[4])

                if stxname not in self.tmap:
                    self.tmap[stxname] = [stxname, bitlength, ntx, list()]
                    current_tmap = self.tmap[stxname]
                else:
                    print('Duplicated Supertranscript name', stxname)
                    current_tmap = self.tmap[stxname]

            else:
                if not current_tmap:
                    print('Wrong line', line)
                    continue

                # bitmask, tid
                fields = line.split('\t')

                bitmask = int(fields[0], 16)

                current_tmap[3].append([bitmask, fields[1]])


        # print(self.tmap)
        return


    def find_offset_in_lookup_table(self, tid, pos):
        # find a rightmost offset less than or equal to the pos
        offset_list = self.offset_lookup_table[tid]

        i = bisect.bisect_right(offset_list, pos)
        if i:
            return offset_list[i - 1]

        raise ValueError

    def map_position_internal(self, tid, tr_pos):

        offset = self.find_offset_in_lookup_table(tid, tr_pos)
        tbl_key = '{},{}'.format(tid, offset)

        trans = self.tbl[tbl_key]

        tr_pos -= offset
        new_chr = trans[0]
        exons = trans[3]

        assert len(exons) >= 1

        e_idx = -1
        new_pos = 0

        for i, exon in enumerate(exons):
            # find first exon
            if tr_pos >= exon[1]:
                tr_pos -= exon[1]
            else:
                new_pos = exon[0] + tr_pos
                e_idx = i
                break

        assert e_idx >= 0

        return new_chr, new_pos, e_idx, trans, offset

    def map_position(self, tid, tr_pos):
        new_chr, new_pos, e_idx, _, _ = self.map_position_internal(tid, tr_pos)
        return new_chr, new_pos

    def translate_pos_cigar(self, tid, tr_pos, cigar_str):
        new_chr, new_pos, e_idx, trans, offset = self.map_position_internal(tid, tr_pos)

        tid_list = list()

        #trans = self.tbl[tid]
        exons = trans[3]

        if cigar_str == '*':
            return new_chr, new_pos, cigar_str, tid_list

        cigars = cigar_re.findall(cigar_str)
        cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]

        read_len = read_len_cigar(cigars)
        tr_pos -= offset
        if tr_pos + read_len > trans[2]:
            # wrong transcript
            return new_chr, new_pos, cigar_str, tid_list

        assert tr_pos + read_len <= trans[2]

        tmp_cigar = list()

        # first exon
        exon_bit_list = list()
        exon_bit_list.append(e_idx)

        r_len = exons[e_idx][1] - (new_pos - exons[e_idx][0])
        for cigar in cigars:
            c_len, c_op = cigar
            c_len = int(c_len)

            if c_op in ['S']:
                tmp_cigar.append([c_len, c_op])
                continue

            while c_len > 0 and e_idx < len(exons):
                if c_len <= r_len:
                    r_len -= c_len
                    tmp_cigar.append([c_len, c_op])
                    break
                else:
                    c_len -= r_len
                    tmp_cigar.append([r_len, c_op])
                    """
                    if e_idx == len(exons):
                        print(tid)
                    """
                    if e_idx == len(exons) - 1:
                        # FIXME: parkch
                        # wrong mapping (across a transcript)
                        gap = c_len
                        e_idx += 1
                    else:
                        gap = exons[e_idx + 1][0] - (exons[e_idx][0] + exons[e_idx][1])
                        e_idx += 1
                        r_len = exons[e_idx][1]

                    if c_op in ['M']:
                        exon_bit_list.append(e_idx)

                tmp_cigar.append([gap, 'N'])


        ##
        if exon_bit_list[-1] == len(exons):
            del exon_bit_list[-1]

        gene_id = trans[4]
        tidmap = self.tmap[gene_id]

        def generate_bitmask(msb, lsb):
            bmask = 2 ** (msb + 1) - 1
            cmask = 2 ** (lsb) - 1
            return bmask & ~cmask


        def generate_bitvalue(bitlist):
            bv = 0
            for i in bitlist:
                bv += 2 ** i
            return bv

        bitmask = generate_bitmask(exon_bit_list[-1], exon_bit_list[0])
        bitvalue = generate_bitvalue(exon_bit_list)

        # find tid
        tid_list = list()
        for t in tidmap[3]:
            # print('{:x}\t{:x}\t{:x}'.format(t[0], bitmask, bitvalue))
            if t[0] & bitmask == bitvalue:
                tid_list.append(t[1])

        # print('{}/{}'.format(len(tid_list), len(tidmap[3])), tid_list)
        ##

        # clean new_cigars
        ni = 0

        for i in range(1, len(tmp_cigar)):
            if tmp_cigar[i][0] == 0:
                pass
            elif tmp_cigar[i][1] == 'S' and tmp_cigar[ni][1] == 'N':
                # replace ni to i
                tmp_cigar[ni] = tmp_cigar[i]
            elif tmp_cigar[i][1] == 'N' and tmp_cigar[ni][1] == 'S':
                # goto next
                pass
            elif tmp_cigar[ni][1] == tmp_cigar[i][1]:
                # merge to ni
                tmp_cigar[ni][0] += tmp_cigar[i][0]
            else:
                ni += 1
                tmp_cigar[ni] = tmp_cigar[i]

        new_cigar = list()
        for i in range(ni+1):
            new_cigar.append('{}{}'.format(tmp_cigar[i][0], tmp_cigar[i][1]))

        return new_chr, new_pos, ''.join(new_cigar), tid_list


class OutputQueue:

    def __init__(self, fp=sys.stdout):
        self.fp = fp

        # read id
        self.current_rid = ''
        # key = [new_chr, new_pos, new_cigar], value = [count, samline_fields]
        #self.alignments = {}
        self.alignments = list()
        # single/paired-end
        self.is_paired = False

        self.NH = 0
        self.left_NH = 0
        self.right_NH = 0

    def reset(self):
        self.current_rid = ''
        #self.alignments = {}
        self.alignments = list()
        self.is_paired = False
        self.NH = 0
        self.left_NH = 0
        self.right_NH = 0
        return

    """
    """
    def updateNH(self, fields, count):
        for i, tag in enumerate(fields):
            if isinstance(tag, str) and tag.startswith("NH:i:"):
                fields[i] = 'NH:i:{}'.format(str(count))
                break

        return

    def updateTITZ(self, fields, tid_list):
        if len(tid_list) == 0:
            return

        fields.append("TI:Z:{}".format(tid_list[0]))

        TO = '|'.join(tid_list[1:])

        if TO:
            fields.append("TO:Z:{}".format(TO))
        return


    def remove_dup(self):
        # not empty
        assert self.alignments

        #self.alignments.append([new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, flag, sam_fields])


        flag = self.alignments[0][5]

        if self.is_paired:
            # paired-end
            tmp_alignments = self.alignments

            #self.alignments = list()

            '''
            if len(tmp_alignments) % 2 > 0:
                print(tmp_alignments, file=sys.stderr)
            '''

            #assert len(tmp_alignments) % 2 == 0
            positions_set = set()

            left_alignments_list = list()
            right_alignments_list = list()

            def check_dup_append(alignments_list, key, value):
                '''
                for item in alignments_list:
                    if key == item[0]:
                        return
                '''
                alignments_list.append([key, value])
                return

            for i in range(0, len(tmp_alignments)):
                alignment = tmp_alignments[i]

                # chr_name, pos, cigar, pair_chr, pair_pos
                key = (alignment[0], alignment[1], alignment[2], alignment[3], alignment[4])
                value = alignment
                flag = alignment[5]

                # check_dup
                if flag & 0x40:
                    # first, left
                    check_dup_append(left_alignments_list, key, value)

                elif flag & 0x80:
                    # last, right
                    check_dup_append(right_alignments_list, key, value)

                else:
                    assert False

            #self.NH = len(self.alignments) // 2
            self.left_NH = len(left_alignments_list)
            self.right_NH = len(right_alignments_list)

            return
            # make paires
            i = 0
            j = 0
            while i < len(left_alignments_list) and j < len(right_alignments_list):
                break
                pass

        else:
            # single-end
            tmp_alignments = self.alignments
            self.alignments = list()
            positions_set = set()

            for alignment in tmp_alignments:
                key = (alignment[0], alignment[1], alignment[2]) # new_chr, new_pos, new_cigar
                if key not in positions_set:
                    positions_set.add(key)
                    self.alignments.append(alignment)

            self.NH = len(self.alignments)
        return

    def print_sam(self):
        if self.is_paired:
            left_count = self.left_NH
            right_count = self.right_NH

            for alignment in self.alignments:
                fields = alignment[6]
                if alignment[0] == alignment[3]:
                    new_pair_chr = '='
                else:
                    new_pair_chr = alignment[3]
                new_pair_pos = alignment[4]

                flag = int(fields[1])
                mapq = 60

                is_left = (flag & 0x40 == 0x40)
                count = left_count if is_left else right_count

                if count == 1:
                    mapq = 60
                else:
                    mapq = 1

                # update NH
                self.updateNH(fields, count)

                assert self.current_rid == '' or self.current_rid == fields[0]

                print('\t'.join([fields[0], str(flag), alignment[0], str(alignment[1] + 1), str(mapq), alignment[2], new_pair_chr, str(new_pair_pos + 1), *fields[8:]]), file=self.fp)

        else:
            count = self.NH
            #        self.alignments.append([new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, flag, sam_fields, tid_list])

            for alignment in self.alignments:
                fields = alignment[6]
                tid_list = alignment[7]
                if alignment[0] == alignment[3]:
                    new_pair_chr = '='
                else:
                    new_pair_chr = alignment[3]
                new_pair_pos = alignment[4]

                flag = int(fields[1])
                mapq = 60

                if count == 1:
                    mapq = 60
                else:
                    mapq = 1

                # update NH
                self.updateNH(fields, count)
                self.updateTITZ(fields, tid_list)

                assert self.current_rid == '' or self.current_rid == fields[0]

                print('\t'.join([fields[0], str(flag), alignment[0], str(alignment[1] + 1), str(mapq), alignment[2], new_pair_chr, str(new_pair_pos + 1), *fields[8:]]), file=self.fp)

    """
    """
    def add(self, flag, new_rid, new_chr, new_pos, new_cigar, sam_fields, new_pair_chr, new_pair_pos, tran_id_list):
        flag_paired = (flag & 0x1 == 0x1)
        flag_primary = (flag & 0x900 == 0)

        if len(self.alignments) == 0:
            # update flag
            self.is_paired = flag_paired
        else:
            assert self.is_paired == flag_paired

        #new_key = (new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, primary)

        assert self.current_rid == '' or self.current_rid == new_rid

        self.alignments.append([new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, flag, sam_fields, tran_id_list])

        """
        if new_key in self.alignments:
            self.alignments[new_key][0] += 1
        else:
            self.alignments[new_key] = [1, sam_fields, new_pair_chr, new_pair_pos]
        """


    """
    """
    def flush(self):
        if self.alignments:
            self.remove_dup()
            self.print_sam()

        self.reset()

    """
    """
    def push(self, flag, new_rid, new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, sam_fields, tran_id_list):
        if self.current_rid != new_rid:
            self.flush()
            self.current_rid = new_rid

        self.add(flag, new_rid, new_chr, new_pos, new_cigar, sam_fields, new_pair_chr, new_pair_pos, tran_id_list)
        return


def main(sam_file, transinfo_file, tmap_file):
    transinfo_table = Transinfo()
    transinfo_table.loadfromfile(transinfo_file)

    transinfo_table.load_tmapfile(tmap_file)


    outq = OutputQueue()

    for line in sam_file:
        line = line.strip()
        if not line:
            continue

        if line.startswith('@'):
            outq.flush()
            print(line)
            continue

        fields = line.split('\t')

        flag = int(fields[1])

        # read id
        rid = fields[0]

        # reference name = chromosome name...
        old_chr = fields[2]
        old_pos = int(fields[3])

        old_pair_chr = fields[6]
        old_pair_pos = int(fields[7])

        '''
        if rid == '989880':
            print(rid)
        '''

        # update pair-read
        if flag & 0x1:

            if old_chr == '*':
                new_chr = '*'
                new_pos = 0
                new_cigar = '*'

                new_pair_chr = '*'
                new_pair_pos = 0
            else:
                # translate locations
                old_pos -= 1
                old_pair_pos -= 1

                cigar_str = fields[5]

                if old_pair_chr == '*':
                    assert False

                if old_pair_chr == "=":
                    old_pair_chr = old_chr

                new_chr, new_pos, new_cigar, tid_list = transinfo_table.translate_pos_cigar(old_chr, old_pos, cigar_str)
                new_pair_chr, new_pair_pos = transinfo_table.map_position(old_pair_chr, old_pair_pos)

            outq.push(flag, rid, new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, fields, tid_list)

        else:
            # single
            if flag & 0x04:
                # unmapped
                outq.flush()
                print(line)
                continue

            old_pos -= 1
            cigar_str = fields[5]

            new_chr, new_pos, new_cigar, tid_list = transinfo_table.translate_pos_cigar(old_chr, old_pos, cigar_str)

            new_pair_chr = '*'
            new_pair_pos = 0

            outq.push(flag, rid, new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, fields, tid_list)

            #print(fields[0], new_chr, new_pos + 1, new_cigar)

    # end
    outq.flush()

    return


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Convert transcriptome-based positions to the corresponding genomic positions')

    parser.add_argument('sam_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input SAM file')

    parser.add_argument('transinfo_file',
                        nargs='?',
                        type=FileType('r'),
                        help='transcript position file')

    parser.add_argument('-t', '--tid-map',
                        dest='tmap_file',
                        type=FileType('r'),
                        help='Transcript Mapping file')

    args = parser.parse_args()

    if not args.sam_file or not args.transinfo_file:
        parser.print_help()
        exit(1)

    main(args.sam_file, args.transinfo_file, args.tmap_file)

