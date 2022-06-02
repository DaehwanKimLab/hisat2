#!/usr/bin/env python3

import os, sys, math, random, re
from collections import defaultdict, Counter
from argparse import ArgumentParser, FileType
import pprint
import bisect

cigar_re = re.compile('\d+\w')

bGeneSAM = False
bPrintGeneHeader = False

def read_len_cigar(cigars):
    read_len = 0
    for cigar in cigars:
        length, cigar_op = cigar
        if cigar_op in "MISH":
            read_len += int(length)

    return read_len


def bitmask_to_list(bitmask):
    bitlist = list()

    i = 0
    while bitmask > 0:
        if bitmask & 0x1:
            bitlist.append(i)

        bitmask = bitmask >> 1
        i += 1

    return bitlist

def generate_bitmask(msb, lsb):
    bmask = 2 ** (msb + 1) - 1
    cmask = 2 ** (lsb) - 1
    return bmask & ~cmask


def generate_bitvalue(bitlist):
    bv = 0
    for i in bitlist:
        bv += 2 ** i
    return bv


class Exon:
    """
    """

    def __init__(self, chr_name='', position=0, length=0):
        self.chr_name = chr_name
        self.position = position
        self.length = length
        return

    def __repr__(self):
        return "<Exon chr_name: {}, position: {}, length: {}>".format(self.chr_name, self.position, self.length)


class SuperTranscript:
    """
    """

    def __init__(self):
        self.stx_name = ''
        self.chr_name = ''
        self.gene_name = ''

        self.stome_name = ''
        self.offset_in_stome = 0

        self.stx_length = 0

        self.gene_length = 0

        self.exons = list()
        self.exon_lookup_table = list()

        self.transcripts = list()  # bitmask of transcript [exon bitmask, transcript_name, transcript_length]
        self.tx_bit_len = 0
        self.tx_bitmask = 0

        return

    def __repr__(self):
        return "<stome_name: {}, stx_name: {}, chr_name: {}, gene_name: {}, offset_in_stome: {}, stx_length: {}, gene_length: {}>" \
            .format(self.stome_name, self.stx_name, self.chr_name, self.gene_name, self.offset_in_stome,
                    self.stx_length, self.gene_length)

    def finish(self):
        # update gene_length
        left = self.exons[0].position
        right = self.exons[-1].position + self.exons[-1].length
        self.gene_length = right - left

        # build exon_lookup_table
        sum_exons = 0
        for e in self.exons:
            self.exon_lookup_table.append(sum_exons)
            sum_exons += e.length

        assert len(self.exons) == len(self.exon_lookup_table)
        assert sum_exons == self.stx_length, [self.stome_name, self.stx_name, sum_exons, self.stx_length]

        return

    def add_exon(self, exon_str):
        """
        chr_name:genomic_position:length
        """

        fields = exon_str.split(':')

        if len(fields) != 3:
            print('Wrong format:', exon_str, file=sys.stderr)
            return

        chr_name = fields[0]
        position = int(fields[1])
        length = int(fields[2])

        exon = Exon(chr_name, position, length)
        self.exons.append(exon)

    def set_tx_bit_len(self, tx_bit_len):
        self.tx_bit_len = tx_bit_len
        # self.tx_bitmask =

    def find_exon(self, stome_position, stome_cigar_str):
        """
        Find exons that match the position with the cigar string
        return:
          list of exon index, offset_in_first_exon, cigars
        """

        offset_in_stx = stome_position - self.offset_in_stome

        i = bisect.bisect_right(self.exon_lookup_table, offset_in_stx)
        if not i:
            print("No exon", file=sys.stderr)
            return None, None, None

        e_i = i - 1
        e = self.exons[e_i]

        exon_offset_in_stx = self.exon_lookup_table[e_i]
        offset_in_exon = offset_in_stx - exon_offset_in_stx

        assert offset_in_exon < e.length
        #
        # chr_position = e.position + offset_in_exon
        # gene_position = e.position - self.exons[0].position + offset_in_exon
        #
        # print(chr_position, gene_position)
        # print(e, offset, offset_in_exon)
        # print(e.position + offset_in_exon)

        #
        cigars = cigar_re.findall(stome_cigar_str)
        cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]
        read_len = read_len_cigar(cigars)

        # print(cigars, read_len)

        assert self.stx_length >= offset_in_stx + read_len

        new_cigars = list()
        exon_list = list()

        # r_len : remains in the exon
        r_len = e.length - offset_in_exon

        # append first exon
        exon_list.append(e_i)

        for cigar in cigars:
            c_len, c_op = cigar
            c_len = int(c_len)

            if c_op in ['S']:
                new_cigars.append([c_len, c_op])
                continue

            while c_len > 0 and e_i < len(self.exons):
                if c_len <= r_len:
                    r_len -= c_len
                    new_cigars.append([c_len, c_op])
                    break
                else:
                    c_len -= r_len
                    new_cigars.append([r_len, c_op])

                    # get next exon

                    if e_i == len(self.exons) - 1:
                        gap = c_len
                        e_i += 1
                    else:
                        gap = self.exons[e_i + 1].position - (self.exons[e_i].position + self.exons[e_i].length)
                        e_i += 1
                        r_len = self.exons[e_i].length

                    if c_op in ['M']:
                        exon_list.append(e_i)

                new_cigars.append([gap, 'N'])

        if exon_list[-1] == len(self.exons):
            del exon_list[-1]

        # clean up
        tmp_cigars = new_cigars
        new_cigars = list()
        ni = 0

        for i in range(1, len(tmp_cigars)):
            if tmp_cigars[i][0] == 0:
                pass
            elif tmp_cigars[i][1] == 'S' and tmp_cigars[ni][1] == 'N':
                # replace ni to i
                tmp_cigars[ni] = tmp_cigars[i]
            elif tmp_cigars[i][1] == 'N' and tmp_cigars[ni][1] == 'S':
                # goto next
                pass
            elif tmp_cigars[ni][1] == tmp_cigars[i][1]:
                # merge
                tmp_cigars[ni][0] += tmp_cigars[i][0]
            else:
                ni += 1
                tmp_cigars[ni] = tmp_cigars[i]

        for i in range(ni + 1):
            new_cigars.append('{}{}'.format(tmp_cigars[i][0], tmp_cigars[i][1]))

        # print(new_cigars, exon_list)

        return exon_list, offset_in_exon, new_cigars

    def find_transcripts(self, exon_index_list, offset_in_exon):
        """
        Find compatible transcripts
        return:
            list of [offset_in_transcript, transcript name]
        """
        transcripts = list()

        # print(self.tx_bit_len)
        # print(self.transcripts)

        bitmask = generate_bitmask(exon_index_list[-1], exon_index_list[0])
        bitvalue = generate_bitvalue(exon_index_list)


        def generate_bitstr(value, width):
            bitstr = '{:b}'.format(value)
            if len(bitstr) < width:

                bitstr = '0'*(width - len(bitstr)) + bitstr
            return bitstr

        for tx in self.transcripts:
            # print('   transcript', generate_bitstr(tx[0], self.tx_bit_len))
            # print('exon bitvalue', generate_bitstr(bitvalue, self.tx_bit_len))
            # print('       bimask', generate_bitstr(bitmask, self.tx_bit_len))

            if tx[0] & bitmask == bitvalue:
                bitlist = bitmask_to_list(tx[0])
                offset = 0
                for i in bitlist:
                    if i >= exon_index_list[0]:
                        break

                    offset += self.exons[i].length

                offset += offset_in_exon
                transcripts.append([offset, tx[1]])

        return transcripts


class MapResult:
    """
    mapping result
        - stome_name, stome_position, ....
        - chromosome name, position in chromosome, cigar, ....
        - gene name, position in gene, cigar, ....
        - list of (transcript name, position in transcript)
    """

    def __init__(self, stome_name='', stome_position=0, stome_cigar=''):
        self.stome_name = stome_name
        self.stome_position = stome_position
        self.stome_cigar = stome_cigar

        self.chr_name = ''
        self.position = 0
        self.cigar = ''

        self.gene_name = ''
        self.gene_position = 0
        self.gene_cigar = ''

        # list of [transcript_name, transcript_position] tuple
        self.transcripts = list()

        return

    def reset(self, stome_name='', stome_position=0, stome_cigar=''):
        self.stome_name = stome_name
        self.stome_position = stome_position
        self.stome_cigar = stome_cigar

        self.chr_name = ''
        self.position = 0
        self.cigar = ''

        self.gene_name = ''
        self.gene_position = 0
        self.gene_cigar = ''

        # list of [transcript_name, transcript_position] tuple
        self.transcripts = list()

        return

    def __repr__(self):
        return "<stome[name: {}, position: {}, cigar: {}], chr[name: {}, position: {}, cigar: {}], gene[name: {}, position: {}, cigar: {}], {}" \
            .format(self.stome_name, self.stome_position, self.stome_cigar, self.chr_name, self.position, self.cigar,
                    self.gene_name, self.gene_position, self.gene_cigar, self.transcripts)


class SuperTranscriptIndex:
    """
    convert supertranscriptome alignment result to Chromosome, Gene, Transcript result
    """

    def __init__(self):
        # dict of SuperTranscript. stx_name as key
        self.supertranscripts = dict()

        # map of offset on SuperTranscriptome(eg, 22_tome, ...)
        # used in looking up SuperTranscript
        self.offset_lookup_table = {}
        self.offset_list = {}

        return

    def loadfromfile(self, map_fp, tmap_fp):

        current_stx = None
        for line in map_fp:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('>'):
                # New supertranscript start
                fields = line.split()

                current_stx = SuperTranscript()

                current_stx.stome_name = fields[0][1:]
                current_stx.stx_name = fields[2]
                current_stx.chr_name = fields[1]
                current_stx.gene_name = fields[2]
                current_stx.offset_in_stome = int(fields[3])
                current_stx.stx_length = int(fields[4])

                stx_name = current_stx.stx_name

                if stx_name in self.supertranscripts:
                    print('Duplicated supertranscript', line, file=sys.stderr)
                    current_stx = None
                    continue

                self.supertranscripts[stx_name] = current_stx

            else:

                if not current_stx:
                    continue

                fields = line.split('\t')
                for e in fields:
                    current_stx.add_exon(e)

        for stx_name, stx in self.supertranscripts.items():
            if stx.stome_name in self.offset_lookup_table:
                self.offset_lookup_table[stx.stome_name].append([stx.offset_in_stome, stx])
            else:
                self.offset_lookup_table[stx.stome_name] = [[stx.offset_in_stome, stx]]

        # sort
        for stome_name in self.offset_lookup_table:
            self.offset_lookup_table[stome_name].sort()
            self.offset_list[stome_name] = [i[0] for i in self.offset_lookup_table[stome_name]]

        # load tmap file
        current_stx = None
        for line in tmap_fp:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('>'):
                fields = line.split()

                stome_name = fields[0][1:]
                chr_name = fields[1]
                stx_name = fields[2]
                gene_name = fields[2]

                tx_bit_len = int(fields[3])
                num_transcripts = int(fields[4])

                # lookup supertranscripts
                if stx_name in self.supertranscripts:
                    current_stx = self.supertranscripts[stx_name]
                else:
                    current_stx = None
                    print('No supertranscript:', stx_name, file=sys.stderr)
                    continue

                # current_stx.tx_bit_len = tx_bit_len
                current_stx.set_tx_bit_len(tx_bit_len)

            else:
                if not current_stx:
                    print('Wrong line:', line)
                    continue

                # bitmask, transcript_id
                fields = line.split('\t')
                bitmask = int(fields[0], 16)
                transcript_id = fields[1]
                bitlist = bitmask_to_list(bitmask)
                transcript_length = 0

                # msb index
                assert bitlist[-1] < len(current_stx.exons)
                for i in bitlist:
                    transcript_length += current_stx.exons[i].length

                # print(current_stx.exons)
                # print(fields[0], bitmask_to_list(bitmask))
                # print(transcript_id, transcript_length)

                current_stx.transcripts.append([bitmask, transcript_id, transcript_length])

        for stx_name, stx in self.supertranscripts.items():
            stx.finish()

        return

    def find_supertranscript(self, stome_name: str, offset: int) -> SuperTranscript:
        """
        find supertranscript

        """

        offset_list = self.offset_list[stome_name]

        i = bisect.bisect_right(offset_list, offset)
        if i:
            return self.offset_lookup_table[stome_name][i - 1][1]

        return None

    def print(self):
        print(self.supertranscripts)

    def mapto(self, stome_name, stome_position, stome_cigar, result: MapResult) -> bool:
        result.reset(stome_name, stome_position, stome_cigar)

        # find supertranscripts
        stx = self.find_supertranscript(stome_name, stome_position)
        if not stx:
            print("Can't fild supertranscript", file=sys.stderr)
            return False

        # print(stx)
        # position_in_stome = stome_position - stx.offset_in_stome

        exon_index_list, offset_in_exon, new_cigars = stx.find_exon(stome_position, stome_cigar)
        if not exon_index_list:
            print("Can't find exon", file=sys.stderr)
            return False

        exon_index = exon_index_list[0]
        e = stx.exons[exon_index]
        transcripts = stx.find_transcripts(exon_index_list, offset_in_exon)

        # print(e, offset_in_exon)
        new_cigar_str = ''.join(new_cigars)

        result.chr_name = e.chr_name
        result.position = e.position + offset_in_exon
        result.cigar = new_cigar_str

        result.gene_name = stx.gene_name
        result.gene_position = e.position - stx.exons[0].position + offset_in_exon
        result.gene_cigar = new_cigar_str

        for transcript in transcripts:
            result.transcripts.append(transcript)

        # print(stx)
        # print(stx.gene_name, stx.gene_length)
        # print(transcripts)

        # gene_position = result.position - stx.offset_in_stome

        # result.gene_name = stx.gene_name
        # result.gene_position =

        return True


class OutputQueue:

    def __init__(self, fp=sys.stdout):
        self.fp = fp

        # read id
        self.current_rid = ''
        # key = [new_chr, new_pos, new_cigar], value = [count, samline_fields]
        # self.alignments = {}
        self.alignments = list()
        # single/paired-end
        self.is_paired = False

        self.NH = 0
        self.left_NH = 0
        self.right_NH = 0

    def reset(self):
        self.current_rid = ''
        # self.alignments = {}
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

    def updateTITO(self, fields, tid_list):
        if len(tid_list) == 0:
            return

        fields.append("TI:Z:{}-{}".format(tid_list[0][0], tid_list[0][1] + 1))

        TO = '|'.join(['{}-{}'.format(t[0], t[1] + 1) for t in tid_list[1:]])

        if TO:
            fields.append("TO:Z:{}".format(TO))
        return

    def remove_dup(self):
        # not empty
        assert self.alignments

        # self.alignments.append([new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, flag, sam_fields])

        flag = self.alignments[0][5]

        if self.is_paired:
            # paired-end
            tmp_alignments = self.alignments

            # self.alignments = list()

            '''
            if len(tmp_alignments) % 2 > 0:
                print(tmp_alignments, file=sys.stderr)
            '''

            # assert len(tmp_alignments) % 2 == 0
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

            # self.NH = len(self.alignments) // 2
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
                key = (alignment[0], alignment[1], alignment[2])  # new_chr, new_pos, new_cigar
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

                print('\t'.join(
                    [fields[0], str(flag), alignment[0], str(alignment[1] + 1), str(mapq), alignment[2], new_pair_chr,
                     str(new_pair_pos + 1), *fields[8:]]), file=self.fp)

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
                self.updateTITO(fields, tid_list)

                assert self.current_rid == '' or self.current_rid == fields[0]

                print('\t'.join(
                    [fields[0], str(flag), alignment[0], str(alignment[1] + 1), str(mapq), alignment[2], new_pair_chr,
                     str(new_pair_pos + 1), *fields[8:]]), file=self.fp)

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

        # new_key = (new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, primary)

        assert self.current_rid == '' or self.current_rid == new_rid

        self.alignments.append(
            [new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, flag, sam_fields, tran_id_list])

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
    stx_index = SuperTranscriptIndex()
    stx_index.loadfromfile(transinfo_file, tmap_file)

    outq = OutputQueue()
    result = MapResult()

    # pass header parts
    for line in sam_file:
        line = line.strip()
        if not line:
            continue

        if not line.startswith('@'):
            break

        print(line)

    if bPrintGeneHeader:
        # Additional header info
        for stx_name, stx in stx_index.supertranscripts.items():
            print('@SQ\tSN:{}\tLN:{}\tAL:'.format(stx_name, stx.stx_length))

        for stx_name, stx in stx_index.supertranscripts.items():
            for tx in stx.transcripts:
                print('@SQ\tSN:{}\tLN:{}\tGE:{}'.format(tx[1], tx[2], stx_name))

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


                # new_chr, new_pos, new_cigar, tid_list = transinfo_table.translate_pos_cigar(old_chr, old_pos, cigar_str)
                # new_pair_chr, new_pair_pos = transinfo_table.map_position(old_pair_chr, old_pair_pos)

            # outq.push(flag, rid, new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, fields, tid_list)

        else:
            # single
            if flag & 0x04:
                # unmapped
                outq.flush()
                print(line)
                continue

            old_pos -= 1
            cigar_str = fields[5]

            stx_index.mapto(old_chr, old_pos, cigar_str, result)

            if bGeneSAM:
                new_chr = result.gene_name
                new_pos = result.gene_position
                new_cigar = result.gene_cigar

            else:

                new_chr = result.chr_name
                new_pos = result.position
                new_cigar = result.cigar

            tid_list = list()

            for tid in result.transcripts:
                tid_list.append([tid[1], tid[0]])

            # new_chr, new_pos, new_cigar, tid_list = transinfo_table.translate_pos_cigar(old_chr, old_pos, cigar_str)

            new_pair_chr = '*'
            new_pair_pos = 0

            outq.push(flag, rid, new_chr, new_pos, new_cigar, new_pair_chr, new_pair_pos, fields, tid_list)

            # print(fields[0], new_chr, new_pos + 1, new_cigar)

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

    parser.add_argument('--gene',
                        dest='bGeneSAM',
                        action='store_true',
                        default=False,
                        help='Generate Gene-based SAM file')

    parser.add_argument('--gene-header',
                        dest='bPrintGeneHeader',
                        action='store_true',
                        default=False,
                        help='Print Gene, Transcript information')

    args = parser.parse_args()

    if not args.sam_file or not args.transinfo_file:
        parser.print_help()
        exit(1)

    bGeneSAM = args.bGeneSAM
    bPrintGeneHeader = args.bPrintGeneHeader
    main(args.sam_file, args.transinfo_file, args.tmap_file)
