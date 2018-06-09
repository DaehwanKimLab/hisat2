#!/usr/bin/env python

import sys
import struct

chr22_seq = ""
for line in open("22.fa"):
    if line.startswith('>'):
        continue
    line = line.strip()
    line = line.replace('N', '')
    chr22_seq += line
chr22_seq += '$'

chr22_sa = []
f = open("22.sa", "rb")
while True:
    fourbytes = f.read(4)
    if fourbytes == "":
        break

    num = struct.unpack('I', fourbytes)[0]
    chr22_sa.append(num)

    if len(chr22_sa) % 5000000 == 0:
        print len(chr22_sa)
f.close()

assert chr22_sa[0] + 1 == len(chr22_sa)
chr22_sa = chr22_sa[1:]
assert len(chr22_seq) == len(chr22_sa)

# for i, num in enumerate(chr22_sa):
#    print "%10d\t%s\t%d" % (i, chr22_seq[num:num+100], num)

seq_len = 100
i = 0
pos_set, pos_seq = None, ""
pos_set2, pos_seq2 = None, ""
while i < len(chr22_sa):
    pos = chr22_sa[i]
    base_seq = chr22_seq[pos:pos+seq_len]
    for j in xrange(i+1, len(chr22_sa)):
        pos2 = chr22_sa[j]
        cmp_seq = chr22_seq[pos2:pos2+seq_len]
        if base_seq != cmp_seq:
            break

    if j - i >= 200:
        if pos_set == None:
            pos_set = sorted(chr22_sa[i:j])
            pos_seq = base_seq
        else:
            pos_set2 = sorted(chr22_sa[i:j])
            pos_seq2 = base_seq
            different = True
            for _pos1 in pos_set:
                for _pos2 in pos_set2:
                    if abs(_pos1 - _pos2) < 1000:
                        different = False
                        break

            if different:
                """
                print pos_set
                print pos_set2
                print pos_seq
                print pos_seq2

                file1 = open("1.fa", "w")
                file2 = open("2.fa", "w")

                pos_seq2 = list(pos_seq2)
                pos_seq2 = pos_seq2[::-1]
                for k in xrange(seq_len):
                    nt = pos_seq2[k]
                    if nt == 'A':
                        nt = 'T'
                    elif nt == 'C':
                        nt = 'G'
                    elif nt == 'G':
                        nt = 'C'
                    else:
                        assert nt == 'T'
                        nt = 'A'
                    pos_seq2[k] = nt
                pos_seq2 = ''.join(pos_seq2)

                for k in xrange(1000000):
                    print >> file1, ">%d" % k
                    print >> file2, ">%d" % k
                    print >> file1, pos_seq
                    print >> file2, pos_seq2
                   
                file1.close()
                file2.close()
                """

                break

    i = j


chr22_seq = ""
for line in open("22.fa"):
    if line.startswith('>'):
        continue
    line = line.strip()
    chr22_seq += line

N_ranges = []
prev_nt = None
for i in xrange(len(chr22_seq)):
    nt = chr22_seq[i]
    if nt == 'N':
        if prev_nt != 'N':
            N_ranges.append([i, i]) # inclusive
        else:
            assert len(N_ranges) > 0
            N_ranges[-1][1] = i
    prev_nt = nt

to_joined_list = []
for N_start, N_end in N_ranges:
    if len(to_joined_list) == 0:
        if N_start > 0:
            to_joined_list.append([0, 0])
        else:
            to_joined_list.append([N_end + 1, 0])
    else:
        N_size = N_end - N_start + 1
        to = N_end + 1 - to_joined_list[-1][0] + to_joined_list[-1][1] - N_size
        assert to > to_joined_list[-1][1]
        to_joined_list.append([N_end + 1, to])

to_genome_list = [[y, x] for x, y in to_joined_list]

N_ranges_tmp = []
for i in xrange(len(to_genome_list)):
    to_genome = to_genome_list[i]
    if i == 0:
        if to_genome[1] > 0:
            N_ranges_tmp.append([0, to_genome[1] - 1])
    else:
        to_genome_before = to_genome_list[i-1]
        N_ranges_tmp.append([to_genome_before[1] + to_genome[0] - to_genome_before[0], to_genome[1] - 1])

assert N_ranges == N_ranges_tmp

file = open("22_rep.info", "w")
def print_rep_info(rep_name, rep_len, pos_set, pos_seq):
    print >> file, ">%s\t%d\t%d" % (rep_name, rep_len, len(pos_set))
    for i in xrange(0, len(pos_set), 10):
        output = ""
        for j in range(i, i + 10):
            if j >= len(pos_set):
                break
            if j > i:
                output += " "

            def convert(pos):
                for i in xrange(len(to_genome_list)):
                    if i + 1 == len(to_genome_list) or (pos >= to_genome_list[i][0] and pos < to_genome_list[i+1][0]):
                        return pos - to_genome_list[i][0] + to_genome_list[i][1]

                assert False
                
            pos = convert(pos_set[j])
            assert chr22_seq[pos:pos+seq_len] == pos_seq
            output += ("22:%d" % pos)
        print >> file, output
print_rep_info("rep1", seq_len, pos_set, pos_seq)
print_rep_info("rep2", seq_len, pos_set2, pos_seq2)
file.close()

chr22_seq = chr22_seq.replace(pos_seq, 'N' * seq_len)
chr22_seq = chr22_seq.replace(pos_seq2, 'N' * seq_len)
file = open("22_rep.fa", "w")
print >> file, ">22_rep"
for i in xrange(0, len(chr22_seq), 60):
    print >> file, chr22_seq[i:i+60]
print >> file, ">rep1"
print >> file, pos_seq
print >> file, ">rep2"
print >> file, pos_seq2
file.close()


