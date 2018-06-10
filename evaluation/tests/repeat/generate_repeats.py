#!/usr/bin/env python

import sys
import struct

chr_name = "20"

chr_seq = ""
for line in open("%s.fa" % chr_name):
    if line.startswith('>'):
        continue
    line = line.strip()
    line = line.replace('N', '')
    chr_seq += line
chr_seq += '$'

chr_sa = []
f = open("%s.sa" % chr_name, "rb")
while True:
    fourbytes = f.read(4)
    if fourbytes == "":
        break

    num = struct.unpack('I', fourbytes)[0]
    chr_sa.append(num)

    if len(chr_sa) % 5000000 == 0:
        print len(chr_sa)
f.close()

assert chr_sa[0] + 1 == len(chr_sa)
chr_sa = chr_sa[1:]
assert len(chr_seq) == len(chr_sa)

# for i, num in enumerate(chr_sa):
#    print "%10d\t%s\t%d" % (i, chr_seq[num:num+100], num)

seq_len = 100
i = 0
repeats = []
while i < len(chr_sa) - 1:
    pos = chr_sa[i]
    base_seq = chr_seq[pos:pos+seq_len]
    for j in xrange(i+1, len(chr_sa)):
        pos2 = chr_sa[j]
        cmp_seq = chr_seq[pos2:pos2+seq_len]
        if base_seq != cmp_seq:
            break

    if j - i >= 200:
        repeats.append([base_seq, sorted(chr_sa[i:j])])

    i = j

    if i % 5000000 == 0:
        print i

found = False
print len(repeats), "repeats"
for i in xrange(len(repeats) - 1):
    for j in xrange(i + 1, len(repeats)):
        num_close = 0
        pos_seq, pos_set = repeats[i]
        pos_seq2, pos_set2 = repeats[j]
        for _pos1 in pos_set:
            for _pos2 in pos_set2:
                if abs(_pos1 - _pos2) < 300:
                    num_close += 1
                
        if num_close == 1:
            found = True
            print pos_set
            print pos_set2
            print pos_seq
            print pos_seq2

            file1 = open("1.fa", "w")
            file2 = open("2.fa", "w")

            pos_seq2_rc = list(pos_seq2)
            pos_seq2_rc = pos_seq2_rc[::-1]
            for k in xrange(seq_len):
                nt = pos_seq2_rc[k]
                if nt == 'A':
                    nt = 'T'
                elif nt == 'C':
                    nt = 'G'
                elif nt == 'G':
                    nt = 'C'
                else:
                    assert nt == 'T'
                    nt = 'A'
                pos_seq2_rc[k] = nt
            pos_seq2_rc = ''.join(pos_seq2_rc)

            for k in xrange(1000000):
                print >> file1, ">%d" % k
                print >> file2, ">%d" % k
                print >> file1, pos_seq
                print >> file2, pos_seq2_rc

            file1.close()
            file2.close()

            break
    if found:
        break

chr_seq = ""
for line in open("%s.fa" % chr_name):
    if line.startswith('>'):
        continue
    line = line.strip()
    chr_seq += line

N_ranges = []
prev_nt = None
for i in xrange(len(chr_seq)):
    nt = chr_seq[i]
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

file = open("%s_rep.info" % chr_name, "w")
def print_rep_info(rep_name, rep_len, pos_set, pos_seq):
    print >> file, ">%s*0\t%d\t%d\t0" % (rep_name, rep_len, len(pos_set))
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
            assert chr_seq[pos:pos+seq_len] == pos_seq
            output += ("%s_rep:%d" % (chr_name, pos))
        print >> file, output
print_rep_info("rep1", seq_len, pos_set, pos_seq)
print_rep_info("rep2", seq_len, pos_set2, pos_seq2)
file.close()

chr_seq = chr_seq.replace(pos_seq, 'N' * seq_len)
chr_seq = chr_seq.replace(pos_seq2, 'N' * seq_len)
file = open("%s_rep.fa" % chr_name, "w")
print >> file, ">%s_rep" % chr_name
for i in xrange(0, len(chr_seq), 60):
    print >> file, chr_seq[i:i+60]
print >> file, ">rep1"
print >> file, pos_seq
print >> file, ">rep2"
print >> file, pos_seq2
file.close()


