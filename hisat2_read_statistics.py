#!/usr/bin/env python3

#
# Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
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

import os, sys, math, gzip, bz2
from argparse import ArgumentParser, FileType
"""
"""
COMPRESSION_NON   = 0
COMPRESSION_GZIP  = 1
COMPRESSION_BZIP2 = 2

SEQUENCE_UNKNOWN  = -1
SEQUENCE_FASTA    = 0
SEQUENCE_FASTQ    = 1

FASTA_EXTENSIONS = ["fa", "fasta", "fna"]
FASTQ_EXTENSIONS = ["fq", "fastq"]

MAX_SKIP_LINES = 10000
"""
"""
def parser_FQ(fp):
    # skip empty line
    skip_line_count = 0
    while skip_line_count < MAX_SKIP_LINES:
        line = fp.readline()

        if line == "":
            # end of file
            return

        if line[0] == '@':
            break

        skip_line_count += 1

    if skip_line_count == MAX_SKIP_LINES:
        raise ValueError("Invalid file format")

    while True:
        id = line[1:].split()[0]
        seq = ""

        line = fp.readline()
        if line == "":
            return

        seq = line.strip()
        yield id, seq

        line = fp.readline()  # '+'
        line = fp.readline()  # quality
        line = fp.readline()  # next ID
        if line == "":
            return

"""
"""
def parser_FA(fp):
    # skip empty line
    skip_line_count = 0
    while skip_line_count < MAX_SKIP_LINES:
        line = fp.readline()

        if line == "":
            # end of file
            return

        if line[0] == '>':
            break

        skip_line_count += 1

    if skip_line_count == MAX_SKIP_LINES:
        raise ValueError("Invalid file format")

    while True:
        id = line[1:].split()[0]
        seq = ""

        while True:
            line = fp.readline()
            if line == "":
                break

            if line[0] == '>':
                break

            seq += line.strip()

        yield id, seq

        if line == "":
            return


"""
"""
def parse_type(fname):
    compression_type = COMPRESSION_NON
    sequence_type = SEQUENCE_UNKNOWN

    ff = fname.split('.')

    ext = ff[-1]
    if ext.lower() == "gz":
        compression_type = COMPRESSION_GZIP
        ext = ff[-2]
    elif ext.lower() == "bz2":
        compression_type = COMPRESSION_BZIP2
        ext = ff[-2]

    if ext.lower() in FASTA_EXTENSIONS:
        sequence_type = SEQUENCE_FASTA
    elif ext.lower() in FASTQ_EXTENSIONS:
        sequence_type = SEQUENCE_FASTQ

    return sequence_type, compression_type

"""
"""
def generate_stats(length_map):
    mn = 0  # minimun read length
    mx = 0  # maximum read length
    cnt = 0  # number of reads
    avg = 0  # average read length

    sum = 0

    if len(length_map) == 0:
        return cnt, mn, mx, avg

    # sort keys
    sorted_map = sorted(length_map)

    mn = sorted_map[0]
    mx = sorted_map[-1]

    for k, v in length_map.items():
        sum += k * v
        cnt += v

    avg = sum // cnt

    return cnt, mn, mx, avg

"""
"""
def reads_stat(read_file, read_count):
    length_map = {}
    try:
        sequence_type, compression_type = parse_type(read_file)

        if compression_type == COMPRESSION_GZIP:
            fp = gzip.open(read_file, 'rt')
        elif compression_type == COMPRESSION_BZIP2:
            fp = bz2.open(read_file, 'rt')
        else:
            assert (compression_type == COMPRESSION_NON)
            fp = open(read_file, 'r')

        if sequence_type == SEQUENCE_FASTA:
            fstream = parser_FA(fp)
        elif sequence_type == SEQUENCE_FASTQ:
            fstream = parser_FQ(fp)
        else:
            raise ValueError("Unsupported file format")

        cnt = 0
        for id, seq in fstream:
            l = len(seq)
            if l in length_map:
                length_map[l] += 1
            else:
                length_map[l] = 1

            cnt += 1
            if read_count > 0 and cnt >= read_count:
                break

        fp.close()

    except BaseException as e:
        print("Warning: {}".format(e), file=sys.stderr)

    cnt, mn, mx, avg =  generate_stats(length_map)
    # sort by (read count, read length)
    length_map = sorted(length_map.items(), key=lambda t: (t[1], t[0]), reverse=True)
    if len(length_map) == 0:
        length_map.append((0, 0))
    print(cnt, mn, mx, avg, ",".join([str(k) for (k,v) in length_map]))


if __name__ == '__main__':

    parser = ArgumentParser(
            description='Compute statistics of reads. Show number of reads and minimum, maximum, average length of reads')

    parser.add_argument('read_file',
                        nargs='?',
                        type=str,
                        help='reads file')

    parser.add_argument('-n',
                        dest='read_count',
                        action='store',
                        type=int,
                        default=10000,
                        help='reads count (default: 10000)')

    args = parser.parse_args()

    if not args.read_file:
        parser.print_help()
        exit(1)

    reads_stat(args.read_file, args.read_count)

