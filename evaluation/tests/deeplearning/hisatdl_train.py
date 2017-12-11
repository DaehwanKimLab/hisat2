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
import tensorflow as tf
import inspect
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common


"""
"""
def train_sequences(base_fname,
                    locus_list,
                    verbose):
    node1 = tf.constant(3.0, dtype=tf.float32)
    node2 = tf.constant(4.0) # also tf.float32 implicitly
    print(node1, node2)
    sess = tf.Session()
    print(sess.run([node1, node2]))
    
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Generate train sequences")
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="hla",
                        help="base filename for backbone sequence, variants, and linking info (Default: hla)")
    parser.add_argument("--locus-list",
                        dest="locus_list",
                        type=str,
                        default="",
                        help="A comma-separated list of gene names (default: empty, all genes)")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.locus_list == "":
        locus_list = []
    else:
        locus_list = args.locus_list.split(',')
             
    if args.base_fname.find('/') != -1:
        elems = args.base_fname.split('/')
        base_fname = elems[-1]
        base_dname = '/'.join(elems[:-1])
    else:
        base_fname = args.base_fname
        base_dname = ""
        
    train_sequences(base_fname,
                    locus_list,
                    args.verbose)

