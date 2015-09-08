#!/usr/bin/env python

import sys, os
import string, re
use_message = '''
'''

def init():
    read_dir_base = "../reads/real/"
    read_dirs = os.listdir(read_dir_base)
    for read_dir in read_dirs:
        if os.path.exists(read_dir):
            continue
        if not os.path.exists(read_dir_base + read_dir + "/1.fq") or \
                not os.path.exists(read_dir_base + read_dir + "/1.fq"):
            continue

        print >> sys.stderr, "Processing", read_dir, "..."

        os.mkdir(read_dir)
        os.chdir(read_dir)

        RNA = (read_dir.find("RNA") != -1)
        tests = [
            ["1M", 1000000],
            ["5M", 5000000],
            ["20M", 20000000],
            ["whole", 0]
            ]

        for dir_name, num_reads in tests:
            if os.path.exists(dir_name):
                continue

            os.mkdir(dir_name)
            os.chdir(dir_name)

            if dir_name == "whole":
                ln_cmd = "ln -s ../../%s%s/[12].fq ." % (read_dir_base, read_dir)
                print >> sys.stderr, ln_cmd
                os.system(ln_cmd)
            else:
                head_cmd = "head -%d ../../%s%s/1.fq > 1.fq" % (num_reads * 4, read_dir_base, read_dir)
                print >> sys.stderr, head_cmd
                os.system(head_cmd)
                head_cmd = "head -%d ../../%s%s/2.fq > 2.fq" % (num_reads * 4, read_dir_base, read_dir)
                print >> sys.stderr, head_cmd
                os.system(head_cmd)

            os.system("ln -s ../../calculate_read_cost.py .")
            os.chdir("..")

        os.chdir("..")
    

if __name__ == "__main__":
    init()
