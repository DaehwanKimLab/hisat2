#!/usr/bin/env python

import sys, os, signal
import string, re

signal.signal(signal.SIGPIPE, signal.SIG_DFL)
use_message = '''
'''

def make_cat_cmd(gzmode, read_dir_base, read_dir, fq_name, num_read):
    cmd = []
    if gzmode:
        cmd += ["zcat"]
    else:
        cmd += ["cat"]

    cmd += ["../../%s%s/%s" % (read_dir_base, read_dir, fq_name)]
    cmd += ["|", "head", "-n", "%d" % (num_read * 4)]

    if gzmode:
        cmd += ["|", "gzip"]

    cmd += [">", fq_name]
    return ' '.join(cmd)


def init():
    read_dir_base = "../reads/real/"
    read_dirs = os.listdir(read_dir_base)
    for read_dir in read_dirs:
        if os.path.exists(read_dir):
            continue

        gz_file = False
        fq_1_name = '1.fq'
        fq_2_name = '2.fq'
        if os.path.exists(read_dir_base + read_dir + "/1.fq.gz") and \
            os.path.exists(read_dir_base + read_dir + "/2.fq.gz"):
            gz_file = True
            fq_1_name = '1.fq.gz'
            fq_2_name = '2.fq.gz'
        else:
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
            ["whole", 0],
            ]

        for dir_name, num_reads in tests:
            if os.path.exists(dir_name):
                continue

            os.mkdir(dir_name)
            os.chdir(dir_name)

            if dir_name == "whole":
                ln_cmd = "ln -s ../../%s%s/%s ." % (read_dir_base, read_dir, fq_1_name)
                print >> sys.stderr, ln_cmd
                os.system(ln_cmd)
                ln_cmd = "ln -s ../../%s%s/%s ." % (read_dir_base, read_dir, fq_2_name)
                print >> sys.stderr, ln_cmd
                os.system(ln_cmd)
            else:
                cmd = make_cat_cmd(gz_file, read_dir_base, read_dir, fq_1_name, num_reads)
                print >> sys.stderr, cmd
                os.system(cmd)

                cmd = make_cat_cmd(gz_file, read_dir_base, read_dir, fq_2_name, num_reads)
                print >> sys.stderr, cmd
                os.system(cmd)

            os.system("ln -s ../../calculate_read_cost.py .")
            os.chdir("..")

        os.chdir("..")
    

if __name__ == "__main__":
    init()
