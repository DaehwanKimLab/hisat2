#!/usr/bin/env python

import sys, os
use_message = '''
'''

def get_data():
    mac = (sys.platform == "darwin")
    if not os.path.exists("aligners"):
        os.mkdir("aligners")
    os.chdir("aligners")
    if not os.path.exists("bin"):
        os.mkdir("bin")
    aligners = ["HISAT", "Bowtie2", "Bowtie", "STAR", "GSNAP"]
    for aligner in aligners:
        if aligner == "HISAT":
            dir = "hisat-0.1.6-beta"
            if os.path.exists(dir):
                continue
            fname = dir + "-source.zip"
            url = "http://www.ccb.jhu.edu/software/hisat/downloads"
            bins =  "hisat-align-s hisat-build-s hisat-inspect-s"
            installs = bins + " hisat hisat-build hisat-inspect"
            cmd = "wget %s/%s; unzip %s; cd %s; make %s; cp %s ../bin; cd .." % \
                (url, fname, fname, dir, bins, installs)
        elif aligner == "Bowtie2":
            dir = "bowtie2-2.2.5"
            if os.path.exists(dir):
                continue
            fname = dir + "-source.zip"
            url = "http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5"
            bins = "bowtie2-align-s bowtie2-build-s bowtie2-inspect-s"
            installs = bins + " bowtie2 bowtie2-build bowtie2-inspect"
            cmd = "wget %s/%s; unzip %s; cd %s; make %s; cp %s ../bin; cd .." % \
                (url, fname, fname, dir, bins, installs)
        elif aligner == "Bowtie":
            dir = "bowtie-1.1.2"
            if os.path.exists(dir):
                continue
            fname = dir + "-src.zip"
            url = "http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2"
            bins = "bowtie-align-s bowtie-build-s bowtie-inspect-s"
            installs = bins + " bowtie bowtie-build bowtie-inspect"
            cmd = "wget %s/%s; unzip %s; cd %s; make %s; cp %s ../bin; cd .." % \
                (url, fname, fname, dir, bins, installs)
        elif aligner == "STAR":
            dir = "STAR_2.4.2a"
            if os.path.exists("STAR-" + dir):
                continue
            fname = dir + ".tar.gz"
            url = "https://github.com/alexdobin/STAR/archive"
            if mac:
                add_cmd = "awk '{if($1 ~ /^CXX/) {print \"CXX =/opt/local/bin/g++-mp-4.8\";} else {print;}}' Makefile > Makefile.tmp; mv Makefile.tmp Makefile"
                make_arg = "STARforMac"
                cmd = "wget %s/%s; tar xvzf %s; cd STAR-%s/source; %s; make; make %s; cp STAR ../../bin; cd ../.." % \
                    (url, fname, fname, dir, add_cmd, make_arg)
            else:
                cmd = "wget %s/%s; tar xvzf %s; cd STAR-%s/source; make; cp STAR ../../bin; cd ../.." % \
                    (url, fname, fname, dir)
        elif aligner == "GSNAP":
            dir = "gmap-2015-07-23"
            dir2 = "gmap-gsnap-2015-07-23"
            if os.path.exists(dir):
                continue
            fname = dir2 + ".tar.gz"
            url = "http://research-pub.gene.com/gmap/src"
            installs = "gmap gmapl get-genome gmapindex iit_store iit_get iit_dump gsnap gsnapl uniqscan uniqscanl snpindex cmetindex atoiindex sam_sort ../util/*.pl"
            cmd = "wget %s/%s; tar xvzf %s; cd %s; ./configure; make; cd src; cp %s ../../bin; cd ../.." % \
                (url, fname, fname, dir, installs)
        else:
            assert False
        print >> sys.stderr, cmd
        os.system(cmd)

    os.chdir("..")
            
    
if __name__ == "__main__":
    get_data()
