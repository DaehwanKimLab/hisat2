# HISAT-bisulfite

We use HISAT2 for graph representation and alignment, which is currently the most practical and quickest program available. We refer to hisat-bisulfite-top as our top directory where all of our programs are located. hisat-bisulfite-top is a place holder that you can change to whatever name you’d like to use.
In order to install HISAT2, please run the following commands.

$ git clone https://github.com/infphilo/hisat2 hisat-bisulfite-top

$ cd hisat-bisulfite-top

hisat-bisulfite-top$ git checkout bisulfite

$ make hisat2-align-s hisat2-build-s hisat2-inspect-s

Add the above directory (hisat-bisulfite-top) to your PATH environment variable (e.g. ~/.bashrc) to make the binaries we just built above and other python scripts available everywhere:

export PATH=hisat-bisulfite-top:hisat-bisulfite-top/hisatgenotype_scripts:$PATH

export PYTHONPATH=hisat-bisulfite-top/hisatgenotype_modules:$PYTHONPATH

$ source ~/.bashrc

$ mkdir bisulfite-analysis

$ cd bisulfite-analysis

Additional program requirements:
 
 SAMtools (version 0.1.19 or later)
 
 bigWigToBedGraph at http://hgdownload.soe.ucsc.edu/admin/exe/

Build epigenome by incorporating known CpG sites into the human reference genome (GRCh38)

$ hisatbisulfite_build_epigenome.py

Test the followign command:
$ hisat-bisulfite-top/evaluation/bisulfite//hisatbisulfite_eval.py --region 10:1000000-2000000

HISAT-bisulfite project is led by Bongsoo Park (genomicspark@gmail.com) and Daehwan Kim (infphilo@gmail.com)
