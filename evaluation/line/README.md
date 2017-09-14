# HISAT-rep

We use HISAT2 for graph representation and alignment, which is currently the most practical and quickest program available. We refer to hisat-rep-top as our top directory where all of our programs are located. hisat-rep-top is a place holder that you can change to whatever name you’d like to use.
In order to install HISAT2, please run the following commands.

$ git clone https://github.com/infphilo/hisat2 hisat-rep-top

$ cd hisat-rep-top

hisat-rep-top$ git checkout line

$ make hisat2-align-s hisat2-build-s hisat2-inspect-s

Add the above directory (hisat-rep-top) to your PATH environment variable (e.g. ~/.bashrc) to make the binaries we just built above and other python scripts available everywhere:

export PATH=hisat-rep-top:hisat-rep-top/hisatgenotype_scripts:$PATH

export PYTHONPATH=hisat-rep-top/hisatgenotype_modules:$PYTHONPATH

$ source ~/.bashrc

$ cd hisat-rep-top/evaluation/line

Additional program requirements:
 
 SAMtools (version 0.1.19 or later)
 
 ClustalW2 at http://www.clustal.org/clustal2/#Download

Test the followign command:
$ hisat-rep-top/evaluation/line/hisatgenotype_locus.py --base rep

HISAT-rep project is led by Hyunmin Kim (human.gim@gmail.com) and Daehwan Kim (infphilo@gmail.com)
