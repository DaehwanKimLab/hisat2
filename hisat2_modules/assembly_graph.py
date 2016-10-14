#!/usr/bin/env python

import sys
import math

try:
    from pyx import *
except ImportError:
    print >> sys.stderr, "pyx module is not available on the system"

class Node:
    # Initialize
    def __init__(self, pos, seq, var):
        self.next = [] # list of next nodes

        self.pos = pos # starting position
        
        self.seq = seq # sequence that node represents
        self.var = var # how sequence is related to backbone
        self.cov = []  # coverage along the sequence by reads
        
    # Check how compatible allele is with read or pair
    def compatible_with_rnode(self, rnode):
        assert rnode.pos + len(rnode.seq) <= len(self.seq)
        score = 0
        for i in range(len(rnode.seq)):
            allele_bp = self.seq[rnode.pos + i]
            read_bp = rnode.seq[i]
            if allele_bp == read_bp:
                score += 1
        return float(score) / len(rnode.seq)
 
    # Display node information
    def help(self):
        print >> sys.stderr, "Pos: %d" % self.pos
        print >> sys.stderr, "\t", self.seq
        print >> sys.stderr, "\t", self.var

                
class Graph:
    def __init__(self, backbone):
        # self.head = Node()

        self.backbone = backbone # backbone sequence

    # Display graph information
    def help(self):
        print >> sys.stderr, "Backbone len: %d" % len(self.backbone)
        print >> sys.stderr, "\t%s" % self.backbone

    # Draw graph
    def draw(self):
        c = canvas.canvas()
        c.text(0, 0, "Hello, world!")
        c.stroke(path.line(0, 0, 2, 0))
        c.stroke(path.line(0, 0, 3, 0))
        c.stroke(path.rect(0, 1, 1, 1))
        c.fill(path.circle(2.5, 1.5, 0.5))
        c.writePDFfile("temp.pdf")
    
