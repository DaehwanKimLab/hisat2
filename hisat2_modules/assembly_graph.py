#!/usr/bin/env python

import sys
import math

class Node:
    # Initialize
    def __init__(self, pos, seq, var):
        self.next = [] # list of next nodes

        self.nodeID = 0 # Node ID

        self.pos = pos # starting position
        
        self.seq = seq # sequence that node represents
        self.var = var # how sequence is related to backbone
        self.cov = []  # coverage along the sequence by reads
        
    # Check how compatible allele is in regard to read or pair
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
        htmlDraw = HtmlDraw()
        htmlDraw.write_html_css("assembly_graph")

        
class HtmlDraw:
    def __init__(self):
        None
            
    def write_html_css(self, base_fname):
        html_file = open("%s.html" % base_fname, 'w')
        print >> html_file, r'<!DOCTYPE html>'
        print >> html_file, r'<html>'
        print >> html_file, r'<head>'
        print >> html_file, r'<title>Smiley Face</title>'
        print >> html_file, r'<link rel="stylesheet" type="text/css" href="%s.css"/>' % base_fname
        print >> html_file, r'</head>'
        print >> html_file, r'<body>'
        print >> html_file, r'<canvas id="a" width="200" height="200">'
        print >> html_file, r'This text is displayed if your browser does not support HTML5 Canvas.'
        print >> html_file, r'</canvas>'
        print >> html_file, r'<script type="text/javascript" src="%s.js"></script>' % base_fname
        print >> html_file, r'</body>'
        print >> html_file, r'</html>'
        html_file.close()

        css_file = open("%s.css" % base_fname, 'w')
        print >> css_file, r'canvas {'
        print >> css_file, r'border: 1px dotted black;'
        print >> css_file, r'}'
        css_file.close()

        js_file = open("%s.js" % base_fname, 'w')
        print >> js_file, r'var a_canvas = document.getElementById("a");'
        print >> js_file, r'var context = a_canvas.getContext("2d");'

        # Draw the face
        print >> js_file, r'context.fillStyle = "yellow";'
        print >> js_file, r'context.beginPath();'
        print >> js_file, r'context.arc(95, 85, 40, 0, 2*Math.PI);'
        print >> js_file, r'context.closePath();'
        print >> js_file, r'context.fill();'
        print >> js_file, r'context.lineWidth = 2;'
        print >> js_file, r'context.stroke();'
        print >> js_file, r'context.fillStyle = "black";'
        
        # Draw the left eye
        print >> js_file, r'context.beginPath();'
        print >> js_file, r'context.arc(75, 75, 5, 0, 2*Math.PI);'
        print >> js_file, r'context.closePath();'
        print >> js_file, r'context.fill();'

        # Draw the right eye
        print >> js_file, r'context.beginPath();'
        print >> js_file, r'context.arc(114, 75, 5, 0, 2*Math.PI);'
        print >> js_file, r'context.closePath();'
        print >> js_file, r'context.fill();'

        # Draw the mouth
        print >> js_file, r'context.beginPath();'
        print >> js_file, r'context.arc(95, 90, 26, Math.PI, 2*Math.PI, true);'
        print >> js_file, r'context.closePath();'
        print >> js_file, r'context.fill();'

        # Write "Hello, World!"
        print >> js_file, r'context.font = "30px Garamond";'
        print >> js_file, r'context.fillText("Hello, World!",15,175);'
        js_file.close()
    
