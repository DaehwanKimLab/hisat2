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
        # print >> sys.stderr, "\t", self.seq
        print >> sys.stderr, "\t", self.var

                
class Graph:
    def __init__(self, backbone):
        # self.head = Node()
        self.backbone = backbone # backbone sequence

        self.nodes = {}
        self.edges = {}

        self.left_margin = 20
        self.right_margin = 20
        self.top_margin = 20
        self.bottom_margin = 20

        self.scale = 5
        self.width = len(self.backbone) * self.scale + self.left_margin + self.right_margin
        self.height = 1000 * self.scale

    def add_node(self, id, node):
        assert id not in self.nodes
        self.nodes[id] = node

    def add_edge(self):
        None

    # Display graph information
    def help(self): 
        print >> sys.stderr, "Backbone len: %d" % len(self.backbone)
        print >> sys.stderr, "\t%s" % self.backbone

    # Draw graph
    def draw(self):
        htmlDraw = HtmlDraw("assembly_graph")
        htmlDraw.write_html_css(self.width, self.height)
        htmlDraw.start_js()
        htmlDraw.draw_smile()

        # DK - debugging purposes
        nodes = [[id, node.pos] for id, node in self.nodes.items()]
        def node_cmp(a, b):
            return a[1] - b[1]
        nodes = sorted(nodes, cmp=node_cmp)
        # Choose font
        js_file = htmlDraw.js_file
        print >> js_file, r'context.font = "10px Serif";'
        for id, pos in nodes:
            read_id, mate = id.split('|')[:2]
            mate = mate.split('_')[0]

            # Draw node
            print >> js_file, r'context.fillStyle = "yellow";'
            print >> js_file, r'context.beginPath();'
            print >> js_file, r'context.arc(%d, %d, %d, 0, 2*Math.PI);' % \
                (self.left_margin + pos * self.scale, self.top_margin + pos * self.scale / 4, 2 * self.scale)
            print >> js_file, r'context.closePath();'
            print >> js_file, r'context.fill();'
            print >> js_file, r'context.lineWidth = 2;'
            print >> js_file, r'context.stroke();'
            print >> js_file, r'context.fillStyle = "black";'

            # Draw label
            print >> js_file, r'context.fillText("%s %s", %d, %d);' % \
                (read_id, mate, self.left_margin + pos * self.scale, self.top_margin + pos * self.scale / 4)
               
        htmlDraw.end_js()

        
class HtmlDraw:
    def __init__(self, base_fname):
        self.base_fname = base_fname

        
    def write_html_css(self, width = 2000, height = 1000):
        base_fname = self.base_fname
        html_file = open("%s.html" % base_fname, 'w')
        print >> html_file, r'<!DOCTYPE html>'
        print >> html_file, r'<html>'
        print >> html_file, r'<head>'
        print >> html_file, r'<title>Smiley Face</title>'
        print >> html_file, r'<link rel="stylesheet" type="text/css" href="%s.css"/>' % base_fname
        print >> html_file, r'</head>'
        print >> html_file, r'<body>'
        print >> html_file, r'<canvas id="a" width="%d" height="%d">' % (width, height)
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

        
    def start_js(self):
        self.js_file = open("%s.js" % self.base_fname, 'w')
        print >> self.js_file, r'var a_canvas = document.getElementById("a");'
        print >> self.js_file, r'var context = a_canvas.getContext("2d");'

        
    def end_js(self):
        self.js_file.close()

        
    def draw_smile(self):
        js_file = self.js_file
        
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
        print >> js_file, r'context.fillText("Hello, World!", 15, 175);'
       
