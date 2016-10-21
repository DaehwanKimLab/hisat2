#!/usr/bin/env python

import sys
import math

class Node:
    # Initialize
    def __init__(self, left, seq, var):
        self.next = [] # list of next nodes

        self.nodeID = 0 # Node ID

        self.left = left # starting position
        
        self.seq = seq # sequence that node represents
        self.var = var # how sequence is related to backbone
        self.cov = []  # coverage along the sequence by reads

        assert 'I' not in seq
        self.right = self.left + len(seq) - 1

        
    # Check how compatible allele is in regard to read or pair
    def compatible_with_rnode(self, rnode):
        assert rnode.left + len(rnode.seq) <= len(self.seq)
        score = 0
        for i in range(len(rnode.seq)):
            allele_bp = self.seq[rnode.left + i]
            read_bp = rnode.seq[i]
            if allele_bp == read_bp:
                score += 1
        return float(score) / len(rnode.seq)

    
    # Check how nodes overlap with each other
    def overlap_with(self, other):
        assert self.left <= other.left
        if self.right < other.left:
            return -1, -1
        
        seq = self.seq.replace('D', '')
        other_seq = other.seq.replace('D', '')
        assert 'I' not in seq and 'I' not in other_seq

        # Currently slow
        max_mm = 2
        for i in range(len(seq)):
            temp_mm = 0
            for j in range(len(other_seq)):
                if i + j >= len(seq):
                    break
                if seq[i+j] != other_seq[j]:
                    temp_mm += 1
                    if temp_mm > max_mm:
                        break
            if temp_mm <= max_mm:
                return i, len(seq) - i

        return -1, -1

    
    # Combine two nodes
    def combine_with(self, other, at):
        seq = self.seq.replace('D', '')
        overlap = len(seq) - at
        assert overlap > 0
        assert self.left <= other.left
        if self.right >= other.right:
            return

        for i in range(len(other.seq)):
            if overlap <= 0:
                break
            c = other.seq[i]
            if c == 'D':
                continue
            assert c in "ACGT"
            overlap -= 1

        assert i < len(other.seq)
        self.seq += other.seq[i:]
        self.var += other.var[i:]
        assert len(self.seq) == len(self.var)

        self.right = self.left + len(self.seq) - 1

        
    # Display node information
    def print_info(self):
        print >> sys.stderr, "Pos: [%d, %d]" % (self.left, self.right)
        print >> sys.stderr, "\t", self.seq
        prev_var = ""
        for var_i in range(len(self.var)):
            var = self.var[var_i]
            if var != "" and var != "unknown" and var != prev_var:
                print >> sys.stderr, "\t%d: %s" % (var_i, var),
            prev_var = var
        print >> sys.stderr

                
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
        if id in self.nodes:
            # print >> sys.stderr, "Warning) multi-mapped read:", id
            return
        assert id not in self.nodes
        self.nodes[id] = node

    def add_edge(self):
        None

    # Display graph information
    def print_info(self): 
        print >> sys.stderr, "Backbone len: %d" % len(self.backbone)
        print >> sys.stderr, "\t%s" % self.backbone

    def assemble(self):
        if len(self.nodes) <= 0:
            return
        nodes = [[id, node.left, node.right] for id, node in self.nodes.items()]
        def node_cmp(a, b):
            return a[1] - b[1]
        nodes = sorted(nodes, cmp=node_cmp)

        # DK - debugging purposes
        queue = [nodes[0]]; nodes.pop(0)
        while len(queue) > 0:
            node_id, node_left, node_right = queue.pop()
            node = self.nodes[node_id]

            node_i = 0
            del_nodes = []
            while node_i < len(nodes):
                node2_id, node2_left, node2_right = nodes[node_i]
                node2 = self.nodes[node2_id]
                at, overlap = node.overlap_with(node2)
                if overlap >= 80:
                    # DK - debugging purposes
                    print >> sys.stderr, node_id; node.print_info()
                    print >> sys.stderr, node2_id; node2.print_info()

                    node.combine_with(node2, at)
                    del_nodes.append(node_i)

                    print >> sys.stderr, "at %d, overlap with %d" % (at, overlap)
                    print >> sys.stderr, "combined:", node_id; node.print_info()
                else:
                    # DK - debugging purposes
                    print >> sys.stderr, node_id; node.print_info()
                    print >> sys.stderr, node2_id; node2.print_info()
                    print >> sys.stderr, "at %d, overlap with %d" % (at, overlap)
                    sys.exit(1)
                    
                    node_i += 1
                    continue

                node_i += 1

            # DK - debugging purposes
            print >> sys.stderr, node_id; node.print_info()
            # sys.exit(1)
            
            new_nodes = []
            for node_i in range(len(nodes)):
                if node_i not in del_nodes:
                    new_nodes.append(nodes[node_i])
            nodes = new_nodes
            if len(nodes) <= 0:
                break
            queue.append(nodes[0]); nodes.pop(0)
        
        
    # Draw graph
    def draw(self):
        htmlDraw = HtmlDraw("assembly_graph")
        htmlDraw.write_html_css(self.width, self.height)
        htmlDraw.start_js()
        htmlDraw.draw_smile()

        # DK - debugging purposes
        nodes = [[id, node.left] for id, node in self.nodes.items()]
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
        print >> html_file, r'<title>HISAT-genotyping HLA</title>'
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
       
