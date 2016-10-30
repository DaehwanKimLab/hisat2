#!/usr/bin/env python

import sys
import math
import random
from copy import deepcopy


#
def get_major_nt(nt_dic):
    nt = ''
    max_count = 0
    for tmp_nt, tmp_count in nt_dic.items():
        assert nt in "ACGTD"
        if tmp_count > max_count:
            max_count = tmp_count
            nt = tmp_nt
    assert nt in "ACGTD"
    return nt                


#
def match_score(nt_dic1, nt_dic2):
    total1, total2 = sum(nt_dic1.values()) * 2.0, sum(nt_dic2.values()) * 2.0
    best = 0.0
    for nt in "ACGT":
        if nt not in nt_dic1 or nt not in nt_dic2:
            continue
        tmp_best = nt_dic1[nt] / total1 + nt_dic2[nt] / total2
        if tmp_best > best:
            best = tmp_best
    return best


#
def get_ungapped_seq(seq):
    ungapped_seq = []
    for nt_dic in seq:
        if get_major_nt(nt_dic) == 'D':
            continue
        ungapped_seq.append(nt_dic)
    return ungapped_seq


class Node:
    # Initialize
    def __init__(self, left, seq, var):
        self.next = [] # list of next nodes

        self.nodeID = 0 # Node ID

        self.left = left # starting position

        # sequence that node represents
        assert 'I' not in seq
        self.seq = []
        for nt in seq:
            assert nt in "ACGTD"
            self.seq.append({nt : 1})

        # how sequence is related to backbone
        self.var = []
        for v in var:
            self.var.append(set([v]) if v != '' else set())

        assert len(self.seq) == len(self.var)
        self.right = self.left + len(seq) - 1

        
    # Check how compatible allele is in regard to read or pair
    def compatible_with_rnode(self, rnode):
        assert False
        assert rnode.left + len(rnode.seq) <= len(self.seq)
        score = 0
        for i in range(len(rnode.seq)):
            allele_bp = self.seq[rnode.left + i]
            read_bp = rnode.seq[i]
            if allele_bp == read_bp:
                score += 1
        return float(score) / len(rnode.seq)


    # Check how nodes overlap with each other without considering deletions
    def overlap_with(self, other):
        assert self.left <= other.left
        if self.right < other.left:
            return -1, -1

        seq, other_seq = get_ungapped_seq(self.seq), get_ungapped_seq(other.seq)
        max_mm = 2.0

        # DK - debugging purposes
        max_mm = 0.1
        
        for i in range(len(seq)):
            tmp_mm = 0.0
            for j in range(len(other_seq)):
                if i + j >= len(seq):
                    break
                other_nt_dic = other.seq[j]
                mismatch = 1.0 - match_score(seq[i+j], other_seq[j])
                assert mismatch >= 0.0
                tmp_mm += mismatch
                if tmp_mm > max_mm:
                    break
            if tmp_mm <= max_mm:
                return i, len(seq) - i

        return -1, -1

    
    # Combine two nodes with considering deletions
    def combine_with(self, other):
        assert self.left <= other.left
        if self.right >= other.right:
            return

        # Merge two sequences
        assert len(other.seq) > 0 and 'D' not in other.seq[0].keys()
        i, j = other.left - self.left, 0
        new_seq, new_var = self.seq[:i], self.var[:i]
        while i < len(self.seq) and j < len(other.seq):
            nt_dic, nt_var = self.seq[i], self.var[i]
            nt_dic2, nt_var2 = other.seq[j], other.var[j]
            new_seq.append(nt_dic)
            new_var.append(nt_var | nt_var2)
            for nt, count in nt_dic2.items():
                if nt in nt_dic:
                    nt_dic[nt] += count
                else:
                    nt_dic[nt] = count
            i += 1
            j += 1
            
        # Append the rest of the other sequence
        assert i == len(self.seq)
        new_seq += other.seq[j:]
        new_var += other.var[j:]

        self.seq, self.var = new_seq, new_var
        assert len(self.seq) == len(self.var)
        self.right = self.left + len(self.seq) - 1


    # Return the length of the ungapped sequence
    def ungapped_length(self):
        return len(get_ungapped_seq(self.seq))

    # Get variants
    def get_vars(self, Vars):
        vars = []
        for var_i in range(len(self.var)):
            for var in self.var[var_i]:
                assert var in Vars
                if len(vars) > 0 and var == vars[-1]:
                    continue
                type, pos, data = Vars[var]
                nt = get_major_nt(self.seq[var_i])
                if data == nt or \
                   type == "deletion" and nt == 'D':
                    vars.append(var)                    
        return vars

        
    # Display node information
    def print_info(self):
        print >> sys.stderr, "Pos: [%d, %d]" % (self.left, self.right)
        seq = ""
        for nt_dic in self.seq:
            seq += get_major_nt(nt_dic)
        print >> sys.stderr, "\t", seq
        prev_var = ""
        for var_i in range(len(self.var)):
            var = self.var[var_i]
            var = '-'.join(list(var))
            if var != "" and var != prev_var:
                print >> sys.stderr, "\t%d: %s" % (var_i, var), self.seq[var_i],
            prev_var = var
        print

                
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

        self.scalex, self.scaley = 5, 2
        self.width = len(self.backbone) * self.scalex + self.left_margin + self.right_margin
        self.height = 2000 * self.scaley


    # Add node, which is an alignment w.r.t. the reference
    def add_node(self, id, node):
        if id in self.nodes:
            # print >> sys.stderr, "Warning) multi-mapped read:", id
            return
        assert id not in self.nodes
        self.nodes[id] = node


    # Generate edges based on the overlapping information between nodes
    def generate_edges(self, overlap_pct = 0.8):
        assert len(self.nodes) > 0
        nodes = [[id, node.left, node.right] for id, node in self.nodes.items()]
        def node_cmp(a, b):
            return a[1] - b[1]
        nodes = sorted(nodes, cmp=node_cmp)

        self.from_node, self.to_node = {}, {}
        for i in range(len(nodes)):
            id1, left1, right1 = nodes[i]
            node1 = self.nodes[id1]
            for j in range(i + 1, len(nodes)):
                id2, left2, right2 = nodes[j]
                if right1 < left2:
                    break
                node2 = self.nodes[id2]
                at, overlap = node1.overlap_with(node2)
                if overlap < node1.ungapped_length() * overlap_pct:
                    continue
                if id1 not in self.to_node:
                    self.to_node[id1] = [[id2, at]]
                else:
                    self.to_node[id1].append([id2, at])
                if id2 not in self.from_node:
                    self.from_node[id2] = [[id1, -at]]
                else:
                    self.from_node[id2].append([id1, -at])

        
    # Display graph information
    def print_info(self): 
        print >> sys.stderr, "Backbone len: %d" % len(self.backbone)
        print >> sys.stderr, "\t%s" % self.backbone


    # Reduce graph
    def reduce(self):
        nodes = [[id, node.left, node.right] for id, node in self.nodes.items()]
        def node_cmp(a, b):
            if a[1] != b[1]:
                return a[1] - b[1]
            return a[2] - b[2]
        nodes = sorted(nodes, cmp=node_cmp)

        # node id to number
        id_to_num, num_to_id, node_len = {}, [], []
        for node_id, _, _ in nodes:
            id_to_num[node_id] = len(id_to_num)
            num_to_id.append(node_id)
            node_len.append(self.nodes[node_id].ungapped_length())

        VACANT, INPLAY, ELIMINATED = range(3)
        EDGE_UNUSED, EDGE_REDUNDANT, EDGE_UNITIG, EDGE_REDUCED, EDGE_USED = range(5)
        
        edge_range = []
        edge_list = []
        contain_list = {}
        contained_list = [False] * len(nodes)
        for node_id, to_node_ids in self.to_node.items():
            id1 = id_to_num[node_id]
            for to_node_id, at in to_node_ids:
                id2 = id_to_num[to_node_id]
                if at == 0:
                    contained_list[id2] = True
                    if id1 not in contain_list:
                        contain_list[id1] = [id2]
                    else:
                        contain_list[id1].append(id2)
                    
                edge_list.append([id1, id2,  at, EDGE_UNUSED])
                edge_list.append([id2, id1, -at, EDGE_UNUSED])

        temp_edge_list = []
        for edge in edge_list:
            if contained_list[edge[0]] or contained_list[edge[1]]:
                continue
            temp_edge_list.append(edge)

        edge_list = temp_edge_list
        edge_list.sort()

        num_node = len(id_to_num)
        count = 0
        begin, end = 0, 0
        for i in range(num_node):
            while count < len(edge_list) and edge_list[count][0] == i:
                end = end + 1
                count = count + 1
            edge_range.append([begin, end])
            edge_list[begin:end] = sorted(edge_list[begin:end], key=lambda x:abs(x[2]))
            begin = end

        mark = [VACANT] * num_node
        # Mark redundant edges v ==> x when v ==> w ==> x
        for v in range(num_node):
            for edge in edge_list[edge_range[v][0]:edge_range[v][1]]:
                w = edge[1]
                mark[w] = INPLAY

            for edge in edge_list[edge_range[v][0]:edge_range[v][1]]:
                v, w, begin = edge[:-1]
                if begin < 0 or mark[w] != INPLAY:
                    continue
                
                for edge2 in edge_list[edge_range[w][0]:edge_range[w][1]]:
                    _, x, begin2 = edge2[:-1]
                    if begin2 < 0 or mark[x] != INPLAY:
                        continue

                    mark[x] = ELIMINATED

            for edge in edge_list[edge_range[v][0]:edge_range[v][1]]:
                w = edge[1]
                if mark[w] == ELIMINATED:
                    edge[3] = EDGE_REDUNDANT

                    for opposite_edge in edge_list[edge_range[w][0]:edge_range[w][1]]:
                        if opposite_edge[1] == v:
                            opposite_edge[3] = EDGE_REDUNDANT
                            break

                mark[w] = VACANT

        edge_list.sort()

        # Unitig
        unitig_list = [[[i, 0]] for i in range(num_node)]
        vertex_included = [i for i in range(num_node)]

        def belong(v):
            v_temp = v
            while True:
                v_parent = vertex_included[v_temp]
                if v_parent == v_temp:
                    vertex_included[v] = v_parent
                    return v_parent

                v_temp = v_parent


        # Merge two unitigs by edge
        def merge_unitig(edge):
            orig_v, orig_w = edge[:2]
            v, w = belong(orig_v), belong(orig_w)
            unitig_v = unitig_list[v][:]
            unitig_w = unitig_list[w][:]

            if unitig_v[0][0] == orig_v:
                v_index = 0
            else:
                v_index = len(unitig_v) - 1

            if unitig_w[0][0] == orig_w:
                w_index = 0
            else:
                w_index = len(unitig_w) - 1

            v_info, w_info = unitig_v[v_index], unitig_w[w_index]

            # DK - debugging purposes
            if v_info[0] != orig_v or w_info[0] != orig_w:
                print >> sys.stderr, "v_info:", orig_v, "w_info:", orig_w
                print >> sys.stderr, "unitig_v:", unitig_v, "unitig_w:", unitig_w
            
            assert v_info[0] == orig_v and w_info[0] == orig_w

            if edge[2] >= 0:
                assert unitig_v[-1][0] == orig_v and unitig_w[0][0] == orig_w
                unitig_w[0][1] = edge[2]
                unitig_v.extend(unitig_w)
            else:
                assert unitig_v[0][0] == orig_v and unitig_w[-1][0] == orig_w
                unitig_v[0][1] = -edge[2]
                unitig_w.extend(unitig_v)
                unitig_v = unitig_w

            # orig_v ==> orig_w
            edge[3] = EDGE_UNITIG
            for edge2 in edge_list[edge_range[orig_w][0]:edge_range[orig_w][1]]:
                # orig_w ==> orig_v
                if edge2[1] == orig_v:
                    edge2[3] = EDGE_UNITIG
                    break

            vertex_included[w] = v
            unitig_list[w] = []
            unitig_list[v] = unitig_v

            return unitig_v

        def find_edges(v):
            result_edge_list = []
            v = belong(v)

            v_list = [unitig_list[v][0][0]]
            if len(unitig_list[v]) > 1:
                v_list.append(unitig_list[v][-1][0])

            result_list = []
            for v_temp in v_list:
                result = []
                for edge in edge_list[edge_range[v_temp][0]:edge_range[v_temp][1]]:
                    if edge[3] == EDGE_UNUSED:
                        result.append(edge)
                result_list.append(result)
                
            return result_list

        def merge_unitigs():
            for orig_v in range(num_node):
                v = belong(orig_v)
                v_edge_llist = find_edges(v)
                for v_edge_list in v_edge_llist:
                    if len(v_edge_list) > 2:
                        continue

                    if len(v_edge_list) == 2:
                        if (v_edge_list[0][2] > 0 and v_edge_list[1][2] > 0) or \
                           (v_edge_list[0][2] < 0 and v_edge_list[1][2] < 0):
                            continue
                    
                    for e in v_edge_list:
                        temp_edge_list = []
                        for e2 in edge_list[edge_range[e[1]][0]:edge_range[e[1]][1]]:
                            if e2[3] == EDGE_UNUSED:
                                temp_edge_list.append(e2)

                        left_count, right_count = 0, 0
                        for e2 in temp_edge_list:
                            if e2[2] < 0:
                                left_count += 1
                            else:
                                assert e2[2] > 0
                                right_count += 1
                        if (e[2] < 0 and right_count > 1) or \
                           (e[2] > 0 and left_count > 1):
                            continue
                            
                        merge_unitig(e)

        merge_unitigs()

        # Incorporate the nodes that are previously inside or identical to other nodes
        new_unitig_list = []
        for unitig in unitig_list:
            if len(unitig) <= 0:
                continue
            if len(unitig) == 1 and contained_list[unitig[0][0]]:
                continue
            new_unitig = []
            for id, at in unitig:
                new_unitig.append([id, at])
                if id in contain_list:
                    for id2 in contain_list[id]:
                        new_unitig.append([id2, 0])
            new_unitig_list.append(new_unitig)
        unitig_list = new_unitig_list


        # DK - debugging purposes
        """
        for i in range(len(unitig_list)):
            if i > 20:
                break
            unitig = unitig_list[i]
            print >> sys.stderr, "unitig", i,
            for node in unitig:
                num, at = node
                id = num_to_id[num]
                id = id.split('_')[0]
                print >> sys.stderr, "[%s, %d]," % (id, nodes[num][1]),
            print >> sys.stderr

        debug_id = 13
        print >> sys.stderr, "%d:" % debug_id,
        for edge in edge_list[edge_range[debug_id][0]:edge_range[debug_id][1]]:
            f, t, at = edge[:3]
            id = num_to_id[f].split('_')[0]
            id2 = num_to_id[t].split('_')[0]
            print >> sys.stderr, "[%s ==> %s at %d]" % (id, id2, at),
        print >> sys.stderr
        """
        

        # Perform the assembly of unitigs into nodes
        new_nodes = {}
        for unitig in unitig_list:
            assert len(unitig) > 0
            node_num = unitig[0][0]
            node_id = nodes[node_num][0]
            node = self.nodes[node_id]
            node.merged_nodes = [node_id.split('|')[0]]
            for node_num2, _ in unitig[1:]:
                node_id2 = nodes[node_num2][0]
                node2 = self.nodes[node_id2]
                node.combine_with(node2)
                node.merged_nodes.append(node_id2.split('|')[0])
            new_nodes[node_id] = node

        self.nodes = new_nodes
        self.generate_edges(0.1) # 10% overlap

        
    # Merge graph using mate pairs
    def assemble_with_mates(self):
        assert len(self.nodes) > 0
        nodes = [[id, node.left, node.right] for id, node in self.nodes.items()]
        def node_cmp(a, b):
            return a[1] - b[1]
        nodes = sorted(nodes, cmp=node_cmp)
        final_nodes = []
        queue = [[self.nodes[nodes[0][0]], nodes[0][0], nodes[0][0]]]; nodes.pop(0)
        while len(queue) > 0 and len(final_nodes) < 10:
            node, node_id, node_id_last = queue.pop()
            if node_id_last not in self.to_node or len(self.to_node[node_id_last]) <= 0:
                final_nodes.append([node, node_id, node_id_last])

                # DK - debugging purposes
                # print >> sys.stderr, "DK1:", node_id, node_id_last
                
                continue
            elif len(self.to_node[node_id_last]) == 1:
                node_id2 = self.to_node[node_id_last][0][0]
                node2 = self.nodes[node_id2]
                new_node = deepcopy(node)
                new_node.combine_with(node2)
                new_node.merged_nodes += node2.merged_nodes
                queue.append([new_node, node_id, node_id2])
            else:
                # DK - debugging purposes
                # print >> sys.stderr, "#merged nodes:", node.merged_nodes
                # print >> sys.stderr, "#edge list:", self.to_node[node_id_last]
                
                node_ids2, max_mate = [], 0
                for tmp_node_id, _ in self.to_node[node_id_last]:
                    tmp_node = self.nodes[tmp_node_id]
                    tmp_mate = len(set(node.merged_nodes) & set(tmp_node.merged_nodes))

                    # DK - debugging purposes
                    # print >> sys.stderr, "tmp_mate:", tmp_mate
                    # print >> sys.stderr, tmp_node_id, tmp_node.merged_nodes

                    if max_mate < tmp_mate:
                        max_mate = tmp_mate
                        node_ids2 = [tmp_node_id]
                    elif max_mate == tmp_mate:
                        node_ids2.append(tmp_node_id)

                if len(node_ids2) > 0:
                    for node_id2 in node_ids2:
                        node2 = self.nodes[node_id2]
                        new_node = deepcopy(node)
                        new_node.combine_with(node2)
                        new_node.merged_nodes += node2.merged_nodes
                        queue.append([new_node, node_id, node_id2])
                else:
                    final_nodes.append([node, node_id, node_id_last])
                    
                    # DK - debugging purposes
                    # print >> sys.stderr, "DK2:", node_id, node_id_last, self.to_node[node_id_last]


        return final_nodes        
    

    # Assemble
    def assemble(self):
        if len(self.nodes) <= 0:
            return
        nodes = [[id, node.left, node.right] for id, node in self.nodes.items()]
        def node_cmp(a, b):
            return a[1] - b[1]
        nodes = sorted(nodes, cmp=node_cmp)

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
                    print node_id; node.print_info()
                    print node2_id; node2.print_info()

                    node.combine_with(node2)
                    del_nodes.append(node_i)

                    print "at %d, overlap with %d" % (at, overlap)
                    print "combined:", node_id; node.print_info()
                else:
                    # DK - debugging purposes
                    print node_id; node.print_info()
                    print node2_id; node2.print_info()
                    print "at %d, overlap with %d" % (at, overlap)
                    sys.exit(1)
                    
                    node_i += 1
                    continue

                node_i += 1

            # DK - debugging purposes
            print node_id; node.print_info()
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
    #   Top left as (0, 0) and Bottom right as (width, height)
    def draw(self, second_allele = sys.maxint):
        assert len(self.nodes) > 0
        nodes = [[id, node.left, node.right] for id, node in self.nodes.items()]
        def node_cmp(a, b):
            return a[1] - b[1]
        nodes = sorted(nodes, cmp=node_cmp)

        # display space
        dspace = [[[0, 1000]]] * (nodes[-1][2] + 100)
        def get_dspace(left, right, height):
            assert left < len(dspace) and right < len(dspace)
            range1 = dspace[left]
            for range2 in dspace[left + 1:right + 1]:
                new_range = []
                # sub range
                for t1, b1 in range1:
                    for t2, b2 in range2:
                        if b1 < t2:
                            break
                        if b2 < t1:
                            continue
                        t, b = max(t1, t2), min(b1, b2)
                        if b - t >= height:
                            new_range.append([t, b])

                range1 = new_range
            if len(range1) <= 0:
                return 0

            t, b = range1[0]
            assert b - t >= height
            b = t + height
            for i in range(left, right+1):
                range1 = dspace[i]
                range2 = []
                found = False
                for j in range(len(range1)):
                    t2, b2 = range1[j]
                    if t2 <= t and b <= b2:
                        found = True
                        if t2 < t:
                            range2.append([t2, t])
                        if b < b2:
                            range2.append([b, b2])
                    else:
                        range2.append([t2, b2])
                dspace[i] = range2
                assert found
            return t

        def get_x(x):
            return self.left_margin + x * self.scalex

        def get_y(y):
            return self.top_margin + y * self.scaley

        # Get scalar
        def get_sx(x):
            return x * self.scalex

        def get_sy(y):
            return y * self.scaley

        htmlDraw = HtmlDraw("assembly_graph")
        htmlDraw.write_html_css(self.width, self.height)
        htmlDraw.start_js()
        # htmlDraw.draw_smile()
        js_file = htmlDraw.js_file

        # Choose font
        print >> js_file, r'ctx.font = "12px Serif";'

        # Draw vertical dotted lines at every 100nt and thick lines at every 500nt
        print >> js_file, r'ctx.fillStyle = "gray";'
        for pos in range(100, nodes[-1][2], 100):
            if pos % 500 == 0:
                print >> js_file, r'ctx.setLineDash([]);'
                print >> js_file, r'ctx.lineWidth = 1;'
            else:
                print >> js_file, r'ctx.setLineDash([5, 15]);'
                print >> js_file, r'ctx.lineWidth = 0.2;'

            print >> js_file, r'ctx.beginPath();'
            print >> js_file, r'ctx.moveTo(%d, %d);' % \
                (get_x(pos), self.top_margin)
            print >> js_file, r'ctx.lineTo(%d, %d);' % \
                (get_x(pos), self.height)
            print >> js_file, r'ctx.stroke();'

            # Draw label
            print >> js_file, r'ctx.fillStyle = "blue";'
            print >> js_file, r'ctx.fillText("%d", %d, %d);' % \
                (pos, get_x(pos+2), get_y(200))

        print >> js_file, r'ctx.setLineDash([]);'

        # Draw nodes
        node_to_y = {}
        for id, left, right in nodes:
            read_id, mate = id.split('|')[:2]
            mate = mate.split('_')[0]

            # Get y position
            y = get_dspace(left, right, 14)
            node_to_y[id] = y

            if int(read_id) < second_allele:
                color = "yellow"
            else:
                color = "green"

            # Draw node
            print >> js_file, r'ctx.beginPath();'
            print >> js_file, r'ctx.rect(%d, %d, %d, %d);' % \
                (get_x(left), get_y(y), get_x(right) - get_x(left), get_sy(10))
            print >> js_file, r'ctx.fillStyle = "%s";' % color
            print >> js_file, r'ctx.fill();'
            print >> js_file, r'ctx.lineWidth = 2;'
            print >> js_file, r'ctx.strokeStyle = "black";'
            print >> js_file, r'ctx.stroke();'

            # Draw label
            print >> js_file, r'ctx.fillStyle = "blue";'
            print >> js_file, r'ctx.fillText("%s %s", %d, %d);' % \
                (read_id, mate, get_x(left + 2), get_y(y + 7))

        # Draw edges
        print >> js_file, r'ctx.lineWidth = 1;'
        line_colors = ["red", "black", "blue"]
        for node_id, to_node_ids in self.to_node.items():
            node = self.nodes[node_id]
            node_x = (get_x(node.left) + get_x(node.right)) / 2
            node_y = get_y(node_to_y[node_id] + 5)
            print >> js_file, r'ctx.strokeStyle = "%s";' % \
                line_colors[random.randrange(len(line_colors))]
            for to_node_id, _ in to_node_ids:
                to_node = self.nodes[to_node_id]
                to_node_x = (get_x(to_node.left) + get_x(to_node.right) + (random.random() * 10 - 5)) / 2
                to_node_y = get_y(node_to_y[to_node_id] + 5)

                jitter1, jitter2 = (random.random() * 10 - 5), (random.random() * 10 - 5)
                jitter1, jitter2 = get_sx(jitter1), get_sx(jitter2)

                print >> js_file, r'ctx.beginPath();'
                print >> js_file, r'ctx.moveTo(%d, %d);' % (node_x + jitter1, node_y)
                print >> js_file, r'ctx.lineTo(%d, %d);' % (to_node_x + jitter2, to_node_y)
                print >> js_file, r'ctx.stroke();'
               
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
        print >> self.js_file, r'var ctx = a_canvas.getContext("2d");'

        
    def end_js(self):
        self.js_file.close()

        
    def draw_smile(self):
        js_file = self.js_file
        
        # Draw the face
        print >> js_file, r'ctx.fillStyle = "yellow";'
        print >> js_file, r'ctx.beginPath();'
        print >> js_file, r'ctx.arc(95, 85, 40, 0, 2*Math.PI);'
        print >> js_file, r'ctx.closePath();'
        print >> js_file, r'ctx.fill();'
        print >> js_file, r'ctx.lineWidth = 2;'
        print >> js_file, r'ctx.stroke();'
        print >> js_file, r'ctx.fillStyle = "black";'
        
        # Draw the left eye
        print >> js_file, r'ctx.beginPath();'
        print >> js_file, r'ctx.arc(75, 75, 5, 0, 2*Math.PI);'
        print >> js_file, r'ctx.closePath();'
        print >> js_file, r'ctx.fill();'

        # Draw the right eye
        print >> js_file, r'ctx.beginPath();'
        print >> js_file, r'ctx.arc(114, 75, 5, 0, 2*Math.PI);'
        print >> js_file, r'ctx.closePath();'
        print >> js_file, r'ctx.fill();'

        # Draw the mouth
        print >> js_file, r'ctx.beginPath();'
        print >> js_file, r'ctx.arc(95, 90, 26, Math.PI, 2*Math.PI, true);'
        print >> js_file, r'ctx.closePath();'
        print >> js_file, r'ctx.fill();'

        # Write "Hello, World!"
        print >> js_file, r'ctx.font = "30px Garamond";'
        print >> js_file, r'ctx.fillText("Hello, World!", 15, 175);'
       
