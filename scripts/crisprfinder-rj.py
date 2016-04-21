import sys
import logging
import argparse
import os
import cPickle
import gzip
import networkx as nx
from string import maketrans
import gc

def load_reads_one_cluster(seqfile):
    def genNextReadPair(fname):
        def genPair(lines):
            orig1 = ''
            orig2 = ''
            sp1 = lines[0].split()
            sp2 = lines[2].split()
            if len(sp1) == 2:
                orig1 = sp1[1]
            if len(sp2) == 2:
                orig2 = sp2[1]
            return sp1[0].strip(">"), lines[1].strip(), orig1, sp2[0].strip(">"), lines[3].strip(), orig2

        with open(fname) as f:
            lines = []
            for line in f:
                lines.append(line.strip())
                if len(lines) == 4:
                    yield genPair(lines)
                    lines = []

    # dict_reads: key: readid,  #value: seq
    dict_reads = {}
    dict_paired = {}
    dict_orig = {}
    for r1, seq1, orig1, r2, seq2, orig2 in genNextReadPair(seqfile):
        dict_reads[r1] =seq1
        dict_reads[r2] =seq2
        dict_paired[r1] = set([r2])
        dict_paired[r2] = set([r1])
        dict_orig[r1] = orig1
        dict_orig[r2] = orig2
    return dict_reads, dict_paired, dict_orig


def load_uniq_readids(fname):
    dict_uniq_ids = set()
    with open(fname) as f:
        for line in f:
            if line.startswith(">"):
                sp = line.strip().split()
                dict_uniq_ids.add(sp[0].lstrip('>'))
    return dict_uniq_ids

def enrich_paired_reads(dict_paired, dict_group2read, dict_read2group):
    for r in dict_paired:
        g = dict_read2group.get(r)
        if not g:  # not in the group
            continue
        to_update = {}
        for p in dict_paired[r]:  # for each paired read of r
            if not p in dict_paired:
                continue
            if r in dict_read2group:
                for r1 in dict_group2read[g]:
                    if not p in to_update:
                        to_update[p] = set()
                    else:
                        to_update[p].add(r1)
        for t in to_update:
            for s in to_update[t]:
                dict_paired[t].add(s)

def simplify_graph(g, confident_thresh, dup=False, printgraph=False):
    for e in g.selfloop_edges():
        g.remove_edge(e[0], e[1])
    for node in g.nodes():
        edges = g.in_edges(node)
        edges.extend(g.out_edges(node))
        plus = []
        minus = []
        for e in edges:
            data = g[e[0]][e[1]]
            if data[node] == '+':
                plus.append(e)
            else:
                minus.append(e)
        if not plus or not minus:
            continue
        if dup: # if a node can act as then end node for several branch, duplicate it.
            if len(plus) > len(minus):
                for e in minus:
                    if g.has_edge(e[0], e[1]):
                        g.remove_edge(e[0], e[1])
            elif len(plus) < len(minus):
                for e in plus:
                    if g.has_edge(e[0], e[1]):
                        g.remove_edge(e[0], e[1])
            else:  # dup
                for e in plus:
                    pass
        else:
            if len(plus) >= len(minus):
                for e in minus:
                    if g.has_edge(e[0], e[1]):
                        g.remove_edge(e[0], e[1])
            if len(plus) <= len(minus):
                for e in plus:
                    if g.has_edge(e[0], e[1]):
                        g.remove_edge(e[0], e[1])
    logging.info("Cleaning graph using confident edge")
    clean_graph_use_confident_edge(g, confident_thresh)
    printGraph(g, clusterid + '.simp1.pdf', sizeLimit, needprint=printGraph)
    logging.info("Checking cycles")
    if not nx.is_directed_acyclic_graph(g):
        #cycles = len(list(nx.simple_cycles(og)))
        logging.info("Cluster {0} is not DAG".format(clusterid))
        #logging.warning("Cycles: {0}".format(cycles))
        back_edges = find_back_edges(g)
        logging.info("Back edges in the graph:")
        for e in back_edges:
            logging.info("Remove back edge {0}".format(e))
            g.remove_edge(e[0], e[1])
    logging.info("Removing bubbles")
    remove_bubbles(g)

    logging.info("Removing tips")
    remove_out_tips(g)
    remove_in_tips(g)
    printGraph(g, clusterid + '.simp2.pdf', sizeLimit, needprint=printGraph)
    logging.info("Correcting directions in linear graphs")
    for c in nx.weakly_connected_component_subgraphs(g):
        if c.number_of_nodes() <= 2:
            continue
        isLinear, ends, source, sink= is_linear_graph(c)
        if isLinear:
            if sink==1 and source == 1:
                continue
            adjust_edge_di(g, c, ends[0], ends[1])

def clean_graph_use_confident_edge(g, confident_thresh):
    for n in g.nodes():
        has_confident_edge = False
        for e in g.out_edges([n], data=True):
            if e[2]['overlap'] >= confident_thresh:
                has_confident_edge = True
                break
        if has_confident_edge:
            for e in g.out_edges([n], data=True):
                if e[2]['overlap'] < confident_thresh:
                    logging.debug("Removing edge {0} -> {1} based on confident threshold.".format(e[0], e[1]))
                    g.remove_edge(e[0], e[1])

def remove_in_tips(g, keepNode=False, onlyRead=True):
    moreWork = False
    for n in g.nodes():
        if onlyRead:
            if g.node[n]['contig']:
                continue
        if len(g.node[n]['seq']) >= 120:
            continue
        if g.in_degree(n) != 0:
            continue
        if g.out_degree(n) > 1:
            continue
        for e in g.out_edges([n]):  # only loop once
            v = e[1]
            if g.in_degree(v) == 1:
                continue
            else:  # in degree > 1
                for e1 in g.in_edges([v]):
                    if e1[0] == n:
                        continue
                    if g.in_degree(e1[0]) > 0:
                        moreWork = True
                        logging.debug("Removing tips: {0}".format(n))
                        if keepNode:
                            g.remove_edge(e[0], e[1])
                        else:
                            g.remove_node(n)
                        break
    return moreWork

def remove_out_tips(g, keepNode=False, onlyRead=True):
    moreWork = False
    for n in g.nodes():
        if onlyRead:
            if g.node[n]['contig']:
                continue
        if len(g.node[n]['seq']) >= 120:
            continue
        if g.out_degree(n) != 0:
            continue
        if g.in_degree(n) > 1:
            continue
        for e in g.in_edges([n]):  # only loop once
            v = e[0]
            if g.out_degree(v) == 1:
                continue
            else:  # in degree > 1
                for e1 in g.out_edges([v]):
                    if e1[1] == n:
                        continue
                    if g.out_degree(e1[1]) > 0:
                        moreWork = True
                        if keepNode:
                            g.remove_edge(e[0], e[1])
                        else:
                            g.remove_node(n)
                        break
    return moreWork

def is_linear_graph(g):
    num_nodes = g.number_of_nodes()
    num_edges = g.number_of_edges()
    if num_nodes != num_edges + 1:
        return False, None, None, None
    source = 0
    sink = 0
    ends = []
    for n in g:
        if g.degree(n) <=2:
            if g.in_degree(n) == 0:
                source += 1
                if g.degree(n) == 1:
                    ends.insert(0, n)
            if g.out_degree(n) == 0:
                sink += 1
                if g.degree(n) == 1:
                    ends.append(n)
        elif g.degree(n) > 2:
            return False, None, None, None
    return True, ends, source, sink


def adjust_edge_di(g, c, start, end):
    cur = start
    preNode = None
    edge2remove = []
    while True:
        nei =  list(nx.all_neighbors(c, cur))
        if preNode:
            nei.remove(preNode)
        if not c.has_edge(cur, nei[0]):
            edge2remove.append((nei[0], cur))
        preNode = cur
        cur = nei[0]
        if nei[0] == end:
            break
    for s, e in edge2remove:
        orig1 = g[s][e]['label'].split("/")[1][1]
        orig2 = g[s][e]['label'].split("/")[1][0]
        label = g[s][e]['label'].split("/")[0] + "/" + orig2 + orig1
        g.add_edge(e, s,  {'weight':g[s][e]['weight'], 'overlap':g[s][e]['overlap'],
                           'label': label, e:orig1, s:orig2})
        g.remove_edge(s, e)

def remove_bubbles(g):
    node_2_remove = []
    for cp in nx.weakly_connected_component_subgraphs(g):
        if cp.number_of_nodes() > 10:
            p = longest_path_one_componet(cp)
            if len(p) < 5:
                return
            for n in p:
                g.node[n]['color'] = 'blue'

            for node in cp.nodes():
                if node in p: # node in path
                    continue
                nei = list(nx.all_neighbors(cp, node))
                if len(nei) == 1 and nei[0] in p:
                    i = p.index(nei[0])
                    if i >=2 and i <len(p)-2:
                        node_2_remove.append(node)
                    continue
                connect_to_path = False
                all_nei_connect_to_path = True
                for n in nei:
                    if n not in p: # node's neighbor is not in the path, then the neighbors of n should all be in the path, except node
                        nnei = nx.all_neighbors(cp, n)
                        for n1 in nnei:
                            if n1 != node:
                                if n1 not in p:
                                    all_nei_connect_to_path = False
                                    break
                                else:
                                    i = p.index(n1)
                                    if i<2 or i >= len(p)-2:
                                        all_nei_connect_to_path = False
                                        break
                    else: # node's neighbor is in the path, need to be in the middle of the path
                        i = p.index(n)
                        if i >=2 and i < len(p)-2:
                            connect_to_path = True
                if connect_to_path and all_nei_connect_to_path:
                    node_2_remove.append(node)
    for node in node_2_remove:
        if node in g.nodes():
            logging.debug("Removing bubble node {0}".format(node))
            g.remove_node(node)


def create_overlap_graphs(clusterid, dict_reads, dict_uniq_ids, graphfile,
                           dict_read2group, dict_group2read, dict_orig,
                           confident_thresh):
    logging.debug("Creating graph for cluster {0}".format(clusterid))
    og = nx.DiGraph()
    for r in dict_uniq_ids:
        og.add_node(r, group=set([r]), seq=dict_reads[r], contig=False, myreads=[r])
        if dict_orig[r] == 'r.n':
            og.node[r]['color'] = 'red'
        g = dict_read2group.get(r)
        if g:
            for r1 in dict_group2read[g]:
                if r1 != r:
                    og.node[r]['group'].add(r1)
    fr_edge = 0
    rf_edge = 0
    ff_edge = 0
    add_later = []
    with open(graphfile) as f:
        for line in f:
            sp = line.split()
            n1 = sp[0]
            n2 = sp[2]
            orig = sp[1]+sp[3]
            data = {'overlap':int(sp[4]), n1: sp[1], n2:sp[3], 'weight':int(sp[4]), 'label':sp[4]+"/"+orig}
            if orig == '++':
                ff_edge += 1
                og.add_edge(n1, n2, data)
            elif orig == '+-':
                add_later.append((n1, n2, data))
                fr_edge += 1
            elif orig == '-+':
                add_later.append((n1, n2, data))
                rf_edge += 1
    for n1, n2, data in add_later:
        inn1 = og.in_degree(n1)
        inn2 = og.in_degree(n2)
        outn1 = og.out_degree(n1)
        outn2 = og.out_degree(n2)
        if inn1 + outn2 < inn2+outn1:
            # need to switch orientation and label
            temp = data[n1]
            data[n1] = data[n2]
            data[n2] = temp
            data['label'] = data['label'].split("/")[0] + '/' + data[n2] + data[n1]
            og.add_edge(n2, n1, data)
        else:
            og.add_edge(n1, n2, data)
    if og.number_of_nodes() > 1000:
        sys.setrecursionlimit(og.number_of_nodes()+10)
    logging.debug("BEFORE: FF: {0}, FR: {1}, RF:{2}".format(ff_edge, fr_edge, rf_edge))
    printGraph(og, clusterid + '.beforesimp.pdf', sizeLimit)

    # fr_edge = 0
    # rf_edge = 0
    # ff_edge = 0
    # for edge in og.edges():
    #     orig = og[edge[0]][edge[1]]['label'].split("/")[1]
    #     if orig == '++':
    #         ff_edge += 1
    #     elif orig == '+-':
    #         fr_edge += 1
    #     elif orig == '-+':
    #         rf_edge += 1
    # logging.debug("AFTER: FF: {0}, FR: {1}, RF:{2}".format(ff_edge, fr_edge, rf_edge))
    return og


def load_group_info(groupfile):
    dict_read2group = {}
    dict_group2read = {}
    with open(groupfile) as f:
        for idx, line in enumerate(f):
            dict_group2read[idx] = set()
            sp = line.strip().split()
            for r in sp:
                dict_group2read[idx].add(r)
                dict_read2group[r] = idx
    return dict_read2group, dict_group2read

def reverse_comp(seq):
    ret = ''
    for c in seq:
        if c == 'A':
            ret = 'T' + ret
        elif c == 'T':
            ret = 'A' + ret
        elif c == 'C':
            ret = 'G' + ret
        else:
            ret = 'C' + ret
    return ret


def generate_seq_from_path(g, path, dict_seq):
    if not path:
        return ''
    outSeq = dict_seq[path[0]]
    if len(path) == 1:
        return outSeq
    data = g[path[0]][path[1]]
    if data[path[0]] == "-":
        outSeq = reverse_comp(outSeq)
    edges = zip(path[:-1], path[1:])
    for e in edges:
        data = g[e[0]][e[1]]
        curSeq = dict_seq[e[1]]
        if data[e[1]] == "-":
            curSeq = reverse_comp(curSeq)
        overlap = data['weight']
        curSeq = curSeq[overlap:]
        outSeq += curSeq
    return outSeq


def combine_line_nodes(g, dict_seq, dict_paired):
    def find_line_pre(g, node):
        cur = node
        ret = []
        while g.in_degree(cur)<=1 and g.out_degree(cur)<=1:
            ret.append(cur)
            pred =  g.predecessors(cur)
            if pred:
                cur = pred[0]
            else:
                break
        return ret[1:]

    def find_line_succ(g, node):
        cur = node
        ret = []
        while g.in_degree(cur)<=1 and g.out_degree(cur)<=1:
            ret.append(cur)
            succ = g.successors(cur)
            if succ:
                cur = succ[0]
            else:
                break
        return ret[1:]


    def combine_nodes(g, nodes2Combine, dict_paired, dict_seq):

        contig = generate_seq_from_path(g, nodes2Combine, dict_seq)

        # add a new node
        newName = nodes2Combine[0] + '/' + nodes2Combine[-1]
        dict_seq[newName] = contig
        g.add_node(newName, group=set(), seq=contig, contig=True, myreads=[])
        dict_paired[newName] = set()
        to_update = {}
        for n in nodes2Combine:
            for read in g.node[n]['myreads']:
                g.node[newName]['myreads'].append(read)
            for r in g.node[n]['group']:
                g.node[newName]['group'].add(r)
            for t in dict_paired[n]:
                dict_paired[newName].add(t)
                if t in dict_paired:
                    to_update[t] = newName
                    #dict_paired[t].add(newName)

        for t in to_update:
            dict_paired[t].add(to_update[t])
        # now connect the new node to the graph
        pred = g.predecessors(nodes2Combine[0])
        succ = g.successors(nodes2Combine[-1])
        for p in pred:
            w = g[p][nodes2Combine[0]]['weight']
            o = g[p][nodes2Combine[0]]['overlap']
            l = g[p][nodes2Combine[0]]['label']
            orig1 = g[p][nodes2Combine[0]][p]
            orig2 = '+'  # always '+' for the new combined node
            #g.add_edge(p, newName,  {'weight':w, 'overlap':o, 'label':l, p:orig1, newName: orig2})
            g.add_edge(p, newName,  {'weight':w, 'overlap':o, 'label':l.split("/")[0] + '/' + orig1 + '+', p:orig1, newName: orig2})
            logging.debug("Adding edge {0} to {1}, {2}".format(p, newName, g[p][nodes2Combine[0]]))
        for s in succ:
            w = g[nodes2Combine[-1]][s]['weight']
            o = g[nodes2Combine[-1]][s]['overlap']
            l = g[nodes2Combine[-1]][s]['label']
            orig1 = '+'
            orig2 = g[nodes2Combine[-1]][s][s]
            #g.add_edge(newName, s, {'weight':w, 'overlap':o, 'label':l, newName:orig1, s:orig2})
            g.add_edge(newName, s, {'weight':w, 'overlap':o, 'label':l.split("/")[0] + '/+' + orig2, newName:orig1, s:orig2})
            logging.debug("Adding edge {0} to {1}, {2}".format(newName, s, g[nodes2Combine[-1]][s]))
        # now remove old node
        for n in nodes2Combine:
            g.remove_node(n)

        return newName

    processed = set()
    moreWork = True
    while moreWork:
        moreWork = False
        for node in g.nodes():
            if node in processed:
                continue
            if g.degree(node)>2:
                processed.add(node)
            pred = find_line_pre(g, node)
            succ = find_line_succ(g, node)
            if not pred and not succ:  # only one node, no need to combine.
                processed.add(node)
                continue
            # need to combine
            moreWork = True
            nodes2Combine = pred[::-1]
            nodes2Combine.append(node)
            nodes2Combine.extend(succ)
            newNode = combine_nodes(g, nodes2Combine, dict_paired, dict_seq)
            processed.add(newNode)


def check_paired_info(g, node1, node2, dict_paired):
    paired1 = dict_paired[node1]
    num_paired =  paired1.intersection(g.node[node2]['group'])
    return num_paired

def check_branch(g, dict_seq, dict_paired, threshold=2):
    def is_branch_node(g, node):
        if g.in_degree(node) > 1 or g.out_degree(node) > 1:
            return True
        return False

    #componets = sorted(nx.connected_components(g), key = len, reverse=True)
    sortedNodes = nx.topological_sort(g)
    moreWork = False
    for node in sortedNodes:
        is_branch = is_branch_node(g, node)
        if is_branch:
            logging.debug("Node {0} is a branching node.".format(node))
            pred = g.predecessors(node)
            succ = g.successors(node)
            max_paired = 0
            max_pred = ''
            max_succ = ''
            for p in pred:
                for s in succ:
                    paired = check_paired_info(g, p, s, dict_paired)
                    num_paired = len(paired)
                    logging.debug("Number of paired end read between {0} and {1}: {2}.".format(p, s, paired))
                    if num_paired > max_paired:
                        max_paired = num_paired
                        max_pred = p
                        max_succ = s
            if max_paired >= threshold:
                moreWork = True
                logging.info("Removing branching edge based on paired end info")
                for p in pred:
                    if p != max_pred:
                        logging.info("Removing edge {0} -> {1}".format(p, node))
                        g.remove_edge(p, node)
                for s in succ:
                    if s != max_succ:
                        logging.info("Removing edge {0} -> {1}".format(node, s))
                        g.remove_edge(node, s)
                g.edge[max_pred][node]['color'] = 'red'
                g.edge[node][max_succ]['color'] = 'red'
    return moreWork


def combine_until_nochange(g, dict_seq, dict_paired, needprint=False):
    cnt = 0

    moreWork = True
    while moreWork:
        logging.info("========================== Round {0} =======================".format(cnt))
        printGraph(og, clusterid + '.' + str(cnt) + '.before.pdf', sizeLimit, needprint=needprint)
        combine_line_nodes(g, dict_seq, dict_paired)
        printGraph(og, clusterid + '.' + str(cnt) + '.after.pdf', sizeLimit, needprint=needprint)
        moreWork = check_branch(g, dict_seq, dict_paired)
        remove_in_tips(g)
        remove_in_tips(g)
        printGraph(og, clusterid + '.' + str(cnt) + '.check.pdf', sizeLimit, needprint=needprint)
        cnt += 1
        #test_node_orig(og)
    logging.info("========================== DONE COMBINING =======================".format(cnt))


def test_node_orig(og):
    for node in og.nodes():
        orig = set()
        for nei in nx.all_neighbors(og, node):
            if og.has_edge(nei, node):
                orig.add(og.edge[nei][node][node])
            if og.has_edge(node, nei):
                orig.add(og.edge[node][nei][node])
        print(node, orig)
        if len(orig) > 1:
            logging.error("Wrong, node {0} has wrong direction...".format(node))
        for nei in nx.all_neighbors(og, node):
            if og.has_edge(nei, node):
                print(og[nei][node])
            else:
                print(og[node][nei])
    print("===================")

def printGraph(g, outputname, sizeLimit, needprint=False):
    if not needprint:
        return
    dict_mapping = {}
    for node in g.nodes():
        dict_mapping[node] =  node + "/s" +str(len(g.node[node]['seq']))
    g1 = nx.relabel_nodes(g, dict_mapping)
    if g1.number_of_nodes() > sizeLimit:
        return
    A=nx.to_agraph(g1)
    A.draw(outputname,prog='dot')



class DFSResult:
    def __init__(self):
        self.parent = {}
        self.start_time = {}
        self.finish_time = {}
        self.edges = {} # Edge classification for directed graph.
        self.order = []
        self.t = 0

def dfs(g):
    results = DFSResult()
    for vertex in g.nodes():
        if vertex not in results.parent:
            dfs_visit(g, vertex, results)
    return results

def dfs_visit(g, v, results, parent = None):
    results.parent[v] = parent
    results.t += 1
    results.start_time[v] = results.t
    if parent:
        results.edges[(parent, v)] = 'tree'

    for n in g.neighbors(v):
        if n not in results.parent: # n is not visited.
            dfs_visit(g, n, results, v)
        elif n not in results.finish_time:
            results.edges[(v, n)] = 'back'
        elif results.start_time[v] < results.start_time[n]:
            results.edges[(v, n)] = 'forward'
        else:
            results.edges[(v, n)] = 'cross'
    results.t += 1
    results.finish_time[v] = results.t
    results.order.append(v)

def find_back_edges(g):
    ret = []
    result = dfs(g)
    for edge in result.edges:
        if result.edges[edge] == 'back':
            ret.append(edge)
    return ret



def longestPathAll2(g, dict_seq, nodePath=True):
    allpaths = []
    if is_linear_graph(g)[0]:
        p = get_path_linear_graph(g)
        return [p]
    while g.number_of_nodes() > 0:
        #print(str(g.number_of_nodes()))
        #sys.stderr.write('\n')
        paths = []
        if nodePath:
            paths = longestPathNode(g, dict_seq)
        else:
            paths = longestPath(g, dict_seq)
        if paths:
            allpaths.extend(paths)
            for p in paths:
                for n in p:
                    g.remove_node(n)
    return allpaths


def get_path_linear_graph(g):
    start = None
    end = None
    for n in g.nodes():
        if g.in_degree(n) == 0:
            start = n
        elif g.out_degree(n) == 0:
            end = n
    cur = start
    path = [start]
    while True:
        o = g.out_edges([cur])
        if o:
            path.append(o[0][1])
            cur = path[-1]
        else:
            return path


def longestPath(g, dict_seq):
    paths = []
    if g.number_of_nodes() == 0:
        return paths
    if g.number_of_nodes() == 1:
        paths.append(g.nodes())
        return paths
    if is_linear_graph(g)[0]:
        p = get_path_linear_graph(g)
        return [p]
    for c in nx.weakly_connected_component_subgraphs(g):
        if c.number_of_nodes() == 1:
            paths.append(c.nodes())
            continue
        dist = {}
        for node in nx.topological_sort(c):
            pairs = [(dist[v][0]+len(dict_seq[node])- g[v][node]['weight'], v) for v in c.pred[node]]
            if pairs:
                dist[node] = max(pairs)
            else:
                dist[node] = (len(dict_seq[node]), node)
        node, (length, _) = max(dist.items(), key=lambda x:x[1])
        path = []
        while length > len(dict_seq[node]):
            path.append(node)
            length, node = dist[node]
        paths.append(list(reversed(path)))
    return paths



def longestPathNode(g, dict_seq):
    paths = []
    if g.number_of_nodes() == 0:
        return paths
    if g.number_of_nodes() == 1:
        paths.append(g.nodes())
        return paths
    if is_linear_graph(g)[0]:
        p = get_path_linear_graph(g)
        return [p]
    for c in nx.weakly_connected_component_subgraphs(g):
        if c.number_of_nodes() == 1:
            paths.append(c.nodes())
            continue
        dist = {}
        for node in nx.topological_sort(c):
            pairs = [(dist[v][0]+1, v) for v in c.pred[node]]
            if pairs:
                dist[node] = max(pairs)
            else:
                dist[node] = (0, node)
        node, (length, _) = max(dist.items(), key=lambda x:x[1])
        path = []
        while length > 0:
            path.append(node)
            length, node = dist[node]
        paths.append(list(reversed(path)))
    return paths


def longest_path_one_componet(g):
    paths = []
    if g.number_of_nodes() == 0:
        return paths
    if g.number_of_nodes() == 1:
        paths.append(g.nodes())
        return paths
    if is_linear_graph(g)[0]:
        p = get_path_linear_graph(g)
        return p
    dist = {}
    for node in nx.topological_sort(g):
        pairs = [(dist[v][0]+1, v) for v in g.pred[node]]
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    node, (length, _) = max(dist.items(), key=lambda x:x[1])
    path = []
    while length > 0:
        path.append(node)
        length, node = dist[node]
    return path


def get_read_in_path(g, path):
    all_reads = []
    for n in path:
        for r in g.node[n]['myreads']:
            all_reads.append(r)
    return all_reads


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='CRISPR gene assembler.'
    )
    parser.add_argument('-d','--debug',
                        help='Print lots of debugging statements',
                        action="store_const",dest="loglevel",const=logging.DEBUG,
                        default=logging.INFO
                    )
    parser.add_argument('-p', '--drawOverlap', help='Save overlap graphs to pdf files', action="store_true")
    parser.add_argument('-s', '--savePath', help='Save path corresponds to contigs', action="store_true")
    parser.add_argument('-o', '--outputFolder', help='Save overlap graphs to pdf files', default="./")
    parser.add_argument('-n', '--namePrefix', help='Prefix for naming output files.')
    parser.add_argument('-l', '--sizeLimit', type=int, help='The size limit when drawing graph. If the number of node in the graph is larger than this, then do not draw.')

    parser.add_argument('clusterfolder', metavar='fasta-folder', help='Folder that contains fasta reads for each cluster.')
    parser.add_argument('groupFile', help='Group file.')
    parser.add_argument('rjfolder', help='The folder contains all the readjoiner results')
    parser.add_argument('initialOverlap', metavar='initial-threshold', help='Initial overlap threshold.', type=int)
    parser.add_argument('overlap', metavar='confident-threshold', help='Confident edge overlap threshold.', type=int)
    args = parser.parse_args()
    #FORMAT = '%(asctime)-15s %(message)s'
    logger = logging.getLogger(__name__)
    FORMAT = '%(asctime)s - %(module)-16s - %(levelname)s - %(message)s'
    logging.basicConfig(level=args.loglevel, format=FORMAT)

    sizeLimit = 1000
    if args.sizeLimit:
        sizeLimit = args.sizeLimit

    overlapthresh = args.initialOverlap
    confidentThresh = args.overlap
    logging.info("Initial overlap threshold in the asqg file: {0}".format(overlapthresh))
    logging.info("Confident edge overlap threshold: {0}".format(confidentThresh))
    if confidentThresh <= overlapthresh:
        logging.error("Confident threhold should be larger than the initial overlap threshold")
        sys.exit(1)
    if not os.path.exists(args.outputFolder):
        os.mkdir(args.outputFolder)
    pefiles = []
    for root, dirs, names in os.walk(args.clusterfolder):
        for name in names:
            if name.endswith(".pe.fa"):
                pefiles.append(os.path.join(root, name))
        break
    logging.info("Number of clusters: {0}".format(len(pefiles)))
    # reads that belong to the same group are loaded.
    dict_read2group, dict_group2read = load_group_info(args.groupFile)
    dict_olg =  {}
    for pefile in pefiles:
        logging.info("\n\n===========================================".format(pefile))
        logging.info("Processing cluster {0}".format(pefile))
        dict_reads, dict_paired, dict_orig = load_reads_one_cluster(pefile)
        enrich_paired_reads(dict_paired, dict_group2read, dict_read2group)
        clusterid = os.path.basename(pefile).split(".")[0]
        graphfile = os.path.join(args.rjfolder, clusterid+".graph.rename.txt")
        undup_fasta = os.path.join(args.rjfolder, clusterid+".reads.filtered.fa")
        logging.info("Number of reads in this cluster: {0}".format(len(dict_reads)))

        if not os.path.exists(graphfile):
            logging.warning("Graph file for cluster{0} does not exist. Read joiner failed on this cluster.".format(clusterid))
            continue
        dict_uniq_ids = load_uniq_readids(undup_fasta)
        logging.info("Number of uniq reads in this cluster: {0}".format(len(dict_uniq_ids)))
        logging.info("Creating overlap graph for cluster {0}".format(pefile))
        og = create_overlap_graphs(clusterid, dict_reads, dict_uniq_ids, graphfile, dict_read2group, dict_group2read, dict_orig, args.overlap)
        simplify_graph(og, confident_thresh, args.drawOverlap)
        printGraph(og, clusterid + '.pdf', sizeLimit, needprint=args.drawOverlap)
        logging.info("Combining node....")
        #combine_until_nochange(og, dict_reads, dict_paired, needprint=args.drawOverlap)
        printGraph(og, clusterid + '.combined.pdf', sizeLimit, needprint=args.drawOverlap)
        cp_index = 0
        all_reads_file = open(clusterid + ".allread.txt", 'w')
        contigfile = open(clusterid + ".contig.fa", 'w')
        for cp in nx.weakly_connected_component_subgraphs(og):
            #print(cp.nodes())
            cp1 = cp.copy()  # copy because we will remove node..
            pathes = longestPathAll2(cp1, dict_reads)
            contig_index = 0
            for p in pathes:
                seq = generate_seq_from_path(cp, p, dict_reads)
                if len(seq) < 105:
                    continue
                seqid = clusterid + '_' + str(cp_index) + "_" + str(contig_index) + '_' + 'l' + str(len(seq))
                contigfile.write('>' + seqid + '\n')
                contigfile.write(seq + '\n')
                all_reads = get_read_in_path(cp, p)
                all_reads_file.write(">" + seqid)
                for r in all_reads:
                    all_reads_file.write('\t' + r)
                all_reads_file.write('\n')
                contig_index += 1
            cp_index += 1
        logging.info("Nodes: {0}, edges: {1}".format(og.number_of_nodes(), og.number_of_edges()))
        gc.collect()
        logging.info("===========================================".format(pefile))
