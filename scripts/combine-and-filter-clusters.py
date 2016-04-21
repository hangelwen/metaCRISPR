import sys
import networkx as nx
import cPickle
import collections
import gc
import logging
import argparse

#adapted from http://courses.csail.mit.edu/6.006/fall11/rec/rec14.pdf
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

def get_which_cluster_for_reads(clusterfile):
    dict_read2cluster = {}
    with open(clusterfile) as f:
        for line in f:
            sp = line.split()
            if len(sp) < 3:
                continue
            for r in sp[1:]:
                dict_read2cluster[r] = sp[0]
    return dict_read2cluster


def get_single_read_cluster(fastaname, dict_read2cluster):
    set_single_read = set()
    with open(fastaname) as f:
        for line in f:
            if line.startswith(">"):
                readid = line.split()[0].strip(">")
                if not readid in dict_read2cluster:
                    set_single_read.add(readid)
    return set_single_read


def enrich_cluster_using_sam(samname, dict_read2cluster):
    single_read_in_sam = {}
    with open(samname) as f:
        for line in f:
            if not line.startswith("@"):
                sp = line.split()
                if sp[2] in dict_read2cluster:
                    dict_read2cluster[sp[0]] = dict_read2cluster[sp[2]]
                else:  # sp[2] is a single read
                    single_read_in_sam[sp[0]] = sp[2]
    return single_read_in_sam


def cluster_to_reads(dict_read2cluster):
    dict_cluster2reads = {}
    for r in dict_read2cluster:
        if dict_read2cluster[r] not in dict_cluster2reads:
            dict_cluster2reads[dict_read2cluster[r]] = set()
        dict_cluster2reads[dict_read2cluster[r]].add(r)
    return dict_cluster2reads



def update_cluster(dict_read2cluster, dict_cluster2reads, old1, old2, new):
    dict_cluster2reads[new] = dict_cluster2reads[old1]
    dict_cluster2reads[new].update(dict_cluster2reads[old2])
    for r in dict_cluster2reads[old1]:
        dict_read2cluster[r] = new
    for r in dict_cluster2reads[old2]:
        dict_read2cluster[r] = new
    dict_cluster2reads.pop(old1)
    dict_cluster2reads.pop(old2)


def read_connected_components(graphfile):
    dict_read_2_overlap = {}
    with open(graphfile) as f:
        for line in f:
            sp = line.split()
            if sp[0] not in dict_read_2_overlap:
                dict_read_2_overlap[sp[0]] = set()
            dict_read_2_overlap[sp[0]].add(line.strip())
            if sp[1] not in dict_read_2_overlap:
                dict_read_2_overlap[sp[1]] = set()
            dict_read_2_overlap[sp[1]].add(line.strip())
    return dict_read_2_overlap


def buildGraphFromReads(reads, dict_read_2_overlap):
    set_other_reads = set()
    g = nx.DiGraph()
    for r in reads:
        if r not in dict_read_2_overlap:  # contained reads
            set_other_reads.add(r)
            continue
        for line in dict_read_2_overlap[r]:
            sp = line.split()
            if not g.has_edge(sp[0], sp[2]):
                g.add_edge(sp[0], sp[2], {sp[0]: sp[1], sp[2]: sp[3]}, weight = int(sp[4]), label=sp[1]+sp[3])
    return g, set_other_reads


def dump_graphs_to_file(dict_cluster2reads, dict_read_2_overlap, outname):
    fout = open(outname, 'w')
    dict_graphs = {}
    for c in dict_cluster2reads:
        g = buildGraphFromReads(dict_cluster2reads[c], dict_read_2_overlap)
        dict_graphs[c] = g
    cPickle.dump(dict_graphs, fout, protocol=2)


def generate_initial_contigs_for_cluster(g, dict_seq, other_reads, kmer):
    dict_cluster_seq = {}
    for r in other_reads:
        dict_cluster_seq[r] = dict_seq[r]
    for cp in nx.weakly_connected_component_subgraphs(g):  # cp is copied from g
        pathes = longestPathAll2(cp, dict_all_reads)
        for p in pathes:
            seq = generate_seq_from_path(g, p, dict_all_reads)
            #TODO: generate_seq_from_path
            if seq.find(kmer) == -1:
                seq = reverse_comp(seq)
                p = p[::-1]
            seqid = "|".join(p)
            dict_cluster_seq[seqid] = seq
    return dict_cluster_seq



def filter_clusters(dict_cluster2reads, dict_read_2_overlap, dict_all_reads, dict_DRReads, dict_dr_reads_in_cluster,
                    disgard, onlyinset, minGraphSize=10):
    cluster_cnt = 0
    for c in dict_cluster2reads:
        gc.collect()
        if onlyinset:
            if c not in onlyinset:
                continue
        logging.info("Processing cluster {0}, number of reads {1}".format(c, len(dict_cluster2reads[c])))
        g, other_reads = buildGraphFromReads(dict_cluster2reads[c], dict_read_2_overlap)
        if g.number_of_nodes() < minGraphSize:
            logging.info("Cluster {0} has less than {1} reads. Filtered".format(c, minGraphSize))
            continue
        simplify_graph(g)
        if not nx.is_directed_acyclic_graph(g): # remove back edges to make it a DAG
            sys.setrecursionlimit(g.number_of_nodes()+10)
            back_edges = find_back_edges(g)
            logging.info("Back edges in the graph:")
            for e in back_edges:
                logging.info("Remove back edge {0}".format(e))
                g.remove_edge(e[0], e[1])

        seqGood = 0
        repeats = []
        seqs = []
        counter = collections.Counter()
        has_bad_cc = False
        for cp in nx.weakly_connected_component_subgraphs(g):  # cp is copied from g
            r = float(cp.number_of_edges())/cp.number_of_nodes()
            if r > 2:
                max_degree = max(g.degree().items(), key=lambda x:x[1])
                logging.info("Cluster {0} has a complex component. {1}, {2}. Filtered".format(c, r, max_degree))
                has_bad_cc = True
                break
        if has_bad_cc:
            continue
        for cp in nx.weakly_connected_component_subgraphs(g):  # cp is copied from g
            pathes = longestPath(cp, dict_all_reads)
            for p in pathes:
                seq = generate_seq_from_path(g, p, dict_all_reads)
                seqs.append(seq)
                repeat= has_good_repeat(seq)
                if repeat:
                    print(repeat)
                    seqGood += 1
                    counter.update(gen_kmer(repeat[0]))
                    counter.update(gen_kmer(repeat[1]))

        if len(counter) == 0:
            logging.info("Cluster {0} has no repeat structure. filtered".format(c))
            continue
        ratio = seqGood/float(len(seqs))
        mc = counter.most_common(5)
        hasgoodKmer = False
        ratio2 = 0
        theGoodKmer = ''
        for e in mc:
            kmer = e[0]
            kmerGood = 0
            for s in seqs:
                all_occ = find_all_occ(s, kmer)
                if len(all_occ) >= 1:
                    kmerGood +=1
                    continue
                all_occ = find_all_occ(s, reverse_comp(kmer))
                if len(all_occ) >= 1:
                    kmerGood +=1
                    continue
            ratio2 = kmerGood/float(len(seqs))
            if ratio2 >=0.5:
                hasgoodKmer = True
                theGoodKmer = kmer
                break
        if hasgoodKmer:
            outnameContig = 'cluster-' + str(cluster_cnt) + ".contig.1.fa"
            outnameAll = 'cluster-' + str(cluster_cnt) + ".pe.fa"
            logging.info("Cluster {0} has good contigs, ratio {1}, ratio2 {2}. File {3} {4} Pass".format(c, ratio, ratio2, outnameContig, outnameAll))
            write_one_cluster_with_all_pe_reads(dict_cluster2reads, dict_DRReads, theGoodKmer, c, outnameAll, disgard)
            cluster_cnt += 1
        else:
                logging.info("Cluster {0} has only a few good contigs, ratio {1}, ratio2 {2}. Filtered".format(c, ratio, ratio2))
                continue

def get_pe_count_clusters(readsBaseNames, dict_read2cluster, dict_DRReads):
    dict_pe_count = {}  # pe reads between two clusters
    dict_pe_count_in_cluster = {}  # pe reads in a cluster
    for r in readsBaseNames:
        r1 = r
        r2 = r
        if r in dict_DRReads:
            if '1' in dict_DRReads[r]:
                r1 = r1 +'.r.1'
            else:
                r1 = r1 + '.1'
            if '2' in dict_DRReads[r]:
                r2 = r2 +'.r.2'
            else:
                r2 = r2 + '.2'
        c1 = dict_read2cluster.get(r1)
        c2 = dict_read2cluster.get(r2)
        if c1 and c2:
            if c1 != c2:
                joint_name = get_joint_name(c1, c2)
                if not joint_name in dict_pe_count:
                    dict_pe_count[joint_name] = 0
                dict_pe_count[joint_name] += 1
            else:
                #logging.info("Two ends of {0} are in the same cluster.\n".format(r))
                if c1 not in dict_pe_count_in_cluster:
                    dict_pe_count_in_cluster[c1] = 0
                dict_pe_count_in_cluster[c1] += 1
        else:
            #logging.info("One end of {0} is not in any cluster.\n".format(r))
            pass
    return dict_pe_count, dict_pe_count_in_cluster


def get_paired_read(r, dict_DRReads):
    bn = get_basename(r)
    e = '1'
    if r[-1] == '1':
        e = '2'
    if bn in dict_DRReads:
        if e in dict_DRReads[bn]:
            return bn + '.r.' + e
    return bn + "." + e


def get_read_basename_and_dr_reads(fastafullname):
    readbasenames = set()
    drreads = {}  # map basename to dr reads
    dict_read_repeatinfo = {}
    with open(fastafullname) as f:
        for line in f:
            if line.startswith(">"):
                sp = line.split()
                namepart = sp[0].strip().split('.')
                bn = namepart[0].lstrip('>')
                which = namepart[-1]
                if namepart[1] == 'r':
                    if bn not in drreads:
                        drreads[bn] = set()
                    drreads[bn].add(which)
                readbasenames.add(bn)
                if len(sp)==2:
                    readname = sp[0].lstrip('>')
                    dict_read_repeatinfo[readname] = {'type':1,
                                                      'repeat':[],
                                                      'pos':[]}
                    sp1 = sp[1].split('/')
                    if len(sp1) == 4:
                        dict_read_repeatinfo[readname]['type'] = 2
                        dict_read_repeatinfo[readname]['repeat'].append(sp1[2])
                        dict_read_repeatinfo[readname]['repeat'].append(sp1[3])
                        dict_read_repeatinfo[readname]['pos'].append(int(sp1[0]))
                        dict_read_repeatinfo[readname]['pos'].append(int(sp1[1]))
                    else:
                        dict_read_repeatinfo[readname]['type'] = 1
                        dict_read_repeatinfo[readname]['repeat'].append(sp1[1])
                        dict_read_repeatinfo[readname]['pos'].append(int(sp1[0]))

    return readbasenames, drreads, dict_read_repeatinfo


def get_basename(r):
    return r.split('.')[0]


def read_full_fasta(fastafull):
    dict_reads = {}
    readid = ''
    seq = ''
    with open(fastafull) as f:
        for line in f:
            if line.startswith(">"):
                if readid:
                    dict_reads[readid] = seq
                readid = line.split()[0].lstrip(">")
            else:
                seq = line.strip()
        dict_reads[readid] = seq
    return dict_reads

def load_overlaps(fname):
    dict_overlaps = {}
    with open(fname) as f:
        for line in f:
            sp = line.split()
            dict_overlaps[(sp[0], sp[2])] = int(sp[4])
            dict_overlaps[(sp[2], sp[4])] = int(sp[4])
    return dict_overlaps


def classify_single_reads(set_single_read, dict_read2cluster, dict_DRReads, output):
    num_disgard = 0
    num_classifed = 0
    for r in set_single_read:
        pr = get_paired_read(r, dict_DRReads)
        if pr in dict_read2cluster:
            dict_read2cluster[r] = dict_read2cluster[pr]
            #dict_cluster2reads[dict_read2cluster[r]] = pr
            #logging.info("Read {0} is classified to cluster {1}.\n".format(r, dict_read2cluster[pr]))
            num_classifed += 1
        else:
            #logging.info("Read {0} can not be classified to any cluster using paired-end info, disgarding.\n".format(r))
            num_disgard += 1
            if output:
                print(r)
    logging.info("{0} single reads can not be classified/disgarded".format(num_disgard))
    logging.info("{0} single reads are classified using PE info".format(num_classifed))


def get_joint_name(c1, c2):
    if c1 < c2:
        return c1 + "@" + c2
    else:
        return c2 + "@" + c1

def get_drread_count_in_clusters(dict_cluster2reads, dict_DRReads):
    dict_dr_reads_in_cluster = {}
    for c in dict_cluster2reads:
        #dict_dr_reads_in_cluster[c] = set()
        processed = set()
        cnt = 0
        for r in dict_cluster2reads[c]:
            b = get_basename(r)
            if b in processed:
                continue
            processed.add(b)
            pr = get_paired_read(r, dict_DRReads)
            if r.find('r') != -1:
                cnt += 1
            if pr.find('r') != -1:
                cnt += 1
        dict_dr_reads_in_cluster[c] = cnt
    return dict_dr_reads_in_cluster


def write_one_cluster_with_all_pe_reads(dict_cluster2reads, dict_DRReads,
                                        kmer, c, fname, disgard):

    drreads = dict_dr_reads_in_cluster[c]
    print("Full PE file: {0}\t{1}".format(fname, c))

    processed = set()
    readnames = []
    reads = []
    for r in dict_cluster2reads[c]:
        b = get_basename(r)
        if b in processed:
            continue
        processed.add(b)
        pr = get_paired_read(r, dict_DRReads)
        r1 = r
        r2 = pr
        if not r1.endswith('.1'):
            r1 ,r2 = r2, r1
        seq1 = dict_all_reads[r1]
        if seq1.find(kmer) == -1:
            seq1 = reverse_comp(seq1)
            r1 = r1 + '\t' + 'r'
            if seq1.find(kmer) == -1:
                r1 = r1 + '.n'

        seq2 = dict_all_reads[r2]
        if seq2.find(kmer) == -1:
            seq2 = reverse_comp(seq2)
            r2 = r2 + '\t' + 'r'
            if seq2.find(kmer) == -1:
                r2 = r2 + '.n'

        readnames.append(r1)
        reads.append(seq1)
        readnames.append(r2)
        reads.append(seq2)

    curKmerSize = 20
    while True:
        need_2_correct = set()
        good = set()
        for idx, seqid in enumerate(readnames):
            if seqid.find('r.n') != -1:
                need_2_correct.add(idx)
            else:
                good.add(idx)
        goodkmers = collections.Counter()
        for idx in good:
            goodkmers.update(gen_kmer(reads[idx],noRC=True, kmersize=curKmerSize))
        corrected = 0
        for idx in need_2_correct:
            positive = 0
            negative = 0
            c = gen_kmer(reads[idx], noRC=True, kmersize=curKmerSize)
            for kmer in c:
                if kmer in goodkmers:
                    positive += 1
            c = gen_kmer(reverse_comp(reads[idx]), noRC=True, kmersize=curKmerSize)
            for kmer in c:
                if kmer in goodkmers:
                    negative += 1
            if positive > negative:
                corrected += 1
                readnames[idx] = readnames[idx].rstrip('.n')
                print("correcting {0}, removing .n".format(readnames[idx]))
            elif positive < negative:
                corrected += 1
                print("correcting {0}, reverse, remove suffix".format(readnames[idx]))
                readnames[idx] = readnames[idx].split()[0]
                reads[idx] = reverse_comp(reads[idx])
            else:
                print("correcting {0} can not be corrected.".format(readnames[idx]))
        print("{0} out of {1} reads were corrected in this round".format(corrected, len(need_2_correct)))
        if corrected == 0:
            curKmerSize = curKmerSize - 1
            if curKmerSize < 18:
                break

    f = open(fname, 'w')
    for i in range(len(readnames)):
        f.write(">"+readnames[i] + "\n")
        f.write(reads[i] + "\n")
    f.close()


def combine_clusters(g, dict_read2cluster, dict_cluster2reads, combine_threshold):
    more_work = True
    while more_work:
        more_work = False
        nodes = set()
        edges = g.edges(data=True)
        edge_can_be_combined = []
        num_need_combine = 0
        for edge in edges:
            if edge[2]['weight'] >= combine_threshold:
                num_need_combine += 1
                if edge[0] in nodes or edge[1] in nodes:
                    continue
                else:
                    nodes.add(edge[0])
                    nodes.add(edge[1])
                    edge_can_be_combined.append(edge)
        logging.info("{0} out of {1} edges can be combined directly in this round. \n".format(len(edge_can_be_combined), num_need_combine))
        logging.info("Nodes in graph: {0}, edges in graph: {1}\n".format(g.number_of_nodes(), g.number_of_edges()))
        if num_need_combine > len(edge_can_be_combined):
            more_work = True
        for edge in edge_can_be_combined:
            if edge[1] == edge[0]:  # same nodes
                logging.info("ERROR: self loop, should never happen. Node: {0}\n".format(edge[0]))
                continue
            v1 = edge[0]
            v2 = edge[1]

            #if len(dict_cluster2reads[v1]) >= cluster_threshold:
                #logging.info("Cluster {0} too large ({1}), leave it along\n".format(v1, len(dict_cluster2reads[v1])))
            #    continue
            #if len(dict_cluster2reads[v2]) >= cluster_threshold:
                #logging.info("Cluster {0} too large ({1}), leave it along\n".format(v2, len(dict_cluster2reads[v2])))
            #    continue
            #logging.info("Combining {0} and {1}, weight {2}\n".format(edge[0], edge[1], edge[2]['weight']))

            cnt_pe_read = g.node[v1]['weight'] + g.node[v2]['weight']
            cnt_pe_read += edge[2]['weight']
            newNode = get_joint_name(v1, v2)
            g.add_node(newNode, weight = cnt_pe_read)
            for v1e in g.edges(v1, data=True):
                if v1e[1] != v2:
                    g.add_edge(newNode, v1e[1], weight=v1e[2]['weight'])
            g.remove_node(v1)

            for v2e in g.edges(v2, data=True):
                if v2e[1] != v1:
                    g.add_edge(newNode, v2e[1], weight=v2e[2]['weight'])
            g.remove_node(v2)


            for n in g.neighbors(newNode):
                edges_two_nodes = g.edges(newNode, data=True)
                w = 0
                for e in edges_two_nodes:
                    if e[1] == n:
                        w += e[2]['weight']
                        #logging.info("Edge {0} is combined....\n".format(e))
                        g.remove_edge(e[0], e[1])
                g.add_edge(n, newNode, weight = w)

            #print(g.edges(newNode, data=True))
            #logging.info("{0}\t{1}\n".format(v1, len(dict_cluster2reads[v1])))
            #logging.info("{0}\t{1}\n".format(v2, len(dict_cluster2reads[v2])))
            update_cluster(dict_read2cluster, dict_cluster2reads, v1, v2, newNode)
            #logging.info("{0}\t{1}\n".format(newNode, len(dict_cluster2reads[newNode])))



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


def has_large_inout_nodes(g, cutoff=5):
    for n in nx.average_neighbor_degree(g).values():
        if n > cutoff:
            return True
    return False


def make_dag(g):
    if nx.is_directed_acyclic_graph(g):
        return
    p = nx.periphery(g)
    for c in nx.weakly_connected_component_subgraphs(g):
        if nx.is_directed_acyclic_graph(g):
            continue
        cycles = nx.simple_cycles(c)
        for c in cycles:
            edges = zip(c[:-1], c[1:])
            edges.append((c[-1], c[0]))
            for e in edges:
                data = g.edges(e[0], e[1])[0][2]
                c.remove_edge(e[0], e[1])

def longestPathAll(g, dict_seq):
    paths = []
    if g.number_of_nodes() == 0:
        return paths
    if g.number_of_nodes() == 1:
        paths.append(g.nodes())
        return paths
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
    for p in paths:
        for node in p:
            g.remove_node(node)
    otherPaths = longestPathAll(g, dict_seq)
    if otherPaths:
        paths.extend(otherPaths)
    return paths


def longestPathAll1(g, dict_seq):
    print(g.number_of_nodes())
    paths = []
    if g.number_of_nodes() == 0:
        return paths
    if g.number_of_nodes() == 1:
        paths.append(g.nodes())
        g.remove_node(g.nodes()[0])
        return paths
    dist = {}
    for node in nx.topological_sort(g):
        pairs = [(dist[v][0]+len(dict_seq[node])- g[v][node]['weight'], v) for v in g.pred[node]]
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
    for p in paths:
        #print(p)
        for node in p:
            g.remove_node(node)
    otherPaths = longestPathAll1(g, dict_seq)
    if otherPaths:
        paths.extend(otherPaths)
    for p in otherPaths:
        for node in p:
            g.remove_node(node)
    return paths

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

def longestPathAll2(g, dict_seq):
    allpaths = []
    if is_linear_graph(g)[0]:
        p = get_path_linear_graph(g)
        return [p]
    while g.number_of_nodes() > 0:
        #logging.info(str(g.number_of_nodes()))
        #logging.info('\n')
        paths = longestPath(g, dict_seq)
        if paths:
            allpaths.extend(paths)
            for p in paths:
                for n in p:
                    g.remove_node(n)
    return allpaths

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


def simplify_graph(g):
    for e in g.selfloop_edges():
        g.remove_edge(e[0], e[1])
    for node in g.nodes():
        neighbors = list(nx.all_neighbors(g, node))
        edges = g.in_edges(node, data=True)
        edges.extend(g.out_edges(node, data=True))
        plus = []
        minus = []
        for e in edges:
            if e[2][node] == '+':
                plus.append(e)
            else:
                minus.append(e)
            if not plus or not minus:
                continue
            if len(plus) >= len(minus):
                for e in minus:
                    if g.has_edge(e[0], e[1]):
                        g.remove_edge(e[0], e[1])
            if len(plus) <= len(minus):
                for e in plus:
                    if g.has_edge(e[0], e[1]):
                        g.remove_edge(e[0], e[1])
    remove_out_tips(g)
    remove_in_tips(g)
    for c in nx.weakly_connected_component_subgraphs(g):
        if c.number_of_nodes() <= 2:
            continue
        isLinear, ends, source, sink= is_linear_graph(c)
        if isLinear:
            if sink==1 and source == 1:
                continue
            adjust_edge_di(g, c, ends[0], ends[1])



def has_good_repeat(seq, minRepeat=20, maxRepeat=55, minSpacer=20, maxSpacer=60):
    seqLen = len(seq)
    searchEnd = seqLen-minRepeat-minRepeat-minSpacer
    for i in range(searchEnd):
        beginSearch = i + minRepeat + minSpacer
        endSearch = i + maxSpacer + maxRepeat
        if endSearch > seqLen:
            endSearch = seqLen
        if endSearch < beginSearch:
            endSearch = beginSearch
        p = seq[i:i+10]
        index = seq.find(p, i+minRepeat+minSpacer)
        if index > 0:
            numExtend = extend(seq, i, index, 10, maxmismatch=1)
            if minSpacer <= index-(i+10+numExtend):
                if minRepeat+numExtend <= maxRepeat:
                    return seq[i:i+10+numExtend], seq[index:index+10+numExtend]
    return ""


def gen_kmer(seq, kmersize=10, noRC=False):
    counter = collections.Counter()
    for i in range(0, len(seq)-kmersize):
        counter.update({seq[i:kmersize+i]:1})
        if not noRC:
            counter.update({reverse_comp(seq[i:kmersize+i]):1})
    return counter

def extend(seq, pos1, pos2, winLen, maxmismatch=1):
    numExtend = 0
    numExtendAfter = 0
    mismatch = 0
    seqLen = len(seq)
    i = pos1+winLen
    j = pos2+winLen
    while j < seqLen:
        if seq[i] == seq[j]:
            i += 1
            j += 1
            if mismatch >=1:
                numExtendAfter += 1
            else:
                numExtend += 1
        else:
            if mismatch < maxmismatch:
                i += 1
                j += 1
                mismatch += 1
                continue
            if mismatch >= maxmismatch:
                break
    if numExtendAfter <= 1:
        return numExtend
    else:
        return numExtend + mismatch + numExtendAfter

def find_all_occ(seq, pattern):
    start = 0
    ret = []
    while True:
        index = seq.find(pattern, start)
        if index >=0:
            ret.append(index)
            start = index + len(pattern)
        else:
            break
    return ret

def remove_in_tips(g, keepNode=True):
    for n in g.nodes():
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
                        if keepNode:
                            g.remove_edge(e[0], e[1])
                        else:
                            g.remove_node(n)
                        break

def remove_out_tips(g, keepNode=True):
    for n in g.nodes():
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
                        if keepNode:
                            g.remove_edge(e[0], e[1])
                        else:
                            g.remove_node(n)
                        break

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
        g.add_edge(e, s,  {'weight':g[s][e]['weight'], 'label':g[s][e]['label'][::-1], e:g[s][e]['label'][1],s:g[s][e]['label'][0]})
        g.remove_edge(s, e)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Combining and filter read clusters.'
    )
    parser.add_argument('-d','--debug',
                        help='Print lots of debugging statements',
                        action="store_const",dest="loglevel",const=logging.DEBUG,
                        default=logging.INFO
                    )

    parser.add_argument('-o', '--outputFolder', help='Save overlap graphs to pdf files', default="./")
    parser.add_argument('-n', '--namePrefix', help='Prefix for naming output files.')
    parser.add_argument('-p', '--pairend',  help='Number of paired end reads support to combine clusters.', type=int, default=1)
    parser.add_argument('-i', '--onlyset',  help='Only check clusters in this file')
    parser.add_argument('-t', '--discard', help='Discard bad reads.', action='store_true')
    parser.add_argument('samfile', help='SAM file.')
    parser.add_argument('clusterfile', help='Cluster file')
    parser.add_argument('fastafile', help='Fasta file without duplicate/contained reads')
    parser.add_argument('graphfile',  help='Overlap graph file')
    parser.add_argument('fastafull',  help='Fasta file after filtration step')
    args = parser.parse_args()
    FORMAT = '%(asctime)s - %(module)-16s - %(levelname)s - %(message)s'
    logging.basicConfig(level=args.loglevel, format=FORMAT)
    rootLogger = logging.getLogger()
    fileHandler = logging.FileHandler("combine.log")
    rootLogger.addHandler(fileHandler)

    onlyset = set()
    if args.onlyset:
        with open(args.onlyset) as f:
            for line in f:
                onlyset.add(line.strip())
    readsBaseNames, dict_DRReads, dict_read_repeatinfo = get_read_basename_and_dr_reads(args.fastafull)
    logging.info("Reading clusters from sga result\n")
    dict_read2cluster = get_which_cluster_for_reads(args.clusterfile)
    logging.info("Get single clustered read...\n")
    set_single_read = get_single_read_cluster(args.fastafile, dict_read2cluster)
    logging.info("Number of single clustered reads: {0}\n".format(len(set_single_read)))
    logging.info("Classifying single clustered read using pe info\n")
    classify_single_reads(set_single_read, dict_read2cluster, dict_DRReads, False)
    single_read_in_sam = enrich_cluster_using_sam(args.samfile, dict_read2cluster)
    classify_single_reads(set_single_read, dict_read2cluster, dict_DRReads, False)
    sam_single_reads = single_read_in_sam.keys()
    classify_single_reads(sam_single_reads, dict_read2cluster, dict_DRReads, False)
    dict_cluster2reads = cluster_to_reads(dict_read2cluster)
    dict_all_reads = read_full_fasta(args.fastafull)
    dict_pe_count, dict_pe_count_in_cluster = get_pe_count_clusters(readsBaseNames, dict_read2cluster, dict_DRReads)
    gc.collect()

    logging.info("Number of clusters before combining: {0}\n".format(len(dict_cluster2reads)))

    cnt_no_pe_reads = 0
    for k in dict_cluster2reads:
        if k not in dict_pe_count_in_cluster:
            cnt_no_pe_reads += 1

    logging.info("Number of clusters without any pe reads: {0}\n".format(cnt_no_pe_reads))

    g = nx.MultiGraph()
    for c in dict_cluster2reads:
        num_pe_reads = 0
        if c in dict_pe_count_in_cluster:
            num_pe_reads = dict_pe_count_in_cluster[c]
        g.add_node(c, weight=num_pe_reads)

    for c in dict_pe_count:
        c1, c2 = c.split('@')
        g.add_edge(c1, c2, weight=dict_pe_count[c])

    logging.info("Number of reads in all clusters: {0}\n".format(len(dict_read2cluster)))
    logging.info("Number of nodes: {0}, number of edges: {1}\n".format(len(g.nodes()), len(g.edges())))

    logging.info("Combining clusters using paired end reads info\n")
    combine_clusters(g, dict_read2cluster, dict_cluster2reads, args.pairend)
    logging.info("Number of nodes: {0}, number of edges: {1}\n".format(len(g.nodes()), len(g.edges())))
    gc.collect()

    dict_dr_reads_in_cluster = get_drread_count_in_clusters(dict_cluster2reads, dict_DRReads)

    for e in g.edges(data=True):
        if e[2]['weight'] < combine_threshold:
            g.remove_edge(e[0], e[1])
    dict_read_2_overlap = read_connected_components(args.graphfile)
    filter_clusters(dict_cluster2reads, dict_read_2_overlap, dict_all_reads, dict_DRReads, dict_dr_reads_in_cluster, args.discard, onlyinset=onlyset)
