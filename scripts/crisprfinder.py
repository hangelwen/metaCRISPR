import graph
import sys
import logging
import argparse
import os
import cPickle
import gzip
import networkx as nx
from string import maketrans


def parse_readID(readID):
    '''
    Parse a readID line, returns the read name, depth, of the
    pairend reads, and other possible information (in a dict)

    readID format:
    ID/ID/ID/ID@ other_info
    '''
    dict_info = {}
    sp0 = readID.strip().split("@")
    sp1 = sp0[0].split("/")
    readname = sp1[0]
    list_readIDs = sp1[1:]
    depth = len(list_readIDs)
    if len(sp0) == 2:
        dict_info['other'] = sp0[1]
    else:
        dict_info['other'] = ""
    return readname, depth, list_readIDs, dict_info

def readInitalOverlapFromFile(fname):
    with gzip.open(fname) as f:
        line = f.readline()
        overlapthresh = int(line.split()[3].split(":")[-1])
        return overlapthresh


def buildGraphFromASQGFile(fnamegzip):
    with gzip.open(fnamegzip) as f:
        line = f.readline()
        overlapthresh = int(line.split()[3].split(":")[-1])
        olGraph = graph.OverlapGraph(overlapthresh)
        for line in f:
            sp = line.split()
            if sp[0] == 'VT':  # node
                readIDStr = sp[1]
                readname, depth, list_readIDs, dict_info = parse_readID(readIDStr)
                readSeq = sp[2]
                node = graph.NodeInfo(readSeq, readname, depth, dict_info)
                node.setPairedNodes(list_readIDs)
                olGraph.addNode(node)
            elif sp[0] == 'ED':  # edge
                readname1, depth1, list_readIDs1, dict_info1 = parse_readID(sp[1])
                readname2, depth2, list_readIDs2, dict_info2 = parse_readID(sp[2])
                orientation = "FF"
                if sp[9] == '1':
                    orientation = "FR"
                    continue
                overlapLen = int(sp[4]) - int(sp[3]) + 1
                dict_other_info = {"orientation": orientation}
                if sp[3] == '0':
                    olGraph.addEdge(readname2, readname1, overlapLen, dict_other_info)
                    if orientation == 'FR':
                        olGraph.FRNodes[readname1] = 1
                else:
                    olGraph.addEdge(readname1, readname2, overlapLen, dict_other_info)
                    if orientation == 'FR':
                        olGraph.FRNodes[readname2] = 1
        return olGraph



def testPath(p1, p2):
    d = set(p1).intersection(set(p2))
    if d:
        print(p1, p2)
        print(d)
        print("Done")

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
    parser.add_argument('asqf', metavar='ASQG-FILE', help='asqf file')
    parser.add_argument('peread', metavar='PairEnd-FILE', help='Fasta file that contain pair end reads.')
    parser.add_argument('overlap', metavar='confident-threshold', help='Confident edge overlap threshold.', type=int)
    args = parser.parse_args()
    #FORMAT = '%(asctime)-15s %(message)s'
    logger = logging.getLogger(__name__)
    FORMAT = '%(asctime)s - %(module)-16s - %(levelname)s - %(message)s'
    logging.basicConfig(level=args.loglevel, format=FORMAT)

    sizeLimit = 1000
    if args.sizeLimit:
        sizeLimit = args.sizeLimit

    overlapthresh = readInitalOverlapFromFile(args.asqf)
    confidentThresh = args.overlap
    logging.info("Initial overlap threshold in the asqg file: {0}".format(overlapthresh))
    logging.info("Confident edge overlap threshold: {0}".format(confidentThresh))
    if confidentThresh <= overlapthresh:
        logging.error("Confident threhold should be larger than the initial overlap threshold")
        sys.exit(1)
    if not os.path.exists(args.outputFolder):
        os.mkdir(args.outputFolder)


    olGraph = buildGraphFromASQGFile(args.asqf)
    logging.info("Size of the initial overlap graph:")
    logging.info("Nodes: {0}".format(olGraph.graph.number_of_nodes()))
    logging.info("Edges: {0}".format(olGraph.graph.number_of_edges()))
    logging.info("Component: {0}".format(olGraph.getNumComponent()))
    logging.info("=======================================================\n\n")

    logging.info("Loading pair end reads information")
    olGraph.addPairEndInfoFromFile(args.peread)


    if not nx.is_directed_acyclic_graph(olGraph.graph):
        logging.info("not DAG")
        #paths = nx.all_pairs_shortest_path(olGraph.graph)
        #for k in paths:
        #    print(k, len(paths[k]))
    if args.drawOverlap:
        logging.info("Drawing overlap graphes before cycle breaking\n\n")
        drawFolderBefore = os.path.join(args.outputFolder + "/og-before-cycle")
        if not os.path.exists(drawFolderBefore):
            os.mkdir(drawFolderBefore)
        olGraph.printCurrentGraph(drawFolderBefore, sizeLimit)


    dict_removed_edges = olGraph.cleanGraphUsingConfidentEdge(confidentThresh)
    logging.info("Size of the overlap graph after cleaning:")
    logging.info("Nodes: {0}".format(olGraph.graph.number_of_nodes()))
    logging.info("Edges: {0}".format(olGraph.graph.number_of_edges()))
    logging.info("Component: {0}".format(olGraph.getNumComponent()))
    logging.info("=======================================================\n\n")


    logging.info("Generating stringent overlap graph")
    dict_removed_edges_stringent = olGraph.genStringentGraph(confidentThresh)
    logging.info("Size of the stringent overlap graph:")
    logging.info("Nodes: {0}".format(olGraph.graph.number_of_nodes()))
    logging.info("Edges: {0}".format(olGraph.graph.number_of_edges()))
    logging.info("Component: {0}".format(olGraph.getNumComponent()))
    logging.info("=======================================================\n\n")


    #logging.info("Cleaning graph using overlap ratio...")
    #olGraph.cleanGraphEvenMore()
    #logging.info("Size of the overlap graph:")
    #logging.info("Nodes: {0}".format(olGraph.graph.number_of_nodes()))
    #logging.info("Edges: {0}".format(olGraph.graph.number_of_edges()))
    #logging.info("Component: {0}".format(olGraph.getNumComponent()))


    logging.info("=======================================================\n\n")
    logging.info("Removing tips...")
    olGraph.removeTips()
    logging.info("Size of the overlap graph:")
    logging.info("Nodes: {0}".format(olGraph.graph.number_of_nodes()))
    logging.info("Edges: {0}".format(olGraph.graph.number_of_edges()))





    logging.info("Checking cycles...")
    cycles = list(nx.simple_cycles(olGraph.graph))
    if cycles:
        logging.info("=======================================================\n\n")
        logging.info("Breaking cycles....")
        logging.info("Breaking cycles using paired end reads info....")
        olGraph.breakCyclesUsingPairedEndInfo()
        logging.info("Breaking cycles using simple heuristic....")
        olGraph.breakCyclesHeuristic()
        logging.info("=======================================================\n\n")
        logging.info("Breaking remaining cycles....")
        olGraph.breakOtherCycles()
        logging.info("=======================================================\n\n")

    # if cycles:
    #     logging.info("Breaking cycles using paired-end reads information.")
    #     for c in cycles:
    #         #logging.info(c)
    #         olGraph.breakCycles(c)
    if args.drawOverlap:
        logging.info("Drawing overlap graphes before cycle breaking\n\n")
        drawFolderAfter = os.path.join(args.outputFolder + "/og-after-cycle")
        if not os.path.exists(drawFolderAfter):
            os.mkdir(drawFolderAfter)
        if not os.path.exists(drawFolderAfter):
            os.mkdir(drawFolderAfter)
        olGraph.printCurrentGraph(drawFolderAfter, sizeLimit)

    cPickle.dump(olGraph, open(os.path.join(args.outputFolder, "olgraph.dump"), 'w'))

    logging.info("=======================================================\n\n")
    logging.info("Generating longest paths...")
    paths = olGraph.findAllLongestPathsFromComponents1()
    #for path in paths:
    #    olGraph.printPathWithInfo(path)
    clusterGraph, dict_path_info = olGraph.clusterPathes(paths)

    logging.info("=======================================================\n\n")
    logging.info("Generating contigs....")
    contigName = "contigs.fa"

    if args.namePrefix:
        contigName = args.namePrefix + "_" + contigName
    contigName = os.path.join(args.outputFolder, contigName)

    seqPrefix = "seq"
    if args.namePrefix:
        seqPrefix = args.namePrefix

    debugPathFile = None
    if args.savePath:
        debugPathFile = os.path.join(args.outputFolder, seqPrefix + ".pathinfo.txt")
    olGraph.genFastaFromClusters(clusterGraph, dict_path_info, contigName, seqPrefix, debugPathFile)
    #for g in nx.connected_component_subgraphs(clusterGraph):
    #    print("=================CLUSTER=========================")
    #    for node in g.nodes():
    #        olGraph.printPathWithInfo(dict_path_info[node]['path'])
    # for g in nx.weakly_connected_component_subgraphs(olGraph.graph):
    #     print(g.number_of_nodes(), g.number_of_edges())
    #     path =  longest_path(g)
    #     olGraph.printPathWithInfo(path)
    #     print(olGraph.genSeqFromPath(path))
