import sys
import os
import argparse
import logging
import string

def run_rj(fname, readsetname, graphprefix, overlap=80, nThread=1):
    command1 = 'gt readjoiner prefilter -db  ' + fname + ' -readset ' + readsetname + ' -q -des'
    print(command1)
    ret = os.system(command1)
    if ret:
        sys.stderr.write("Error occurred when calling gt readjoiner prefilter\n")
        sys.exit(-1)
    command2 = 'gt -j ' + str(nThread) + ' readjoiner overlap -readset ' + readsetname + ' -l ' + str(overlap) + ' -memlimit 1GB'
    print(command2)
    ret = os.system(command2)
    if ret:
        sys.stderr.write("Error occurred when calling gt readjoiner overlap\n")
        sys.exit(-1)
    dirname = os.path.dirname(readsetname)
    command3 = 'gt readjoiner spmtest -readset ' + readsetname + ".0 -test showlist |awk  '$1 != $3' > " + dirname + '/' + graphprefix + '.graph.txt'
    print(command3)
    ret = os.system(command3)
    if ret:
        sys.stderr.write("Error occurred when calling gt readjoiner spmtest\n")
        sys.exit(-1)


def valid_name(name):
    if not name:
        return False
    for c in name:
        if c in string.printable:
            continue
        else:
            return False
    return True

def extract_info_from_gt_result(graphfile, desfile, fastafullfile, fasta_nodup, rjfolder, prefix):

    cnt = 0
    dict_name_mappings = {}
    set_nodup_readnames = set()
    # mapping between readjoiner read names and true read names
    with open(desfile) as f:
        for line in f:
            if not valid_name(line.strip()):
                break
            dict_name_mappings[cnt] = line.strip()
            set_nodup_readnames.add(line.strip())
            cnt += 1
    # save non duplicate reads to a file
    fout_nodup = open(fasta_nodup, 'w')
    with open(fastafullfile) as f:
        seq = ''
        line = f.readline()
        sp = line.strip().split()
        readid = sp[0].lstrip(">")
        for line in f:
            if line.startswith('>'):
                if readid in set_nodup_readnames:
                    fout_nodup.write('>' + readid + '\n')
                    fout_nodup.write(seq)
                sp = line.strip().split()
                readid = sp[0].lstrip(">")
            else:
                seq = line
        if readid in set_nodup_readnames:
            fout_nodup.write('>' + readid + '\n')
            fout_nodup.write(seq)

    num_nodup = len(set_nodup_readnames)
    print(num_nodup)
    generate_conected_components(graphfile, rjfolder, prefix, dict_name_mappings)

def generate_conected_components(graphfile, rjfolder, prefix, dict_name_mapping):

    num_nodes = len(dict_name_mapping)
    rank = [0] * num_nodes
    parent = range(num_nodes)

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        x = find(x)
        y = find(y)
        if x == y:
            return
        if rank[x] < rank[y]:
            parent[x] = y
        elif rank[x] > rank[y]:
            parent[y] = x
        else:
            parent[y] = x
            rank[x] += 1

    with open(graphfile) as f:
        for line in f:
            sp = line.split()
            union(int(sp[0]), int(sp[2]))
    dict_componet2reads = {}
    for n in range(num_nodes):
        c = find(n)
        if c  not in dict_componet2reads:
            dict_componet2reads[c] = set()
        dict_componet2reads[c].add(dict_name_mapping[n])

    clustername = os.path.join(rjfolder, prefix+".clusters.0.txt")
    fcluster = open(clustername, 'w')
    for c in dict_componet2reads:
        fcluster.write(str(c))
        for r in dict_componet2reads[c]:
            fcluster.write('\t' + r)
        fcluster.write('\n')
    fcluster.close()

    renamedGraph = os.path.join(rjfolder, prefix+".graph.rename.txt")
    fgraph = open(renamedGraph, 'w')
    with open(graphfile) as f:
        for line in f:
            sp = line.split()
            fgraph.write(dict_name_mapping[int(sp[0])] + '\t')
            fgraph.write(sp[1] + '\t')
            fgraph.write(dict_name_mapping[int(sp[2])] + '\t')
            fgraph.write(sp[3] + '\t')
            fgraph.write(sp[4] + '\n')
    fgraph.close()


def run_bowtie(fullfasta, nodupfasta, bowtiefolder, readgroupname, nThread=1):
    command = 'bowtie-build -f ' + nodupfasta + ' ' + bowtiefolder + '/index > ' + bowtiefolder + '/build.log 2> '+ bowtiefolder + '/build.err'
    print(command)
    ret = os.system(command)
    ret = 0
    if ret:
        sys.stderr.write("Error occurred when running bowtie-build, exit.\n")
        sys.exit(1)
    command = 'bowtie -v 0 -a -f -p ' + str(nThread) + ' -S ' + bowtiefolder + '/index ' + fullfasta
    command += ' |grep -v "@" |awk'
    command += " '$2!=4' > " + bowtiefolder + "/bowtie.result.sam"
    print(command)
    ret = os.system(command)
    if ret:
        sys.stderr.write("Error occurred when running bowtie, exit.\n")
        sys.exit(1)
    extract_contained_readgroup(bowtiefolder+"/bowtie.result.sam", readgroupname)

def extract_contained_readgroup(samfile, outputname):
    dict_readgroups = {}
    with open(samfile) as f:
        for line in f:
            if not line.startswith("@"):
                sp = line.split()
                if sp[-1].split(":")[-1] != '0':
                    continue
                if sp[2] not in dict_readgroups:
                    dict_readgroups[sp[2]] = set()
                dict_readgroups[sp[2]].add(sp[0])
    fout = open(outputname, 'w')
    for r in dict_readgroups:
        fout.write(r)
        for r1 in dict_readgroups[r]:
            if r1 != r:
                fout.write('\t' + r1)
        fout.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Script to run readjoiner on filtered dataset.'
    )
    parser.add_argument('-d','--debug',
                        help='Print lots of debugging statements',
                        action="store_const",dest="loglevel",const=logging.DEBUG,
                        default=logging.INFO
                    )
    parser.add_argument('-o', '--out', help='Outputfolder. If not specified, using the current folder')
    parser.add_argument('-b', '--bowtie', help='Run bowtie.', action="store_true")
    parser.add_argument('-t1', '--threadRJ', help='Number of threads to run Readjoiner', type=int, default=1)
    parser.add_argument('-t2', '--threadB', help='Number of threads to run Bowtie', type=int, default=1)

    #parser.add_argument('-p', '--prefix', help='Prefix for output names')
    parser.add_argument('fasta', help='Filtered fasta file.')
    parser.add_argument('readset', help='Name of the read set for Readjoiner')
    parser.add_argument('overlap', help='Overlap threashold for building the overlap graph.')
    args = parser.parse_args()
    #FORMAT = '%(asctime)-15s %(message)s'
    logger = logging.getLogger(__name__)
    FORMAT = '%(asctime)s - %(module)-16s - %(levelname)s - %(message)s'
    logging.basicConfig(level=args.loglevel, format=FORMAT)

    outfolder = os.getcwd()
    rjfolder = ''
    bowtiefolder = ''
    if args.out:
        outfolder = args.out
        if not os.path.exists(outfolder):
            logging.info("Creating output folder")
            os.mkdir(outfolder)
        logging.info("Working directory: {0}".format(outfolder))

    rjfolder = os.path.join(outfolder, 'rj/')
    if not os.path.exists(rjfolder):
        logging.info("Creating rj folder")
        os.mkdir(rjfolder)
    if args.bowtie:
        bowtiefolder = os.path.join(outfolder, 'bowtie')
        if not os.path.exists(bowtiefolder):
            logging.info("Creating bowtie folder")
            os.mkdir(bowtiefolder)
    logging.info("Output folder: {0}".format(outfolder))
    logging.info("RJ folder: {0}".format(rjfolder))
    logging.info("Readset name: {0}".format(args.readset))
    logging.info("Thread RJ: {0}".format(args.threadRJ))
    logging.info("Thread Bowtie: {0}".format(args.threadB))
    run_rj(args.fasta, os.path.join(rjfolder, args.readset), args.readset, overlap=args.overlap, nThread=args.threadRJ)
    graphfile = os.path.join(rjfolder, args.readset+".graph.txt")
    rjdesfile = os.path.join(rjfolder, args.readset+".des")
    filteredFasta = os.path.join(rjfolder, args.readset+".reads.filtered.fa")
    logging.info("Graph file name: {0}".format(graphfile))
    logging.info("RJ des file name: {0}".format(rjdesfile))
    logging.info("Filtered fasta file name: {0}".format(filteredFasta))
    extract_info_from_gt_result(graphfile, rjdesfile, args.fasta, filteredFasta, rjfolder, args.readset)
    if args.bowtie:
        readgroupname = os.path.join(outfolder, args.readset+'.readgroups.txt')
        run_bowtie(args.fasta, filteredFasta, bowtiefolder, readgroupname, nThread=args.threadB)
