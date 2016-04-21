# metaCRISPR

metaCRISPR is a tool to assemble CRISPRs (Clustered Regularly Interspaced Short
Palindromic Repeats and Associated Proteins) from metagenomic sequencing data without relying on generic assembly, which is
error-prone and computationally expensive for complex data. It
can run on commonly available machines in small labs. It employs
properties of CRISPRs to decompose generic assembly into local
assembly.

# 1. How to install

## 1.1 Dependencies.
To run MetaCRISPR, one needs to have the following tools/packages installed:

1. Java.
2. [Maven](https://maven.apache.org/). This is needed to make the java package for read identification.
3. Python
4. [networkx](https://networkx.github.io/)
5. [Genometools](https://github.com/genometools/genometools). MetaCRISPR uses readjoiner to construct overlap graph from reads.
6. [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml).

For Linux users these tools/packages can be easily installed using package management systems of the Linux distribution (e.g. apt-get on Ubuntu). For MacOS users, they can be installed using [HomeBrew](http://brew.sh/).

## 1.2 Make Java package.
In the `ReadRecuiter` folder, run:

    mvn package

This will generate a jar file in the ReadRecuiter/target folder. The jar can be used to identify CRISPR reads from metagenomic reads data.


# 2. Assemble CRISPRs.
MetaCRIPSRs pipeline has three major steps: identify CRISPR reads, cluster CRISPR reads, and assembly.

## 2.1 CRISPR reads identification

    java -cp ReadRecuiter-0.0.1-SNAPSHOT-jar-with-dependencies.jar crispr.ReadRecuiter  -input test.fa -minRepeat 23 -maxRepeat 60 -minSpacer 20 -
maxSpacer 80 -threads 4 -prefix  test -maxMismatch 1  -repeats repeats.db..txt

This command identify CRISPR reads from input read file test.fa. Command line options are:
`minRepeat`: minimum repeat size
`maxRepeat`: maximum repeat size
`minSpacer`: minimum spacer size
`maxSpacer`: maximum spacer size
`threads`: number of threads to use to do read identification.
`maxMismatch`: maximum mismatch allowed in repeat sequences.
`prefix`: prefix for output files.
`repeats`: a file that contains known CRISPR direct repeats. Each line contains the sequence of a known repeat. This is optional.


If one uses N threads to run the identification program, it will generate N sets of files whose names all starts with the prefix set by user. To run the next step, concatenate N files with name prefix.{0-N-1}.filtered.fa together to form one file:

    cat prefix.*.filtered.fa > filtered.reads.fa
~/CRISPRFinder-rsycn
This file contains all filtered reads.

## 2.2 CRISPR reads clustering

    python scripts/run-rj-on-filtered-reads.py -b -t1 4 -t2 4 -o rj filtered.reads.fa mock  80

This will generate all clustered result data to folder `rj` (-o option). `-t1` option specifies the number of threads to use to run readjoiner. `-t2` specifies the number of threads to use to run bowtie.
There are three required parameters:
1. Filtered reads fasta file (`filterd.reads.fa` in the example)
2. Name of the read set for Readjoiner. One can think of this as an identifier to identify the sample of the data. (`mock` in the example).
3. Overlap size to build the overlap graph.

After generating the clusters, using the following command to filter the clusters:

    mkdir clusters
    cd clusters
    python scripts/combine-and-filter-clusters.py ../rj/bowtie/bowtie.result.sam ../rj/rj/mock.clusters.0.txt ../rj/rj/mock.reads.filtered.fa ../rj/rj/mock.graph.rename.txt .../filtered.reads.fa

This produces reads for each belongs to each clusters to the clusters folder.

## 2.3 Assembly

    python scripts/run-rj-each-cluster.py   clusters   scripts/run-rj-on-filtered-reads.py cluster-each-rj
    python scripts/crisprfinder-rj.py  clusters rj/mock.readgroups.txt cluster-each-rj/rj 30 50

Here `clusters` is the folder that contains all the fasta sequences of each cluster generated in the previous step. `30` is the size of overlap to generate the overlap graph in the initial assemble step. `50` is the confident edge size.
