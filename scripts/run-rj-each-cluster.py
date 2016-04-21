import sys
import os

scriptfile = sys.argv[2]
pefiles = []
outputfolder = sys.argv[3]
for root, dirs, names in os.walk(sys.argv[1]):
    for name in names:
        if name.endswith(".pe.fa"):
            pefiles.append(os.path.join(root, name))
    break

for name in pefiles:
    base = os.path.basename(name).split(".")[0]
    #command = 'python ~/Dropbox/research/CRISPR/CRISPRFinder/src/run-rj-on-filtered-reads.py ' + name + ' ' + base + ' 40  -t1 3  -o cluster'
    command = 'python ' + scriptfile +  ' ' + name + ' ' + base + ' 30  -t1 2  -o ' + sys.argv[3]
    print(command)
    os.system(command)
