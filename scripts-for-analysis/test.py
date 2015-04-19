import sys

for line in open(sys.argv[1]):
    dat = line.rstrip().split('\t')
    print len(dat)
