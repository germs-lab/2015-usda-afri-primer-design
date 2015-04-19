import sys

for line in open(sys.argv[1]):
    dat = line.rstrip().split('\t')
    for n, x in enumerate(dat):
        if n == 2:
            taxon = dat[n].split(';')
            for n2, y in enumerate(taxon):
                dat.append(y.strip())
    print '\t'.join(dat)
