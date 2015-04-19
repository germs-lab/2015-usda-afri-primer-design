import sys

d = {}
for n, line in enumerate(open(sys.argv[1])):
    if n > 0:
        dat = line.rstrip().split(',')
        bp = dat[-6]
        run_id = dat[0]
        d[run_id] = int(bp)

d2 = {}
for n, line in enumerate(open(sys.argv[2])):
    if n == 0:
        dat = line.rstrip().split('\t')
        for n2, y in enumerate(dat[1:]):
            d2[n2] = y

d3 = {}
for key in d2:
    d3[key] = d[d2[key]]

for n, line in enumerate(open(sys.argv[2])):
    if n == 0:
        print line.rstrip()
    else:
        dat = line.rstrip().split('\t')
        for n3, i in enumerate(dat):
            if n3 == 0:
                continue
            else:
                #print n3, i, d3[n3-1], float(i)/d3[n3-1]
                dat[n3] = str(float(i)/d3[n3-1])
        print '\t'.join(dat)
