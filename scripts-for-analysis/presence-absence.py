import sys

for n, line in enumerate(open(sys.argv[1])):
    if n == 0:
        print line.rstrip()
    else:
        dat = line.rstrip().split('\t')
        for n2, x in enumerate(dat):
            if n2 > 0:
                if int(x) > 0:
                    dat[n2] = str(1)
                else:
                    dat[n2] = str(0)
        print '\t'.join(dat)
                    
    
