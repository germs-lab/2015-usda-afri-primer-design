import sys

for n, line in enumerate(open(sys.argv[1])):
    print ">" + sys.argv[2].split('_')[0] + "v" + str(n) + "_" + sys.argv[2].split('_')[1]
    print line.rstrip()
