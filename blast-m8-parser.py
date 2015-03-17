import sys

d = {}

for line in open(sys.argv[1]):
    data = line.rstrip().split('\t')
    query = data[0]
    hit = data[1]
    identity = float(data[2])
    length = int(data[3])
    evalue = float(data[-2])

    if d.has_key(query):
        d[query].append(hit)
    else:
        d[query] = [hit]

for x in d.keys():
    #print sorted(d[x])
    d2  = {}
    for y in d[x]:
        y = y.split('_')[0]
        d2[y] = d2.get(y,0) + 1
    for x2 in d2:
        if x2 != "1":
            if d2[x2] > 1:
                forward = x2 + "_f"
                reverse = x2 + "_r"
                if forward in d[x] and reverse in d[x]:
                    print x, forward, reverse
