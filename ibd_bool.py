#!/usr/bin/env python

import sys
from itertools import combinations

matchfile, indfile, position, outfile = sys.argv[1:]
position = int(position)

delim = '\t'
minsegmentlength = (.5 * 10**6)
markerdensitylimit = (100 / float(10**6))

with open(indfile) as indf:
    inds = {'.'.join(x.strip().split()) for x in indf}

with open(matchfile) as sharef:
    shared = {}
    for i,line in enumerate(sharef):
        if i % 10000 == 0:
            print 'reading line %s' % i
        l = line.strip().split()
        ind1 = '.'.join(l[0:2])
        ind2 = '.'.join(l[2:4])
        markersinshare = int(l[9])
        pair = frozenset([ind1,ind2])
        start,stop = [int(x) for x in l[5:7]]
        
        if (stop - start) < minsegmentlength:
            continue
        
        if (markersinshare / float(stop - start)) < markerdensitylimit:
            continue
        
        if pair not in shared:
            shared[pair] = []
        shared[pair].append((start, stop))

def does_share(a,b, position):
    try:
        for start, stop in shared[frozenset([a,b])]:
            if start <= position <= stop:
                return True
        else:
            return False
    except KeyError:
        return False


with open(outfile,'w') as of:
    of.write(delim.join(['IND1', 'IND2', 'IBD_STATUS']))
    of.write('\n')
    for a,b in combinations(inds,2):
        of.write(delim.join([a, b, 'IBD' if does_share(a,b, position) else 'NOT_IBD']))
        of.write('\n')
