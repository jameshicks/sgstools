#!/usr/bin/env python

# Usage: <script> sharefile mapfile affectedlist fullinds outfile
import sys
from itertools import combinations, izip
from random import sample
from multiprocessing import Pool

import numpy as np 

nrep=1000
nthreads=10


def share(pair):
    retval = np.zeros(nmark)
    shares = shared[pair]
    lshares = len(shares)
    sidx = 0
    start, stop = shares[sidx]
    for i,p in enumerate(positions):
        if p > stop:
            sidx += 1
            if (sidx + 1) > lshares:
                break
            start, stop = shares[sidx]
        if start <= p <= stop:
            retval[i] = 1
    return retval
    
def shares(inds):
    print ';call to shares'
    ninds = len(inds)
    tmaxshares = 0.5 * ninds * (ninds - 1)

    nshare = np.zeros(nmark)
    validpairs = {frozenset(x) for x in combinations(inds, 2)} & keyset
    for pair in validpairs:
        pair = frozenset(pair)
        nshare += share(pair)
    return nshare / tmaxshares


with open(sys.argv[3]) as f:
    affinds = set('.'.join(x.strip().split()) for x in f)
    naff = len(affinds)
with open(sys.argv[4]) as f:
    fullinds = set('.'.join(x.strip().split()) for x in f)
with open(sys.argv[2]) as f:
    map = [x.strip().split() for x in f]
    positions = [int(x[3]) for x in map]
    nmark = len(positions)
print ';reading shares'
with open(sys.argv[1]) as sharef:
    shared = {}
    for i,line in enumerate(sharef):
        if i % 100000 == 0:
            print ';reading line %s' % i
        l = line.strip().split()
        ind1 = '.'.join(l[0:2])
        ind2 = '.'.join(l[2:4])
        pair = frozenset([ind1,ind2])
        #if not pair <= affinds: continue
        start,stop = [int(x) for x in l[5:7]]
        if pair not in shared:
            shared[pair] = []
        shared[pair].append([start, stop])
    print ';sorting'
    for k in shared:
        shared[k] = sorted(shared[k])
    keyset = frozenset(shared.keys())

print ';calc affectedshares'
affshare = shares(affinds)

print ';calc nullshares'

def nsharehelper(x):
    print ';sharecall %s' % x 
    return shares(sample(fullinds, naff))

pool = Pool(processes=nthreads)

nullshares = pool.map(nsharehelper, xrange(nrep))
nullshares = zip(*nullshares)



print ';calc pvals'
pvals = [sum(1 for ns in n if ns >= a) / float(nrep)
         for a,n in izip(affshare, nullshares)]


with open(sys.argv[5],'w') as f:
    
    for m,a,p in izip(map, affshare, pvals):
        chr,snp,cm,pos = m
        f.write(','.join(str(x) for x in [chr,snp,cm,pos,a,p]))
        f.write('\n')
