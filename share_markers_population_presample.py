#!/usr/bin/env python

# Usage: <script> sharefile mapfile affectedlist fullinds outfile
import sys
from itertools import combinations, izip
from random import sample
from multiprocessing import Pool

import numpy as np 

nrep = 1000

parallel = True
nthreads = 8


def share(pair):
    retval = np.zeros(nmark)
    shares = shared[pair]

    for start, stop in shares:
        a,b = posd[start],posd[stop]
        retval[a:(b+1)] = 1 

    return retval


def shares(inds):
    ninds = len(inds)
    tmaxshares = 0.5 * ninds * (ninds - 1)

    nshare = np.zeros(nmark)
    validpairs = {frozenset(x) for x in combinations(inds, 2)} & keyset
    for pair in validpairs:
        pair = frozenset(pair)
        nshare += share(pair)
    return nshare / tmaxshares

print 'Reading individual lists'
with open(sys.argv[3]) as f:
    affinds = set('.'.join(x.strip().split()) for x in f)
    naff = len(affinds)
with open(sys.argv[4]) as f:
    fullinds = set('.'.join(x.strip().split()) for x in f)

print 'Reading map'
with open(sys.argv[2]) as f:
    gmap = [x.strip().split() for x in f]
    positions = [int(x[3]) for x in gmap]
    nmark = len(positions)
    posd =  dict([(y,x) for x,y in enumerate(positions)])

print 'Reading Share file'
with open(sys.argv[1]) as sharef:
    shared = {}
    for i,line in enumerate(sharef):
        if i % 100000 == 0:
            print 'reading line %s' % i
        l = line.strip().split()
        ind1 = '.'.join(l[0:2])
        ind2 = '.'.join(l[2:4])
        pair = frozenset([ind1,ind2])
        start,stop = [int(x) for x in l[5:7]]
        if pair not in shared:
            shared[pair] = []
        shared[pair].append([start, stop])
    keyset = frozenset(shared.keys())

print 'Calculating sharing from affecteds'
affshare = shares(affinds)


print 'Calculating sharing from %s draws of %s individuals' % (nrep, naff)
def nsharehelper(x):
    print 'sharecall %s' % x 
    return shares(sample(fullinds, naff))

if parallel:
    pool = Pool(processes=nthreads)
    nullshares = pool.map(nsharehelper, xrange(nrep))
else:
    nullshares = map(nsharehelper, xrange(nrep))

nullshares = zip(*nullshares)



print 'Calculating Empirical P-values'
pvals = [sum(1 for ns in n if ns >= a) / float(nrep)
         for a,n in izip(affshare, nullshares)]

print 'Minimum observed P: %s' % min(pvals)
print
print 'Writing output'
with open(sys.argv[5],'w') as f:
    f.write(','.join(['chr','snp','cm','pos','pctshares', 'p%s' % nrep]) + '\n')
    for m,a,p in izip(gmap, affshare, pvals):
        chr,snp,cm,pos = m
        f.write(','.join(str(x) for x in [chr,snp,cm,pos,a,p]))
        f.write('\n')
