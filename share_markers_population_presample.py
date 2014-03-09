#!/usr/bin/env python

# Usage: <script> sharefile mapfile affectedlist fullinds outfile
import sys
from itertools import combinations, izip
from random import sample
from multiprocessing import Pool

import numpy as np 

nrep = 10000

parallel = True
nthreads = 8

# Unsigned 16 bit integer
datatype = 'u4'

def shares(inds):
    ninds = len(inds)
    tmaxshares = numpairs(ninds)
    s = np.zeros(nmark, dtype=datatype)
    validpairs = {frozenset(x) for x in combinations(inds, 2)} & keyset
    for pair in validpairs:
        for start, stop in shared[pair]:
            s[start:(stop+1)] += 1
    return s / tmaxshares

def empirical_p(observed, nullvalues):
    return nullvalues[nullvalues >= observed].shape[0] / float(nullvalues.shape[0])

def numpairs(n):
    return n * (n-1) * 0.5

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
        start,stop = [posd[int(x)] for x in l[5:7]]
        if pair not in shared:
            shared[pair] = []
        shared[pair].append([start, stop])
    keyset = frozenset(shared.keys())

if numpairs(naff) > np.iinfo(datatype).max:
    print 'More affected pairs than can fit into datatype %s' % datatype
    print 'Number of affected pairs: %s' % numpairs(naff)
    print 'Max possible with datatype %s: %s' % (datatype,
                                                 np.iinfo(datatype).max)
    exit(1)
print 'Calculating sharing from affecteds'
affshare = shares(affinds)


print 'Calculating sharing from %s draws of %s individuals' % (nrep, naff)

def nsharehelper(x):
    if x % 1000 == 0:
        print 'Random draw %s' % x 
    return shares(sample(fullinds, naff))

if parallel:
    pool = Pool(processes=nthreads)
    nullshares = pool.map(nsharehelper, xrange(nrep))
else:
    nullshares = map(nsharehelper, xrange(nrep))

nullshares = np.array(nullshares).T

print 'Calculating Empirical P-values'
pvals = np.array([empirical_p(observed, nulls)
                  for observed,nulls in izip(affshare, nullshares)])

print 'Minimum observed P: %s' % min(pvals)
print
print 'Writing output'
with open(sys.argv[5],'w') as f:
    f.write(','.join(['chr','snp','cm','pos','pctshares', 'p%s' % nrep]) + '\n')
    for m,a,p in izip(gmap, affshare, pvals):
        chr,snp,cm,pos = m
        f.write(','.join(str(x) for x in [chr,snp,cm,pos,a,p]))
        f.write('\n')
