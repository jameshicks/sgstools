#!/usr/bin/env python

# Usage: <script> sharefile mapfile affectedlist fullinds outfile
import sys
from itertools import combinations, izip, imap
from random import sample
from multiprocessing import Pool

import numpy as np 

nrep = 100

parallel = False
nthreads = 8

markerdensitylimit = 100 / float(10**6)
minsegmentlength = .5 * 10**6

def shares(inds):
    ninds = len(inds)
    tmaxshares = numpairs(ninds)
    s = np.zeros(nmark, dtype=datatype)
    for pair in combinations(inds,2):
        try:
            for start, stop in shared[frozenset(pair)]:
                s[start:(stop+1)] += 1
        except KeyError:
            pass
    return s / tmaxshares

def numpairs(n):
    return n * (n-1) * 0.5

print 'Minimum segment length: %sMb' % (minsegmentlength / float(10**6))
print 'Minimum marker density: %s markers/Mb' % (markerdensitylimit *  float(10**6))

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

if naff == 0:
    print 'ERROR: No affected individuals!'
    exit(1)

if affinds - fullinds:
    print "Some individuals in the affected list weren't in the population list"
    print "These individuals were removed"
    affinds = affinds - fullinds

print 'Reading Share file'
with open(sys.argv[1]) as sharef:
    shared = {}
    for i,line in enumerate(sharef):
        if i % 100000 == 0:
            print 'reading line %s' % i
        l = line.strip().split()
        ind1 = '.'.join(l[0:2])
        ind2 = '.'.join(l[2:4])
        markersinshare = int(l[9])
        pair = frozenset([ind1,ind2])
        start,stop = [int(x) for x in l[5:7]]
        istart,istop = [posd[int(x)] for x in l[5:7]]

        if (stop - start) < minsegmentlength:
            continue

        if (markersinshare / float(stop - start)) < markerdensitylimit:
            continue
        
        if pair not in shared:
            shared[pair] = []
        shared[pair].append([istart, istop])
    keyset = frozenset(shared.keys())

for dtype in ['u1','u2','u4','u8']:
    datatype = dtype
    if numpairs(naff) < np.iinfo(datatype).max:
        break
else:
    print 'More affected pairs than can fit into datatype %s' % datatype
    print 'Number of affected pairs: %s' % numpairs(naff)
    print 'Max possible with datatype %s: %s' % (datatype,
                                                 np.iinfo(datatype).max)
    exit(1)

print 'Calculating sharing from affecteds'
affshare = shares(affinds)


print 'Calculating emperical p-values from %s draws of %s individuals' % (nrep, naff)

def nsharehelper(x):
    if x % 1000 == 0:
        print 'Random draw %s' % x 
    return shares(sample(fullinds, naff))

if parallel:
    pool = Pool(processes=nthreads)
    nullshares = pool.imap_unordered(nsharehelper, xrange(nrep), chunksize=250)
else:
    nullshares = imap(nsharehelper, xrange(nrep))

pvals = sum(n >= affshare for n in nullshares) / float(nrep)

print 'Minimum observed P: %s' % min(pvals)
print
print 'Writing output'

with open(sys.argv[5],'w') as f:
    f.write(','.join(['chr','snp','cm','pos','pctshares', 'p']) + '\n')
    for m,a,p in izip(gmap, affshare, pvals):
        chr,snp,cm,pos = m
        f.write(','.join(str(x) for x in [chr,snp,cm,pos,a,p]))
        f.write('\n')
