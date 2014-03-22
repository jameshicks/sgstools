#!/usr/bin/env python

import argparse
import random

from itertools import combinations, izip, imap
from multiprocessing import Pool

import numpy as np 


parser = argparse.ArgumentParser()

group = parser.add_argument_group('Files')
group.add_argument('--match', metavar='FILE', required=True, dest='matchfile',
                    help='GERMLINE match file')
group.add_argument('--map', metavar='FILE', required=True, dest='mapfile',
                    help='PLINK format map file')
group.add_argument('-o', '--out', metavar='FILE', required=True, dest='outfile',
                    help='Output file')
group.add_argument('--affecteds', metavar='FILE', required=True, dest='afffile',
                    help='File with list of affecteds')
group.add_argument('--population', metavar='FILE', required=True, dest='popfile',
                    help='File with list of population individuals')
group.add_argument('--kinship', default=None, help='Optional kinship file for matching')

group = parser.add_argument_group('Options for resampling based p-values')
group.add_argument('--matchkinship', action='store_true', default=False,
                   help='Match for mean kinship when drawing null population samples')
group.add_argument('-j','--nproc', metavar='n', default=1, action='store', type=int,
                    help='Parallelize P value computation by spreading across'
                    'n processes', dest='njobs')
group.add_argument('--nrep', default=10**5, type=int, metavar='n', dest='nrep',
                    help='Number of replications for emperical p-values')
group.add_argument('--seed', type=int, action='store', default=0,
                    help='Seed value for random number generator')

group = parser.add_argument_group('QC for shared segments')
group.add_argument('--density', default=(100 / float(10**6)), metavar='d',
                    dest='markerdensitylimit', type=float,
                    help='Minimum number of markers per basepair to include'
                    'segment in analysis, default 100 markers/Mb')
group.add_argument('--minsegment', default=(.5 * 10**6), type=float, metavar='l',
                    dest='minsegmentlength',
                    help='Minimum size (in Mb) for segment to be included in analysis') 


args = parser.parse_args()


# Python implementation of shares function if pydigree is not available
def shares_py(inds, shared, nmark):
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

try:
    from pydigree.sgs import proportion_shares as shares
except:
    print "Could not find module 'pydigree', using slower pure python implementation"
    shares = shares_py


### Kinship functions for matching

def get_kinship(a,b):
    try:
        return kinship[frozenset({a,b})]
    except KeyError:
        return 0.0

def mean_kinship(inds):
    ninds  = len(inds)
    return sum(get_kinship(a,b) for a,b in combinations(inds,2)) / (ninds * (ninds-1) * 0.5)

def match_kinship_sample(affs, population, threshold=0.001):
    ninds = len(affs)
    ak = mean_kinship(affs)
    for x in xrange(10000):
        ns = random.sample(population, ninds)
        if abs(mean_kinship(ns) - ak) < threshold:
            return ns
    else:
        print 'ERROR: Took to many tries to match!'
        exit(1)

def numpairs(n):
    return n * (n-1) * 0.5

if args.seed:
    random.seed(args.seed)
if args.matchkinship  and not args.kinship:
    print 'ERROR: No kinship file specified'
    exit(1)

print 'Minimum segment length: %sMb' % (args.minsegmentlength / float(10**6))
print 'Minimum marker density: %s markers/Mb' % (args.markerdensitylimit *  float(10**6))

print 'Reading individual lists'
with open(args.afffile) as f:
    affinds = set('.'.join(x.strip().split()) for x in f)
    naff = len(affinds)
with open(args.popfile) as f:
    fullinds = set('.'.join(x.strip().split()) for x in f)

print 'Reading map'
with open(args.mapfile) as f:
    gmap = [x.strip().split() for x in f]
    positions = [int(x[3]) for x in gmap]
    nmark = len(positions)
    posd =  dict([(y,x) for x,y in enumerate(positions)])

if naff == 0:
    print 'ERROR: No affected individuals!'
    exit(1)

if affinds - fullinds:
    print "%d individuals in the affected list weren't in the population list" % len(affinds - fullinds)
    print "These individuals were removed"
    affinds = affinds & fullinds

print 'Reading Share file'
with open(args.matchfile) as sharef:
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

        if (stop - start) < args.minsegmentlength:
            continue

        if (markersinshare / float(stop - start)) < args.markerdensitylimit:
            continue
        
        if pair not in shared:
            shared[pair] = []
        shared[pair].append((istart, istop))
    keyset = frozenset(shared.keys())

kinship = {}
if args.kinship and args.matchkinship:
    print 'Reading kinship'
    with open(args.kinship) as kinshf:
        for line in kinshf:
            fam, ida, idb, phi = line.strip().split()
            ida = (fam, ida)
            idb = (fam, idb)
            if not {ida, idb} <= fullinds:
                continue
            else:
                kinship[frozenset({ida, ibd})] = float(phi)


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
affshare = shares(affinds, shared, nmark)


print 'Calculating emperical p-values from %s draws of %s individuals' % (args.nrep, naff)
print 'Using %s processes' % args.njobs

def nsharehelper(x):
    if x % 1000 == 0:
        print 'Random draw %s' % x 
    if args.matchkinship:
        ns = match_kinship_sample(affinds, fullinds, threshold=0.001)
    else:
        ns = random.sample(fullinds, naff)
    return shares(ns, shared, nmark)

if args.njobs > 1:
    pool = Pool(processes=args.njobs)
    nullshares = pool.imap_unordered(nsharehelper, xrange(args.nrep), chunksize=500)
else:
    nullshares = imap(nsharehelper, xrange(args.nrep))

pvals = sum(n >= affshare for n in nullshares) / float(args.nrep)

print 'Minimum observed P: %s' % min(pvals)
print
print 'Writing output'

with open(args.outfile,'w') as f:
    f.write(','.join(['chr','snp','cm','pos','pctshares', 'p']) + '\n')
    for m,a,p in izip(gmap, affshare, pvals):
        chr,snp,cm,pos = m
        f.write(','.join(str(x) for x in [chr,snp,cm,pos,a,p]))
        f.write('\n')
