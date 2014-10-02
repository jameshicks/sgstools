#!/usr/bin/env python

import sys
import itertools

import pandas as pd
import numpy as np

from pydigree.sgs import nshares as shares


matchf,mapf, affectedf, ndf = sys.argv[1:]

#def ecdf(value, dist): 
#    less = sum(x <= value for x in dist)
#    return len(less) / float(len(dist))

def ecdf(v, d):
    vvals = set(v)
    vd = {}
    for val in vvals:
        vd[val] = len([x for x in d if x < val])
    e = [vd[val] for val in v]
    return np.array(e) / float(d.shape[0])
    #return (v < d).sum() / float(d.shape[0])

print 'Reading null distribution'
with open(ndf) as f:
    nulldist = {}
    for line in f:
        fam, vals = line.strip().split(None,1)
        vals = [float(x) for x in vals.split()]
        nulldist[fam] = np.array([x for x in vals])

with open(affectedf) as f:
    affecteds = { tuple(x.strip().split()) for x in f }
    fams = {x[0] for x in affecteds}

print 'Reading map'
with open(mapf) as f:
    gmap = [x.strip().split() for x in f]
    positions = [int(x[3]) for x in gmap]
    nmark = len(positions)
    posd =  dict([(y,x) for x,y in enumerate(positions)])

print '{} markers'.format(nmark)
print 

def haplocheck(a,b):
    return a.endswith('.0') or a.endswith('.1') and \
           b.endswith('.0') or b.endswith('.1')

def read_germline():
    with open(matchf) as sharef:
        shared = {}
        
        # Test the first line to see if we're in a haploid file
        line = sharef.readline()
        l = line.strip().split()
        ind1 = '.'.join(l[0:2])
        ind2 = '.'.join(l[2:4])
#        haploid = haplocheck(ind1, ind2)

#        if args.model == 'Spairs' and not haploid:
#            print 'Spairs option requires output from `germline --haploid`'
#            print 'This file (%s) is not formatted properly.' % args.matchfile
#            exit(1)

 #       sharef.seek(0)
            
        for i,line in enumerate(sharef):
            if i % 100000 == 0:
                print 'reading line %s' % i
            l = line.strip().split()

            # Germline haploid output gives which haplotype is shared.
            # If we're just counting up (like we would for spairs) we
            # can just lop off the haplotype identifiers, and the algorithm
            # will just add another 1 when it comes over an overlapping region
#            if haploid:
#                l[1] = l[1][:-2]
#                l[3] = l[3][:-2]


            ind1 = tuple(l[0:2])
            ind2 = tuple(l[2:4])
            
            pair = frozenset({ind1,ind2})

            start = int(l[5])
            stop = int(l[6])

            istart = posd[start]
            istop = posd[stop]
            markersinshare = istop - start

            minsegmentlength = (.5 * 10**6)
            if (stop - start) < minsegmentlength:
                continue
            
            markerdensitylimit = (100 / float(10**6))
            if (markersinshare / float(stop - start)) < markerdensitylimit:
                continue

            if pair not in shared:
                shared[pair] = []
            shared[pair].append((istart, istop))
    return shared

shared = read_germline()


#gmap = pd.read_table(mapfile, sep='\t', names=['chr','snp','cm','pos'])
#nmark = gmap.pos.shape[0]


pvals = {}
sharestats = {}

#from pydigree.io.fastio import fast_ecdf as ecdf

for fam in fams:
    if fam not in nulldist:
        print 'No null distribution for pedigree {}, skipping...'.format(fam)
        continue
    else:
        print 'Calculating pedigree {}'.format(fam)
    affs_in_fam = [x for x in affecteds if x[0] == fam]
    s = shares(affs_in_fam, shared, nmark)
    p = 1 - np.array(ecdf(s, nulldist[fam]))
    sharestats[fam] = s
    pvals[fam] = p

if True:
    o = pd.DataFrame()
    o['pos'] = positions
    for fam in pvals.keys():
        o['p_fam{}'.format(fam)] = pvals[fam]
    o.to_csv('output.pvals', index=False)
