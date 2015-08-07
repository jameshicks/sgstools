#!/usr/bin/env python

import sys
import itertools

from math import factorial

import pandas as pd
import numpy as np

from scipy.stats import chi2

from pydigree.sgs import nshares as shares
from pydigree.io import read_ped # read_germline

matchf,mapf, affectedf, ndf = sys.argv[1:]

#def ecdf(value, dist): 
#    less = sum(x <= value for x in dist)
#    return len(less) / float(len(dist))

def npairs(i):
    return len(i) * (len(i) - 1) / 2
def pvals(v, d):
    vvals = set(v)
    vd = {}
    for val in vvals:
        vd[val] = (d >= val).sum() / float(d.shape[0])
    e = [vd[val] for val in v]
    import ipdb; ipdb.set_trace()
    return np.array(e) 

def shares2(inds, shared, nmark):
    s = np.zeros(nmark)
    for a,b in itertools.combinations(inds,2):
        pair = frozenset([a,b])
        if pair not in shared:
            continue
        for interval in shared[pair]:
            start, stop = interval
            s[start:(stop+1)] += 1
    return s

shares = shares2

print 'Reading null distribution'
with open(ndf) as f:
    nulldist = {}
    for line in f:
        l = line.strip().split()
        fam = l[0]
        values = np.array([float(x) for x in l[1:]])
        nulldist[fam] = values

#with open(affectedf) as f:
#    affecteds = { tuple(x.strip().split()) for x in f }
#    fams = {x[0] for x in affecteds}
peds = read_ped(affectedf)
affecteds = {x for x in peds.individuals if x.phenotypes['affected']}
print '{} affected individuals'.format(len(affecteds))

for ind in peds.individuals:
    if ind.phenotypes['affected'] and ind.is_marryin_founder():
        ind.phenotypes['affected'] = None
#for a in affecteds.copy():
#    if a.is_marryin_founder():
#        print 'Removed affected marry-in founder %s' % a
#        affecteds.remove(a)

affecteds = {x.full_label for x in peds.individuals if x.phenotypes['affected']}
print '{} affecteds after removing marry-in founders'.format(len(affecteds))
fams = {x[0] for x in affecteds}



print 'Reading map'
with open(mapf) as f:
    gmap = [x.strip().split() for x in f]
    positions = [int(x[3]) for x in gmap]
    nmark = len(positions)
    posd =  dict([(y,x) for x,y in enumerate(positions)])

print '{} markers'.format(nmark)
print 

def spairs(affs, shared, nmark):
    return shares(affs, shared, nmark)

def sbool(affs, shared, nmark):
    return shares(affs, shared, nmark) / float(npairs(affs))

score = spairs


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
        haploid = haplocheck(ind1, ind2)

        #if score == spairs and not haploid:
        #    print 'Spairs option requires output from `germline --haploid`'
        #    print 'This file (%s) is not formatted properly.' % args.matchfile
        #    exit(1)

        sharef.seek(0)
            
        for i,line in enumerate(sharef):
            if i % 100000 == 0:
                print 'reading line %s' % i
            l = line.strip().split()

            # Germline haploid output gives which haplotype is shared.
            # If we're just counting up (like we would for spairs) we
            # can just lop off the haplotype identifiers, and the algorithm
            # will just add another 1 when it comes over an overlapping region
            if haploid:
                l[1] = l[1][:-2]
                l[3] = l[3][:-2]


            ind1 = tuple(l[0:2])
            ind2 = tuple(l[2:4])
            
            if ind1[0] != ind2[0]:
                continue
            pair = frozenset({ind1,ind2})

            start = int(l[5])
            stop = int(l[6])

            istart = posd[start]
            istop = posd[stop]
            markersinshare = istop - istart
            #import pdb; pdb.set_trace()
            #print start, stop, istart, istop, markersinshare

            #minsegmentlength = (5 * 10**6)
            #if (stop - start) < minsegmentlength:
            #    continue
            
            #markerdensitylimit = (100 / float(10**6))
            #if (markersinshare / float(stop - start)) < markerdensitylimit:
            #    continue

            if pair not in shared:
                shared[pair] = []
            shared[pair].append((istart, istop))
    return shared

shared = read_germline()#matchf, mapf, haploid=True)


#gmap = pd.read_table(mapfile, sep='\t', names=['chr','snp','cm','pos'])
#nmark = gmap.pos.shape[0]



pvalues = {}
sharestats = {}

#from pydigree.io.fastio import fast_ecdf as ecdf

for fam in fams:
    if fam not in nulldist:
        #print 'No null distribution for pedigree {}, skipping...'.format(fam)
        continue
    affs_in_fam = [x for x in affecteds if x[0] == fam]
    print 'Calculating pedigree {} ({} affecteds)... '.format(fam, len(affs_in_fam)),
    s = score(affs_in_fam, shared, nmark) #/ float(npairs(affs_in_fam))
    #import pdb; pdb.set_trace()
    p = pvals(s, nulldist[fam])
    print 'Maximum sharing: {} (p={})'.format(s.max(),p.min())
    #import pdb; pdb.set_trace()
    sharestats[fam] = s
    pvalues[fam] = p

fishercombine = True
if fishercombine:
    print "Combining p-values by Fisher's method"
    def fisher(pvector):
        k = pvector.shape[0]
        chisq = -2 * (np.log(pvector)).sum()
        return 1 - chi2.cdf(chisq, 2*k)
    o = pd.DataFrame()
    for fam in pvalues.keys():
        o['p_fam{}'.format(fam)] = pvalues[fam]
    p_meta = o.apply(fisher,axis=1)


if True:
    o = pd.DataFrame()
    o['pos'] = positions
    if fishercombine:
        o['p_meta'] = p_meta
    for fam in pvalues.keys():
        o['p_fam{}'.format(fam)] = pvalues[fam]
    o.to_csv('output.pvals', index=False)
if True:
    o = pd.DataFrame()
    o['pos'] = positions
    for fam in sharestats.keys():
        o['s_fam{}'.format(fam)] = sharestats[fam]
    o.to_csv('output.scores', index=False)
