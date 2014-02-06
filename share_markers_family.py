#!/usr/bin/env python

# Usage: <script> sharefile mapfile affectedlist
import sys
from math import factorial
from itertools import combinations,izip

sep = ','
def pair_shares(pair,position):
    if pair not in shares: return False
    for x,y in shares[pair]:
        if x <= position <= y: return True
    return False
def nCk(n,k): return factorial(n) / (factorial(k) * factorial(n-k))
def wmean(values,weights):
    return sum(val*weight for val,weight in izip(values,weights))/sum(weights)
def mean(values): return wmean(values,[1]*len(values))


with open(sys.argv[3]) as f:
    affinds = set(tuple(x.strip().split()) for x in f)
    aff_pairs = [x for x in combinations(affinds,2)]
    num_aff_pairs = len(aff_pairs)
    families = sorted(frozenset(x[0] for x in affinds))
with open(sys.argv[2]) as f:
    map = [x.strip().split() for x in f]
    positions = [int(x[3]) for x in map]

with open(sys.argv[1]) as sharef:
    shares = {}
    for line in sharef:
        l = line.strip().split()
        # Skip pairs not from the same family
        if not l[0] == l[2]: continue
        ind1 = tuple(l[0:2])
        ind2 = tuple(l[2:4])
        pair = frozenset([ind1,ind2])
        if not pair <= affinds: continue
        start,stop = [int(x) for x in l[5:7]]
        if pair not in shares:
            shares[pair] = []
        shares[pair].append([start,stop])

print sep.join(['chr','SNP','cM','pos','meanshareprecent','wmeansharepercent'] + ['FAM'+x for x in families])
for chr,snp,cm,pos in map:
    pos = int(pos)
    per_family_shares = []
    for family in families:
        faminds = [x for x in affinds if x[0] == family]
        fam_pairs = [frozenset(x) for x in combinations(faminds,2)]
        does_share = [pair_shares(x,pos) for x in fam_pairs]
        sharecounts = sum(1 for x in does_share if x)
        numpairs = float(len(fam_pairs))
        per_family_shares.append([numpairs,sharecounts/numpairs])
    weights,sharep = zip(*per_family_shares)
    mshare = mean(sharep)
    wmshare = wmean(sharep,weights)
    print sep.join(str(x) for x in [chr,snp,cm,pos,mshare,wmshare]+list(sharep))
