#!/usr/bin/env python

# Usage: <script> sharefile mapfile affectedlist
import sys
from itertools import combinations
from math import factorial

def nCk(n,k):
    return factorial(n) / (factorial(k) * factorial(n-k))
def pair_shares(pair,position):
    for x,y in shares[pair]:
        if x <= position <= y: return True
    return False

with open(sys.argv[3]) as f:
    affinds = set('.'.join(x.strip().split()) for x in f)
    num_aff_pairs = nCk(len(affinds),2)

with open(sys.argv[2]) as f:
    map = [x.strip().split() for x in f]
    positions = [int(x[3]) for x in map]

with open(sys.argv[1]) as sharef:
    shares = {}
    for line in sharef:
        l = line.strip().split()
        ind1 = '.'.join(l[0:2])
        ind2 = '.'.join(l[2:4])
        pair = frozenset([ind1,ind2])
        if not pair <= affinds: continue
        start,stop = [int(x) for x in l[5:7]]
        if pair not in shares:
            shares[pair] = []
        shares[pair].append([start,stop])
share_keys = [x for x in shares.keys()]

print ','.join(['chr','SNP','cM','pos','nshares','percentshares'])
for chr,snp,cm,pos in map:
    pos = int(pos)
    does_share = [pair_shares(x,pos) for x in share_keys]
    sharecounts = sum(1 for x in does_share if x)
    print ','.join([str(x) for x in [chr,snp,cm,pos,sharecounts,(sharecounts/float(num_aff_pairs))]])
