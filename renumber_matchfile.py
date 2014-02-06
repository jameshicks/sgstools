#!/usr/bin/env python

import sys

matchf,analysisindsf = sys.argv[1:3]

with open(analysisindsf) as f:
    d = {}
    for line in f:
        fam,ind = line.strip().split()[0:2]
        d[ind]=fam

with open(matchf) as f:
    for line in f:
        l = line.strip().split()
        if l[1] not in d or l[3] not in d: continue
        l[0] = d[l[1]]
        l[2] = d[l[3]]
        if l[0] != l[2]: continue
        print '\t'.join(l)
