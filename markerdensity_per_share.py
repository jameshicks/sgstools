#!/usr/bin/env python

import sys

mindensity = float(sys.argv[2])

with open(sys.argv[1]) as f:
    for line in f:
        l = line.strip().split()
        start,stop = [int(x) for x in l[5:7]]
        nsnps = int(l[9])
        density = nsnps / ((stop-start) / float(10**6))
        if density >= mindensity: print line.strip()
