#!/usr/bin/env python

import sys

family = sys.argv[2]

with open(sys.argv[1]) as f:
    for line in f:
        l = line.strip().split()
        if l[0] == family and l[2] == family:
            print line.strip()
