#!/usr/bin/env python

import sys

with open(sys.argv[1]) as f:
    for line in f:
        print ' '.join(line.strip().split()[0:7])
