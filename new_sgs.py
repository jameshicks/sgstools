import argparse
from itertools import combinations
from math import factorial

import pydigree as pyd
from pydigree.sgs import SGSAnalysis
from pydigree.io.sgs import write_sgs

parser = argparse.ArgumentParser()
parser.add_argument('--ped', required=True, dest='plinkped',
                    help='Prefix of plink formatted PED file')
parser.add_argument('--map', required=True, dest='plinkmap',
                    help='Prefix of plink formatted MAP file')
parser.add_argument('--out', required=True, help='File name for outfit')
parser.add_argument('--only-within-peds', action='store_true',
                    dest='onlywithin',
                    help='Only perform SGS within pedigrees')
parser.add_argument('--minsize', type=float, default=1.0,
                    help='Minimum size of segment (in megabases)')
parser.add_argument('--seedsize', dest='seedsize', default=1000,
                    help='Seed size required for shared segment')
parser.add_argument('--njobs', '-j', dest='njobs', type=int, default=1,
                    help='Number of processes to use for parallel computation')
args = parser.parse_args()

print 'Reading data'
peds = pyd.io.plink.read_plink(pedfile=args.plinkped, mapfile=args.plinkmap)
print 'Performing SGS analysis'

analysis = SGSAnalysis.direct_to_disk(args.out, peds, njobs=args.njobs,
                                      onlywithin=args.onlywithin,
                                      seed_size=args.seedsize,
                                      min_length=args.minsize,
                                      size_unit='mb')


