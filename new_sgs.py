import argparse
from itertools import combinations
from math import factorial

import pydigree as pyd
from pydigree.sgs import sgs_population
from pydigree.io.sgs import write_sgs

parser = argparse.ArgumentParser()
parser.add_argument('--ped', required=True, dest='plinkped',
                    help='Prefix of plink formatted PED file')
parser.add_argument('--map', required=True, dest='plinkmap',
                    help='Prefix of plink formatted MAP file')
parser.add_argument('--out', required=True, help='File name for outfit')
parser.add_argument('--only-within-peds', action='store_true', dest='onlywithin',
                    help='Only perform SGS within pedigrees')
parser.add_argument('--minsize', type=float, default=1.0,
                    help='Minimum size of segment (in megabases)')
parser.add_argument('--seedsize', dest='seedsize', default=1000,
                    help='Seed size required for shared segment')
args = parser.parse_args()

print 'Reading data' 
peds = pyd.io.plink.read_plink(pedfile=args.plinkped, mapfile=args.plinkmap)
print 'Performing SGS analysis'

analysis = sgs_population(peds,
                          onlywithin=args.onlywithin,
                          seed_size=args.seedsize,
                          min_length=args.minsize,
                          size_unit='mb')

print 'Writing output to {}'.format(args.out)
pyd.io.sgs.write_sgs(analysis, args.out)
