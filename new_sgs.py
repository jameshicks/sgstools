import argparse
import pydigree as pyd
from pydigree.sgs import sgs_population
from pydigree.io.sgs import write_sgs

parser = argparse.ArgumentParser()
parser.add_argument('--plink', required=True,
                    help='Prefix of plink formatted PED/MAP files')
parser.add_argument('--out', required=True, help='File name for outfit')
parser.add_argument('--only-within-peds', action='store_true', dest='onlywithin',
                    help='Only perform SGS within pedigrees')
parser.add_argument('--minsize', type=float, default=1.0,
                    help='Minimum size of segment (in megabases)')
parser.add_argument('--seedsize', dest='seedsize', default=1000,
                    help='Seed size required for shared segment')
args = parser.parse_args()

peds = pyd.io.plink.read_plink(prefix=args.plink)

analysis = sgs_population(peds,
                          onlywithin=args.onlywithin,
                          seed_size=args.seedsize,
                          min_length=args.minsize,
                          size_unit='mb')
import ipdb; ipdb.set_trace()
pyd.io.sgs.write_sgs(analysis, args.out)
