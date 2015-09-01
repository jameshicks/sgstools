#!/usr/bin/env python

import argparse
import sys
from itertools import combinations

import pydigree as pyd
from pydigree.sgs import sgs_unphased, make_intervals
from pydigree.io.plink import write_map

parser = argparse.ArgumentParser()
parser.add_argument('--vcf')
parser.add_argument('--ped')
parser.add_argument('--plink')
parser.add_argument('--minmb', default=.5, type=float)
parser.add_argument('--maxmissrate', default=0.0, type=float)
parser.add_argument('--seedsize', default=None, type=int)
parser.add_argument('--out', default=sys.stdout)

args = parser.parse_args()
seedsize = args.seedsize


def seg2germline(seg):
    oline = list(seg.ind1.full_label) + \
        list(seg.ind2.full_label) + [seg.chromosome.label]
    oline += list(seg.physical_location)
    oline += [seg.start, seg.stop, seg.nmark, seg.physical_size / 1e6, 'MB']
    oline += ['X'] * 3
    return '\t'.join([str(x) for x in oline])


print 'Reading pedigree'
ped = pyd.io.read_ped(args.ped)
print 'Done'
if len([x.label for x in ped.individuals]) < len({x.label for x in ped.individuals}):
    print 'Individual identifiers not unique'

if args.vcf:
    print 'Reading VCF'
    genodata = pyd.io.read_vcf(
        args.vcf, sparse=True, geno_missrate=args.maxmissrate)
    if not seedsize:
        seedsize = 10000
elif args.plink:
    print 'Reading plink'
    genodata = pyd.io.read_plink(prefix=args.plink,
                                 onlyinds={x.full_label for x in ped.individuals})
    if not seedsize:
        seedsize = 1000
print 'Done'

ped.merge(genodata)
genoedinds = [x for x in ped.individuals if x.has_genotypes()]
ngenoed = len([x for x in ped.individuals if x.has_genotypes()])
print '{} genotyped individuals'.format(ngenoed)
print '{} combinations'.format(ngenoed * (ngenoed-1) / 2)
wi_combos = len([(x, y) for x, y in combinations(genoedinds, 2)
                 if x.pedigree.label == y.pedigree.label])
print '{} within-pedigree combinations'.format(wi_combos)
bw_combos = len([(x, y) for x, y in combinations(genoedinds, 2)
                 if x.pedigree.label != y.pedigree.label])
print '{} between-pedigree combinations'.format(bw_combos)

for chrom in ped.chromosomes:
    print 'Chromosome {}: {} markers'.format(chrom.label, chrom.nmark())

print 'Writing filtered map file'
write_map(ped, args.out + '.map')

with open(args.out + '.match', 'w') if args.out is not sys.stdout else sys.stdout as o:
    for i, pair in enumerate([(x, y) for x, y in combinations(genoedinds, 2)]):
        ind1, ind2 = pair
        if i % 1000 == 0:
            print 'Pair {}'.format(i)
        for chridx, chrom in enumerate(ped.chromosomes):
            res = sgs_unphased(ind1, ind2, chridx, seed_size=seedsize,
                               min_length=args.minmb, size_unit='mb', min_density=.1)
            for result in res:
                o.write(seg2germline(result) + '\n')
