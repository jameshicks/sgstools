import argparse

from bisect import bisect_left
from math import log, log10, e

import numpy as np

import pydigree

from pydigree.common import log_base_change
from pydigree.mixedmodel import MixedModel


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--ped')
parser.add_argument('--map')
parser.add_argument('--phen')
parser.add_argument('--sgs')
parser.add_argument('--outcome')
parser.add_argument('--fixefs')
parser.add_argument('--every', default=0.5, type=float, 
                    help='Distance between evaluations (in Mb)')
args = parser.parse_args()


# Read pedigree data
peds = pydigree.io.read_ped(args.ped)

# Read genotype map data
chroms = pydigree.io.plink.read_map(args.map)
for chrom in chroms:
    peds.add_chromosome(chrom)

# Read phenotype data
pydigree.io.read_phenotypes(peds, args.phen)

# Get valid individuals from phenotypes
analysis_individuals = [x for x in peds.individuals
                        if args.outcome in x.phenotypes]

# Read SGS data
sgs = pydigree.io.sgs.read_germline(args.sgs)
sgs.update_segment_references(peds)

print 'Fitting polygenic model'
null_model = MixedModel(peds, outcome=args.outcome, fixefs=args.fixefs)
null_model.add_genetic_effect()
null_model.fit_model()
null_model.maximize()
llik_null = null_model.full_loglikelihood()
print 'Done'


def vc_linkage(null_model, locus):
    ibd_model = null_model.copy()
    ibdmat = sgs.ibd_matrix(sgs, individuals, locus, location_type='physical')
    ibd_model.add_random_effect('IBD', analysis_individuals, ibdmat)

    ibd_model.fit()
    ibd_model.maximize()
    return ibd_model

for chromidx, chromosome in enumerate(peds.chromosomes):
    pstart, pstop = chromosome.physical_map[0], chromosome.physical_map[-1]
    evaluation_sites = np.arange(pstart, pstop, args.every)

    for evaluation_site in evaluation_sites:
        markidx = chromosome.closest_marker(evaluation_site)
        locus = chromidx, markidx
        ibd_model = vc_linkage(null_model, locus)
        llik_ibd = ibd_model.full_loglikelihood()

        lod = (log_base_change(llik_ibd, e, 10) -
               log_base_change(llik_null, e, 10))

        print '\t'.join([chromosome.label, evaluation_site, lod])
