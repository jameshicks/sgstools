import argparse

from math import log, log10, e

import numpy as np

import pydigree

from pydigree.common import log_base_change
from pydigree.io.sgs import read_germline
from pydigree.mixedmodel import MixedModel, RandomEffect


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--ped')
parser.add_argument('--map')
parser.add_argument('--phen')
parser.add_argument('--sgs')
parser.add_argument('--outcome', required=True)
parser.add_argument('--fixefs')
parser.add_argument('--every', default=0.5, type=float,
                    help='Distance between evaluations (in Mb)')
args = parser.parse_args()


# Read pedigree data
print 'Reading pedigree'
peds = pydigree.io.read_ped(args.ped)

# Read genotype map data
print 'Reading map'
chroms = pydigree.io.plink.read_map(args.map)
for chrom in chroms:
    peds.add_chromosome(chrom)

# Read phenotype data
print 'Reading phenotype data'
pydigree.io.read_phenotypes(peds, args.phen)

# Get valid individuals from phenotypes
analysis_individuals = [x for x in peds.individuals
                        if args.outcome in x.phenotypes]

# Read SGS data
print 'Reading SGS data'
sgs = read_germline(args.sgs)

print 'Updating references'
sgs.update_segment_references(peds)

print 'Fitting polygenic model'
null_model = MixedModel(peds, outcome=args.outcome, fixed_effects=args.fixefs)
null_model.add_genetic_effect()
null_model.fit_model()
null_model.maximize()
null_model.summary()
llik_null = null_model.loglikelihood()
print 'Done'


def vc_linkage(null_model, locus):
    ibd_model = MixedModel(peds, outcome=args.outcome, fixed_effects=args.fixefs)
    ibd_model.add_genetic_effect()
    ibdmat = sgs.ibd_matrix(analysis_individuals, locus, location_type='index')
    ranef = RandomEffect(analysis_individuals,
                         'IBD',
                         incidence_matrix='eye',
                         covariance_matrix=ibdmat)
    ibd_model.add_random_effect(ranef)

    ibd_model.fit_model()
    ibd_model.maximize()
    return ibd_model

for chromidx, chromosome in enumerate(peds.chromosomes):
    pstart, pstop = chromosome.physical_map[0], chromosome.physical_map[-1]
    evaluation_sites = np.arange(pstart, pstop, args.every * 1e6)

    for evaluation_site in evaluation_sites:
        markidx = chromosome.closest_marker(evaluation_site)
        locus = chromidx, markidx
        ibd_model = vc_linkage(null_model, locus)
        llik_ibd = ibd_model.loglikelihood()

        lod = (log_base_change(llik_ibd, e, 10) -
               log_base_change(llik_null, e, 10))

        import ipdb; ipdb.set_trace()
        output = [chromosome.label, evaluation_site, lod]
        print '\t'.join(str(x) for x in output)
