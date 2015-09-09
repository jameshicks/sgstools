import argparse

from math import log, log10, e

import numpy as np
from scipy.stats import chi2

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
parser.add_argument('--onlywithin', action='store_true')
parser.add_argument('--sites', nargs='*', type=int, dest='evalsites')
parser.add_argument('--method', default='FS', dest='maxmethod')
parser.add_argument('--verbose', action='store_true', help='Show progress of maximizer')
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
null_model.maximize(method=args.maxmethod, verbose=args.verbose)
null_model.summary()
llik_null = null_model.loglikelihood()
print 'Done'


def vc_linkage(locus):
    ibd_model = MixedModel(peds, outcome=args.outcome, fixed_effects=args.fixefs)
    ibd_model.add_genetic_effect()
    ibdmat = sgs.ibd_matrix(analysis_individuals, locus, location_type='index', onlywithin=args.onlywithin)

    ranef = RandomEffect(analysis_individuals,
                         'IBD',
                         incidence_matrix='eye',
                         covariance_matrix=ibdmat)
    ibd_model.add_random_effect(ranef)

    ibd_model.fit_model()
    ibd_model.maximize(verbose=args.verbose, method=args.maxmethod)
    return ibd_model

class VCLResult(object):
    def __init__(self, alternative, null):
        self.llik_alt = alternative.loglikelihood()
        self.llik_null = null.loglikelihood()

    @property
    def chisq(self):
        return -2.0 * self.llik_null + 2.0 * self.llik_alt 
    
    @property
    def pvalue(self):
        return 1 - chi2.cdf(self.chisq, 1)
    
    @property
    def lod(self):
        return self.chisq / (2.0 * log(10.0))

outputlist = []
for chromidx, chromosome in enumerate(peds.chromosomes):
    pstart, pstop = chromosome.physical_map[0], chromosome.physical_map[-1]
    if not args.evalsites:
        evaluation_sites = np.arange(pstart, pstop, args.every * 1e6)
    else:
        evaluation_sites = args.evalsites
    for evaluation_site in evaluation_sites:
        markidx = chromosome.closest_marker(evaluation_site)
        locus = chromidx, markidx
        
        ibd_model = vc_linkage(locus)
        llik_ibd = ibd_model.loglikelihood()
        vc = VCLResult(ibd_model, null_model)
        lod = (log_base_change(llik_ibd, e, 10) -
               log_base_change(llik_null, e, 10))
    
       # lod = np.nan
    
        output = [chromosome.label, chromosome.physical_map[markidx], vc.lod, vc.pvalue]
        outputlist.append(output)
        print '\t'.join(str(x) for x in output)


print 'Done!'