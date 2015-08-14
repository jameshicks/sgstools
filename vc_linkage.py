import argparse

from bisect import bisect_left
from math import log, log10, e

import pydigree

from pydigree.common import log_base_change
from pydigree.mixedmodel import MixedModel
from pydigree.sgs import ibd_matrix

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--ped')
parser.add_argument('--map')
parser.add_argument('--phen')
parser.add_argument('--sgs')
parser.add_argument('--outcome')
parser.add_argument('--fixefs')
parser.add_argument('--every', help='Distance between evaluations')
args = parser.parse_args()

# Support functions


def read_map(filename):
    pass


# Read pedigree data
peds = pydigree.io.read_ped(args.ped)
# TODO: Read genotype map data

# Read phenotype data
pydigree.io.read_phenotypes(peds, args.phen)

# Get valid individuals from phenotypes
analysis_individuals = [x for x in peds.individuals
	if args.outcome in x.phenotypes]

# TODO: Read SGS data
sgs = pydigree.io.sgs.read_germline(args.sgs)
sgs.update_segment_references(peds)

print 'Fitting polygenic model'
null_model = MixedModel(peds, outcome=args.outcome, fixefs=args.fixefs)
null_model.add_genetic_effect()
null_model.fit_model()
null_model.maximize()
llik_null = null_model.restricted_loglikelihood()
print 'Done'


# TODO: Find evaluation sites
def find_site(target, values):
    pass


def vc_linkage(null_model, locus):
	# TODO: Find evaluation site's locus
    eval_site = None
    ibd_model = null_model.copy()
    ibdmat = sgs.ibd_matrix(sgs, individuals, locus, location_type='physical')
    ibd_model.add_random_effect('IBD', analysis_individuals, ibdmat)

    ibd_model.fit()
    ibd_model.maximize()
    return ibd_model

for x in xrange(evaluation_sites):

	ibd_model = vc_linkage(null_model, x)
    llik_ibd = ibd_model.restricted_loglikelihood()

    lod = (log_base_change(llik_ibd, e, 10) - log_base_change(llik_null, e, 10)
    
    print '\t'.join([chrom, pos, lod])
