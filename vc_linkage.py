import argparse

from bisect import bisect_left

import pydigree

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

# TODO: Read pedigree data
# TODO: Read genotype map data
# TODO: Read phenotype data
# TODO: Get valid individuals from phenotypes
analysis_individuals = None
# TODO: Read SGS data

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

for x in xrange(evaluation_sites):
	# TODO: Find evaluation site's locus 
	eval_site = None
	ibd_model = null_model.copy()
	ibdmat = ibd_matrix(sgs, individuals)
	ibd_model.add_random_effect('IBD', analysis_individuals, eval_site)

	ibd_model.fit()
	ibd_model.maximize()
	llik_ibd = ibd_mode.restricted_loglikelihood()

	lod = llik_ibd - llik_null
	print '\t'.join([chrom, pos, lod])


