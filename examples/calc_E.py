# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 14:58:44 2014

@author: eladn
"""
import sys, logging, os
import numpy as np
from python.component_contribution import ComponentContribution
from python.kegg_reaction import KeggReaction
from python.compound_cacher import CompoundCacher
from python.training_data import TrainingData

#pH = 7
I = 0.2
T = 298.15
F = 96.48 / 1000.0 # kJ/mol / mV
ccache = CompoundCacher()
CC_CACHE_FNAME = 'cache/component_contribution.mat'
if not os.path.exists(CC_CACHE_FNAME):
    #TrainingData.FNAME_DICT['FORMATION'] = \
    #    ('../data/formation_energies_transformed.tsv', 10.0)
    td = TrainingData()
    cc = ComponentContribution(td)
    cc.save_matfile(CC_CACHE_FNAME)
else:
    cc = ComponentContribution.from_matfile(CC_CACHE_FNAME)

#formula = 'C00033 + C00282 <=> C00084 + C00001'
formula = 'C00067 + C00001 <=> C00058 + C00010'
#formula = 'C00003 + C00282 <=> C00004'
#############################################################################
reaction = KeggReaction.parse_formula(formula)
reaction_atom_bag = reaction._get_reaction_atom_bag()
n_e = 0
if 'e-' in reaction_atom_bag:
    n_e = reaction_atom_bag['e-']
    del reaction_atom_bag['e-']
if len(reaction_atom_bag) != 0:
    raise Exception('Reaction is not balanced (not only with respect to e-)')

dG0_r, u_r = cc.get_dG0_r(reaction)

for pH in xrange(6, 9):
    ddG0_r = reaction.get_transform_ddG0(pH=pH, I=I, T=T)
    dG0_r_prime = dG0_r+ddG0_r
    E0_prime = -dG0_r_prime / (n_e*F)
    print 'pH = %4.1f, E\'0 = %6.0f' % (pH, E0_prime)


"""
dG0_r_prime = 0
for cid, coeff in reaction.iteritems():
    print '****', cid, '****'
    comp = ccache.get_compound(cid)
    dG0_f = cc.get_major_ms_dG0_f(cid)
    print cc.get_compound_json(cid)

    #print comp._dG0_prime_vector(pH, I, T)
    #print dG0_f + comp.transform_pH7(pH=pH, I=I, T=T)
    #print list(comp.get_species(dG0_f, T))
    ddG0_f = comp.transform_pH7(pH=pH, I=I, T=T)
    #print cc.get_compound_json(cid)
    print 'dG0 = %.1f, dG0_prime = %.1f' % (dG0_f, dG0_f+ddG0_f)
    dG0_r_prime += coeff*(dG0_f+ddG0_f)
print 'total dG0_prime = %.1f' % dG0_r_prime
"""


