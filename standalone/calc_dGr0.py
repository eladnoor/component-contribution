#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 15:36:47 2017

@author: noore

A stand-alone version of Component Contribution that can calculate the 
Delta-Gr'0 of any reaction (with KEGG notation, i.e. whose reactants are
already cached in our database), at a given pH and I.
"""

import argparse
import logging
from preprocessing import Preprocessing, Reaction

def MakeParser():
    parser = argparse.ArgumentParser(
        description=('Estimate the Gibbs energy of a reaction. For example,'
                     'the following calculates dGr0 for ATP hydrolysis '
                     'at pH 6: calc_dGr0.py --ph 6 "C00002 + C00001 = '
                     'C00008 + C00009"'))
    parser.add_argument('--ph', type=float, help='pH level', default=7.0)
    parser.add_argument('--i', type=float,
                        help='ionic strength in M',
                        default=0.1)
    parser.add_argument('reaction', type=str, help='reaction in KEGG notation')
    return parser


###############################################################################
parser = MakeParser()
args = parser.parse_args()

logging.getLogger().setLevel(logging.WARNING)

print 'pH: %.1f' % args.ph
print 'I: %.1f M' % args.i
print 'Reaction: ' + args.reaction

# parse the reaction
reaction = Reaction.parse_formula(args.reaction)

p = Preprocessing()

# use the preprocessing class to calculate the estimated dG0 and uncertainty
dG0_prime, U = p.dG0_prime(reaction, pH=args.ph, I=args.i)
print u'dGr0 = %.2f \u00B1 %.2f kJ/mol' % (dG0_prime[0, 0], U[0, 0])
