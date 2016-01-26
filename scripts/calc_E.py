# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25th 2015

@author: flamholz
"""

from component_contribution.component_contribution_trainer import ComponentContribution
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.thermodynamic_constants import F, default_T

import argparse
import csv
import numpy as np

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=
				'Calculate reduction potentials for a number of reactions.')
	parser.add_argument('infile', type=argparse.FileType(),
				help='path to input file containing a list of reactions')
	parser.add_argument('outfile', type=argparse.FileType('w'),
				help='path to output file')
	parser.add_argument('--ionic_strength', default=0.2, type=int,
				help='ionic strength in molar units.')
	parser.add_argument('--pH_min', default=5, type=int,
				help='lowest pH to produce E0 for.')
	parser.add_argument('--pH_max', default=9, type=int,
				help='highest pH to produce E0 for.')
	parser.add_argument('--pH_step', default=0.05, type=float,
				help='pH increment.')

	args = parser.parse_args()

	I = args.ionic_strength
	T = default_T

	cc = ComponentContribution.init()
	pHs = np.arange(args.pH_min, args.pH_max + args.pH_step,
					args.pH_step)
	
	reactions_and_energies = []
	reader = csv.reader(args.infile)
	for row in reader:
		formula = row[0].strip()
		reaction = KeggReaction.parse_formula(formula)
		reaction_atom_bag = reaction._get_reaction_atom_bag()
		n_e = reaction_atom_bag.pop('e-', 0)
		if len(reaction_atom_bag) != 0:
			raise ValueError('This is not a half-reaction'
							 ' (i.e. cannot be balanced by adding e-)')
	
		dG0_r, u_r = cc.get_dG0_r(reaction)
		E0s = []
		for pH in pHs:
			ddG0_r = reaction.get_transform_ddG0(pH=pH, I=I, T=T)
			dG0_r_prime = dG0_r + ddG0_r
			E0_prime = 1000 * -dG0_r_prime / (n_e*F) # mV
			E0s.append(E0_prime)

		reactions_and_energies.append((row, E0s))
	
	header = ['reaction', 'reaction_description', 'E\'m', 'Type', 'Source']
	pH_header = ['pH %.1f (mV)' % pH for pH in pHs]
	header += pH_header
	writer = csv.writer(args.outfile)
	writer.writerow(header)
	for rxn_data, pH_E0 in reactions_and_energies:
		energies_fmted = ['%.2f' % e for e in pH_E0]	
		writer.writerow(rxn_data + energies_fmted)


