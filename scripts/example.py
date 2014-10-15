# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 18:32:46 2014

@author: noore
"""
from scripts.max_min_driving_force import KeggFile2ModelList, MaxMinDrivingForce
from python.component_contribution import ComponentContribution
from scripts.html_writer import HtmlWriter

REACTION_FNAME = 'scripts/example_pathway.txt'
HTML_FNAME = 'scripts/example_results.html'

html_writer = HtmlWriter(HTML_FNAME)
pathways = KeggFile2ModelList(REACTION_FNAME)
p = pathways[0]
cc = ComponentContribution.init()


p['model'].add_thermo(cc)

mdf = MaxMinDrivingForce(p['model'], p['fluxes'], p['bounds'],
                         pH=p['pH'], I=p['I'], T=p['T'],
                         html_writer=html_writer)

mdf_solution, dG_r_prime = mdf.Solve(uncertainty_factor=1.0)
print mdf_solution