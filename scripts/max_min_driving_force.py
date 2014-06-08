# -*- coding: utf-8 -*-
"""
Created on Wed May 28 13:13:37 2014

@author: eladn
"""

import numpy as np
import matplotlib.pyplot as plt
from python.thermodynamic_constants import R, default_c_range, default_c_mid
from python.component_contribution import ComponentContribution
from python.kegg_model import KeggModel

from scripts.pathway_modeling import KeggPathway, DeltaGNormalization
from scripts.html_writer import HtmlWriter
    
class ConcentrationConstraints(object):
    pass

class MaxMinDrivingForce(object):
    
    def __init__(self, model, fluxes, pH, I, T):
        self.fluxes = np.matrix(fluxes)
        self.pH, self.I, self.T = pH, I, T
        self.c_range = default_c_range
        self.c_mid = default_c_mid

        self.html_writer = HtmlWriter
        self.html_writer.write('Parameters:</br>\n')
        condition_list = ['pH = %g' % self.pH,
                          'Ionic strength = %g M' % self.I,
                          'Temperature = %g K' % self.T,
                          'Concentration range = %g - %g M' % self.c_range,
                          'Default concentration = %g M' % self.c_mid]
        self.html_writer.write_ul(condition_list)

    def GetBounds(self):
        cid2bounds = {cid:self.c_range for cid in self.model.cids}
        cid2bounds['C00001'] = (1, 1) # the default for H2O is 1

    def Solve(self, prefix='mdf'):
        S = self.model.S
        rids = self.model.rids
        cids = self.model.cids
        dG0_prime, dG0_std = model.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
        cid2bounds = self.GetBounds()

        keggpath = KeggPathway(S, rids, self.fluxes, cids, None, dG0_prime,
                               cid2bounds=cid2bounds, c_range=self.c_range)
        keggpath.normalization = DeltaGNormalization.SIGN_FLUX
        _, mdf = keggpath.FindMDF()
        ln_conc, total_dG_prime = keggpath.GetTotalReactionEnergy(mdf, maximize=True)
        
        odfe = 100 * np.tanh(mdf / (2*R*self.thermo.T))
        concentrations = np.exp(ln_conc)
        
        profile_fig = keggpath.PlotProfile(concentrations)
        plt.title('ODFE = %.1f%%' % odfe, figure=profile_fig)
        self.html_writer.embed_matplotlib_figure(profile_fig, name=prefix+"_prfl")
        keggpath.WriteProfileToHtmlTable(self.html_writer, concentrations)
        
        concentration_fig = keggpath.PlotConcentrations(concentrations)
        plt.title('ODFE = %.1f%%' % odfe, figure=concentration_fig)
        self.html_writer.embed_matplotlib_figure(concentration_fig, name=prefix+"_conc")
        keggpath.WriteConcentrationsToHtmlTable(self.html_writer, concentrations)

        self.write_metabolic_graph(prefix+"_grph", S, rids, cids)
        
        average_dG_prime = total_dG_prime/np.sum(fluxes)
        average_dfe = 100 * np.tanh(-average_dG_prime / (2*R*self.thermo.T))
        
        print ','.join([prefix, '%.1f' % mdf, '%.1f' % -average_dG_prime, 
                        '%.1f' % odfe, '%.1f' % average_dfe, 
                        '%.1f' % total_dG_prime, '%g' % np.sum(fluxes)])
        return "MTDF = %.1f (avg. = %.1f) kJ/mol, ODFE = %.1f%% (avg. = %.1f%%), Total &Delta;<sub>r</sub>G' = %.1f kJ/mol, no. steps = %g" %\
            (mdf, -average_dG_prime, odfe, average_dfe, total_dG_prime, np.sum(fluxes))

if __name__ == '__main__':
    REACTION_FNAME = 'examples/wolf_reactions.txt'
    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    model = KeggModel.from_formulas(reaction_strings)    

    cc = ComponentContribution()
    model.add_thermo(cc)

    fluxes = [1] * len(model.rids)
    mdf = MaxMinDrivingForce(model, fluxes, pH=7.5, I=0.2, T=298.15)
    print mdf.Solve()