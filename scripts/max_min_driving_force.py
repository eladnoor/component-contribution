# -*- coding: utf-8 -*-
"""
Created on Wed May 28 13:13:37 2014

@author: eladn
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from python.thermodynamic_constants import R, default_c_range, default_c_mid
from python.component_contribution import ComponentContribution
from python.kegg_model import KeggModel

from scripts.mdf_dual import KeggPathway
from scripts.html_writer import HtmlWriter
    
class ConcentrationConstraints(object):
    pass

class MaxMinDrivingForce(object):
    
    def __init__(self, model, fluxes, pH, I, T):
        self.model = model
        self.fluxes = np.matrix(fluxes)
        self.pH, self.I, self.T = pH, I, T
        self.c_range = default_c_range
        self.c_mid = default_c_mid

        self.html_writer = HtmlWriter('res/max_min_driving_force.html')
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
        return cid2bounds

    def Solve(self, prefix='mdf'):
        S = self.model.S
        rids = self.model.rids or ['R%05d' % i for i in xrange(S.shape[1])]
        cids = self.model.cids
        dG0_prime, dG0_std = model.get_transformed_dG0(pH=self.pH, I=self.I, T=self.T)
        cid2bounds = self.GetBounds()

        keggpath = KeggPathway(S, rids, self.fluxes, cids,
                               formation_energies=None,
                               reaction_energies=dG0_prime.T,
                               cid2bounds=cid2bounds, c_range=self.c_range)
        _mdf, params = keggpath.FindMDF()
        concentrations = params['concentrations']
        total_dG_prime = params['maximum total dG']
        
        odfe = 100 * np.tanh(_mdf / (2*R*self.T))
        
        profile_fig = keggpath.PlotProfile(concentrations)
        plt.title('ODFE = %.1f%%' % odfe, figure=profile_fig)
        #self.html_writer.embed_matplotlib_figure(profile_fig, name=prefix+"_prfl")
        keggpath.WriteProfileToHtmlTable(self.html_writer, concentrations)
        
        #concentration_fig = keggpath.PlotConcentrations(concentrations)
        #plt.title('ODFE = %.1f%%' % odfe, figure=concentration_fig)
        #self.html_writer.embed_matplotlib_figure(concentration_fig, name=prefix+"_conc")
        keggpath.WriteConcentrationsToHtmlTable(self.html_writer, concentrations)

        #self.write_metabolic_graph(prefix+"_grph", S, rids, cids)
        
        average_dG_prime = total_dG_prime/np.sum(fluxes)
        average_dfe = 100 * np.tanh(-average_dG_prime / (2*R*self.T))
        
        print ','.join([prefix, '%.1f' % _mdf, '%.1f' % -average_dG_prime, 
                        '%.1f' % odfe, '%.1f' % average_dfe, 
                        '%.1f' % total_dG_prime, '%g' % np.sum(fluxes)])
        return "MDF = %.1f (avg. = %.1f) kJ/mol, " % (_mdf, -average_dG_prime) + \
               "ODFE = %.1f%% (avg. = %.1f%%), " % (odfe, average_dfe) + \
               "Total &Delta;<sub>r</sub>G' = %.1f kJ/mol, " % total_dG_prime + \
               "no. steps = %g" % np.sum(fluxes)

if __name__ == '__main__':
    REACTION_FNAME = 'examples/wolf_reactions.txt'
    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    model = KeggModel.from_formulas(reaction_strings)
    if model is None:
        sys.exit(-1)

    cc = ComponentContribution()
    model.add_thermo(cc)

    fluxes = [1] * model.S.shape[1]
    mdf = MaxMinDrivingForce(model, fluxes, pH=7.5, I=0.2, T=298.15)
    print mdf.Solve()