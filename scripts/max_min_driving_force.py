# -*- coding: utf-8 -*-
"""
Created on Wed May 28 13:13:37 2014

@author: eladn
"""

import logging, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from python.thermodynamic_constants import R, default_c_range, default_c_mid
from python.component_contribution import ComponentContribution
from python.kegg_model import KeggModel

from scripts.mdf_dual import KeggPathway
from scripts.html_writer import HtmlWriter, NullHtmlWriter
from scripts.kegg_parser import ParsedKeggFile
    
class ConcentrationConstraints(object):
    pass

class MaxMinDrivingForce(object):
    
    def __init__(self, model, fluxes, bounds, pH, I, T, html_writer=None):
        self.model = model
        self.fluxes = np.matrix(fluxes)
        self.pH, self.I, self.T = pH, I, T
        self.c_range = default_c_range
        self.c_mid = default_c_mid
        self.bounds = bounds

        if html_writer is None:
            self.html_writer = NullHtmlWriter()
        else:
            self.html_writer = html_writer
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
        for cid, b in self.bounds.iteritems():
            cid2bounds[cid] = b
        return cid2bounds

    def Solve(self, prefix='mdf'):
        S = self.model.S
        rids = self.model.rids or ['R%05d' % i for i in xrange(S.shape[1])]
        cids = self.model.cids
        dG0_prime, dG0_std = self.model.get_transformed_dG0(pH=self.pH, I=self.I, T=self.T)
        cid2bounds = self.GetBounds()

        keggpath = KeggPathway(S, rids, self.fluxes, cids,
                               formation_energies=None,
                               reaction_energies=dG0_prime.T,
                               cid2bounds=cid2bounds, c_range=self.c_range)
        _mdf, params = keggpath.FindMDF()
        total_dG_prime = params['maximum total dG']
        odfe = 100 * np.tanh(_mdf / (2*R*self.T))
        average_dG_prime = total_dG_prime/np.sum(self.fluxes)
        average_dfe = 100 * np.tanh(-average_dG_prime / (2*R*self.T))        
        res =  ["MDF = %.1f (avg. = %.1f) kJ/mol" % (_mdf, -average_dG_prime),
               "ODFE = %.1f%% (avg. = %.1f%%)" % (odfe, average_dfe),
               "Total &Delta;<sub>r</sub>G' = %.1f kJ/mol" % total_dG_prime,
               "no. steps = %g" % np.sum(self.fluxes)]
        self.html_writer.write_ul(res)
        
        profile_fig = keggpath.PlotProfile(params)
        plt.title('ODFE = %.1f%%' % odfe, figure=profile_fig)
        self.html_writer.embed_matplotlib_figure(profile_fig, width=320, height=320)
        keggpath.WriteProfileToHtmlTable(self.html_writer, params)
        keggpath.WriteConcentrationsToHtmlTable(self.html_writer, params)

        concentrations = keggpath.GetMillimolarConcentrations()
        dGm_prime = keggpath.CalculateReactionEnergiesUsingConcentrations(concentrations)
        
        return _mdf, dGm_prime, dG0_std

def KeggFile2ModelList(pathway_file):
    kegg_file = ParsedKeggFile.FromKeggFile(pathway_file)
    entries = kegg_file.entries()
    pathways = []
    for entry in entries:
        fields = kegg_file[entry]
        rids, fluxes, reactions = ParsedKeggFile.ParseReactionModule(fields)
        bounds = ParsedKeggFile.ParseBoundModule(fields)
        model = KeggModel.from_formulas(reactions)
        model.rids = rids
        pH = fields.GetFloatField('PH', 7.5)
        I = fields.GetFloatField('I', 0.2)
        T = fields.GetFloatField('T', 298.15)
        pathways.append({'entry': entry, 'model': model, 'fluxes': fluxes,
                         'bounds': bounds, 'pH': pH, 'I': I, 'T': T})
    return pathways

if __name__ == '__main__':
    fname = sys.argv[1]
    
    REACTION_FNAME = 'scripts/%s.txt' % fname
    pathways = KeggFile2ModelList(REACTION_FNAME)
    html_writer = HtmlWriter('res/%s.html' % fname)

    cc = ComponentContribution()
    cc.train()
    for p in pathways:
        html_writer.write('<h2>%s</h2>' % p['entry'])

        p['model'].add_thermo(cc)

        dG0, u = p['model'].get_transformed_dG0(pH=p['pH'], I=p['I'], T=p['T'])
        mdf = MaxMinDrivingForce(p['model'], p['fluxes'], p['bounds'],
                                 pH=p['pH'], I=p['I'], T=p['T'],
                                 html_writer=html_writer)


        #print 'dG\'0 = ' + ', '.join(map(lambda x:'%.2f' % x, dG0[:,0]))
        mdf_solution, dGm_prime, dG0_std = mdf.Solve()
        logging.info('Pathway %s: MDF = %.1f' % (p['entry'], mdf_solution))
        savemat('res/%s.mat' % fname, {'dGm_prime':dGm_prime, 'dG0_std':dG0_std})
