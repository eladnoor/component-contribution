# -*- coding: utf-8 -*-
"""
Created on Wed May 28 13:13:37 2014

@author: eladn
"""

import logging, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from component_contribution.thermodynamic_constants import R, default_c_range, default_c_mid
from component_contribution.component_contribution import ComponentContribution
from component_contribution.kegg_model import KeggModel

from scripts.mdf_dual import KeggPathway
from scripts.html_writer import HtmlWriter, NullHtmlWriter
from scripts.kegg_parser import ParsedKeggFile
    
class ConcentrationConstraints(object):
    pass

class MaxMinDrivingForce(object):
    
    def __init__(self, model, fluxes, bounds, pH, I, T, html_writer=None,
                 cid2name=None):
        """
            model    - a KeggModel object
            fluxes   - a list of fluxes, should match the model in length
            bounds   - a dictionary mapping KEGG compound IDs to tuples of 
                       (low,up) bound
            pH, I, T - the aqueous conditions for the thermodynamics
            html_writer - (optional) write a report to HTML file
        """
        self.model = model
        self.fluxes = fluxes
        self.cid2name = cid2name
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

    def MapCIDs(self, cid_dict):
        try:
            self.cids = map(cid_dict.get, self.cids)
        except KeyError as e:
            raise KeyError('Error: not all CIDs are mapped to their new names '
                           + str(e))
    def GetBounds(self):
        cid2bounds = {cid:self.c_range for cid in self.model.cids}
        cid2bounds['C00001'] = (1, 1) # the default for H2O is 1
        for cid, b in self.bounds.iteritems():
            cid2bounds[cid] = b
        return cid2bounds

    def Solve(self, uncertainty_factor=3.0, diagonal_covariance=False):
        S = self.model.S
        f = self.fluxes
        cids = self.model.cids
        rids = self.model.rids or ['R%05d' % i for i in xrange(S.shape[1])]
        dG0_prime, dG0_std, sqrt_Sigma = self.model.get_transformed_dG0(pH=self.pH, I=self.I, T=self.T)
        if diagonal_covariance:
            sigma = uncertainty_factor*np.matrix(np.diag(dG0_std))
        else:
            sigma = uncertainty_factor*sqrt_Sigma

        cid2bounds = self.GetBounds()

        keggpath = KeggPathway(S, rids, f, cids, dG0_prime, sigma,
                               cid2bounds=cid2bounds,
                               c_range=self.c_range,
                               cid2name=self.cid2name)
        _mdf, params = keggpath.FindMDF(calculate_totals=True)
        total_dG_prime = params.get('maximum total dG', np.nan)
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

        return _mdf, params['gibbs energies raw'], dG0_std

    def SolveIterative(self, uncertainty_factor=3.0, diagonal_covariance=False):
        S = self.model.S
        f = self.fluxes
        rids = self.model.rids or ['R%05d' % i for i in xrange(S.shape[1])]
        cids = self.model.cids
        dG0_prime, dG0_std, sqrt_Sigma = self.model.get_transformed_dG0(pH=self.pH, I=self.I, T=self.T)
        if diagonal_covariance:
            #sigma = uncertainty_factor * np.matrix(np.diag(dG0_std))
            sigma = 0.0 * sqrt_Sigma
            dG0_prime -= uncertainty_factor * dG0_std
        else:
            sigma = uncertainty_factor * sqrt_Sigma
        
        cid2bounds = self.GetBounds()
        
        rid2bounds = {}
        total_active_reactions = len(filter(None, f.flat))
        
        iter_counters = [-1] * len(rids)
        params_list = []
        for i in xrange(len(rids)):
            keggpath = KeggPathway(S, rids, f, cids, dG0_prime, sigma,
                                   rid2bounds=rid2bounds,
                                   cid2bounds=cid2bounds,
                                   c_range=self.c_range,
                                   cid2name=self.cid2name)
            _mdf, params = keggpath.FindMDF(calculate_totals=False)
            params_list.append(params)
            
            # is not the same and maybe there is a mixup
            tmp = zip(rids,
                      map(lambda x:'%.1f' % x, params['gibbs energies'].flat),
                      map(lambda x:'%.1f' % x, params['reaction prices'].flat))
            
            logging.debug('\n'.join(map(', '.join, tmp)))
        
            # fix the driving force of the reactions that have shadow prices
            # to the MDF value, and remove them from the optimization in the
            # next round
            shadow_prices = params['reaction prices']
            
            print '\rIterative MDF: %3d%%' % \
                (len(rid2bounds) * 100 / total_active_reactions),
            for rid, s_p in zip(rids, shadow_prices):
                if rid not in rid2bounds and s_p > 1e-5:
                    rid2bounds[rid] = -_mdf + 1e-4 # add 'epsilon' for numerical reasons
                    iter_counters[rids.index(rid)] = i

            if len(rid2bounds) == total_active_reactions: 
                break
        
        print '\rIterative MDF: [DONE]'
        self.html_writer.write("<p>MDF = %.1f kJ/mol</p>\n" % params_list[0]['MDF'])
        
        params_list[-1]['profile figure'] = keggpath.PlotProfile(params_list[-1])
        self.html_writer.embed_matplotlib_figure(params_list[-1]['profile figure'],
                                                 width=320, height=320)
        self.WriteIterativeReport(keggpath, params_list, iter_counters)
        return params_list[-1]

    def WriteIterativeReport(self, keggpath, params_list, iter_counters):
        headers = ["reaction", "formula", "flux",
                   "&Delta;<sub>r</sub>G' [kJ/mol]"] + ['I%02d' % i for i in xrange(len(params_list))]
        dict_list = []
        for r, rid in enumerate(keggpath.rids):
            if keggpath.fluxes[0, r] == 0:
                continue
            d = {'reaction'  : rid,
                 'flux'      : keggpath.fluxes[0, r],
                 'formula'   : keggpath.GetReactionString(r),
                 headers[3]  : params_list[-1]['gibbs energies'][r, 0],
                 'iteration' : iter_counters[r]}
            for i, p in enumerate(params_list):
                if i < iter_counters[r]:
                    d['I%02d' % i] = '%.3f' % p['gibbs energies'][r, 0]
                else:
                    d['I%02d' % i] = '<b>%.3f</b>' % p['gibbs energies'][r, 0]
            
            dict_list.append(d)

        dict_list.sort(key=lambda d:d['iteration'], reverse=False)

        d = {'reaction' : 'Total',
             'flux'     : '1',
             'formula'  : keggpath.GetTotalReactionString(),
             headers[3] : float(keggpath.fluxes * params_list[-1]['gibbs energies'])}
        dict_list.append(d)
        
        self.html_writer.write_table(dict_list, headers=headers, decimal=1)

        concentrations = params_list[-1]['concentrations']

        headers = ['compound', 'Concentration LB [M]',
                   'Concentration [M]', 'Concentration UB [M]']
        dict_list = []
        for c, cid in enumerate(keggpath.cids):
            d = {}
            d['compound'] = keggpath.c_names[c]
            lb, ub = keggpath.GetConcentrationBounds(cid)
            d['Concentration LB [M]'] = '%.2e' % lb
            d['Concentration [M]'] = '%.2e' % concentrations[c, 0]
            d['Concentration UB [M]'] = '%.2e' % ub
            dict_list.append(d)
       
        self.html_writer.write_table(dict_list, headers=headers)

###############################################################################
            
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
    #fname = sys.argv[1]
    fname = 'mdf_pathways'
    
    REACTION_FNAME = 'scripts/%s.txt' % fname
    pathways = KeggFile2ModelList(REACTION_FNAME)
    html_writer = HtmlWriter('res/%s.html' % fname)
    cc = ComponentContribution.init()

    matdict = {}
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
        
        matdict[p['entry'] + '.dGm_prime'] = dGm_prime
        matdict[p['entry'] + '.dG0_std'] = dG0_std
        
    savemat('res/%s.mat' % fname, matdict)
