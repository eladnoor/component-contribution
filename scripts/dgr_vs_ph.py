
import sys
import numpy as np
from scipy.io import savemat
from python.component_contribution import ComponentContribution
from python.kegg_model import KeggModel

from scripts.kegg_parser import ParsedKeggFile

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

    pH_range = np.arange(5, 9.0001, 0.1)
    I = 0.2
    T = 298.15

    cc = ComponentContribution.init()

    cc_dict = {'pH' : np.matrix(pH_range).T}
    for p in pathways:
        entry = p['entry']
        p['model'].add_thermo(cc)
        
        dG0_primes = []
        dG0_std = []
        for pH in pH_range:
            dG0_prime, dG0_std = p['model'].get_transformed_dG0(pH=pH, I=I, T=T)
            dG0_primes.append(dG0_prime)
            
        dG0_std = np.matrix(dG0_std)
        cc_dict['%s.dG0_primes' % entry] = dG0_primes
        cc_dict['%s.dG0_std' % entry] = dG0_std.T
        
    savemat('res/%s.mat' % fname, cc_dict)
