import sys, logging, os
import numpy as np

REACTION_FNAME = 'examples/wolf_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = 'python/component_contribution.py'

def python_main():
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    from python.component_contribution import ComponentContribution
    from python.kegg_model import KeggModel
    from python.thermodynamic_constants import default_RT

    cc = ComponentContribution.init()
    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    model = KeggModel.from_formulas(reaction_strings)

    model.add_thermo(cc)
    dG0_prime, dG0_std = model.get_transformed_dG0(7.0, 0.2, 298.15)
    
    cid2c_range = {cid : (1e-6, 1e-2) for cid in model.cids}
    cid2c_range['C00001'] = (1.0, 1.0) # water must be always set to 1 (equivalent to 55M)
    
    ln_conc = np.log(np.matrix([cid2c_range[cid] for cid in model.cids]))
    ln_conc_min = ln_conc[:, 0]
    ln_conc_max = ln_conc[:, 1]
    
    S_minus, S_plus = model.get_unidirectional_S()
    
    dG_prime_min = dG0_prime + default_RT * (S_minus.T * ln_conc_max + S_plus.T * ln_conc_min)
    dG_prime_max = dG0_prime + default_RT * (S_minus.T * ln_conc_min + S_plus.T * ln_conc_max)
    
    for i, r in enumerate(reaction_strings):
        print '-'*50
        print r.strip()
        print "dG'0      = %8.1f +- %5.1f" % (dG0_prime[i, 0], dG0_std[i, i] * 1.96)
        print "dG' range = %8.1f  - %8.1f" % (dG_prime_min[i, 0], dG_prime_max[i, 0])
        if dG_prime_min[i, 0] < 0 and dG_prime_max[i, 0] > 0:
            print "REVERSIBLE!"
        else:
            print "IRREVERSIBLE!"

if __name__ == '__main__':
    pwd = os.path.realpath(os.path.curdir)
    _, directory = os.path.split(pwd)
    if directory == 'examples':
        sys.path.append(os.path.abspath('..'))
        os.chdir('..')
        
    python_main()
