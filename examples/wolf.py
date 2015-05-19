import sys, logging, os
import numpy as np

REACTION_FNAME = 'examples/wolf_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = 'python/component_contribution.py'

def python_main():
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    import component_contribution
    from component_contribution.thermodynamic_constants import default_RT

    cc = component_contribution.component_contribution.ComponentContribution.init()
    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    model = component_contribution.kegg_model.KeggModel.from_formulas(reaction_strings)

    model.add_thermo(cc)
    dG0_prime, dG0_std, sqrt_Sigma = model.get_transformed_dG0(7.0, 0.1, 298.15)
    
    mM_conc = 1e-3 * np.matrix(np.ones((len(model.cids), 1)))
    if 'C00001' in model.cids:
        mM_conc[model.cids.index('C00001'), 0] = 1.0
    dGm_prime = dG0_prime + default_RT * model.S.T * np.log(mM_conc)
    
    for i, r in enumerate(reaction_strings):
        print '-'*50
        print r.strip()
        print "dG0 = %8.1f +- %5.1f" % (dG0_prime[i, 0], dG0_std[i, 0] * 1.96)
        print "dGm = %8.1f +- %5.1f" % (dGm_prime[i, 0], dG0_std[i, 0] * 1.96)

if __name__ == '__main__':
    pwd = os.path.realpath(os.path.curdir)
    _, directory = os.path.split(pwd)
    if directory == 'examples':
        sys.path.append(os.path.abspath('..'))
        os.chdir('..')
        
    python_main()
