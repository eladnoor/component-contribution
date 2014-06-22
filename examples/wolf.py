import sys, logging, os
import numpy as np

REACTION_FNAME = 'examples/wolf_reactions.txt'
CC_CACHE_FNAME = 'cache/component_contribution.mat'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = 'python/component_contribution.py'

def python_main():
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    from python.component_contribution import ComponentContribution
    from python.kegg_reaction import KeggReaction

    reaction_strings = open(REACTION_FNAME, 'r').readlines()

    logging.info('Calculating the component-contributions from raw data')
    cc = ComponentContribution()

    reactions = map(KeggReaction.parse_formula, reaction_strings)
    ddG0 = np.array([r.get_transform_ddG0(pH=7.5, I=0.2, T=298.15)
                     for r in reactions])
    
    res = np.array([cc.get_dG0_r(r) for r in reactions])
    
    dG0_prime = res[:, 0] + ddG0
    std = res[:, 1]

    for i, r in enumerate(reactions):
        print "dG0' = %8.1f +- %5.1f" % (dG0_prime[i], std[i] * 1.96)

if __name__ == '__main__':
    pwd = os.path.realpath(os.path.curdir)
    _, directory = os.path.split(pwd)
    if directory == 'examples':
        sys.path.append(os.path.abspath('..'))
        os.chdir('..')
        
    python_main()
