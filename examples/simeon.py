import sys, logging, os
import numpy as np
from scipy.io import savemat

REACTION_FNAME = 'examples/simeon_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = 'python/component_contribution.py'
OUTPUT_FNAME = 'res/simeon.mat'

def python_main():
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    from python.component_contribution import ComponentContribution
    from python.kegg_model import KeggModel

    cc = ComponentContribution.init()
    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    model = KeggModel.from_formulas(reaction_strings)
    model.add_thermo(cc)

    dG0_prime, dG0_std = model.get_transformed_dG0(7.0, 0.2, 298.15)

    print "For a linear problem, define two vector variables 'x' and 'y', each of length Nr (i.e. " + \
          "the same length as the list of reactions). Then add these following " + \
          "constraints: \n" + \
          "-1 <= y <= 1\n" + \
          "x = dG0_prime + 3 * dG0_std * y\n" + \
          "Then use 'x' as the value of the standard Gibbs energy for all reactions."
    print "The results are writteng to: " + OUTPUT_FNAME
    
    savemat(OUTPUT_FNAME,
            {'dG0_prime': dG0_prime, 'dG0_std': dG0_std},
            oned_as='row')
    
if __name__ == '__main__':
    pwd = os.path.realpath(os.path.curdir)
    _, directory = os.path.split(pwd)
    if directory == 'examples':
        sys.path.append(os.path.abspath('..'))
        os.chdir('..')
        
    python_main()
