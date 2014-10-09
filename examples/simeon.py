import sys, logging, os
import numpy as np

REACTION_FNAME = 'examples/simeon_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = 'python/component_contribution.py'

def python_main():
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    from python.component_contribution import ComponentContribution
    from python.kegg_model import KeggModel

    cc = ComponentContribution.init()
    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    Nr = len(reaction_strings)
    model = KeggModel.from_formulas(reaction_strings)
    model.add_thermo(cc)

    dG0_prime, dG0_std = model.get_transformed_dG0(7.0, 0.2, 298.15)

    y = np.matrix(np.random.multivariate_normal(np.zeros(Nr), np.eye(Nr))).T
    dG0 = dG0_prime + dG0_std * y
    print "A random sample of the standard Gibbs energies: " +\
          ", ".join(['%.1f' % x for x in dG0.flat])
    
    print "For a linear problem, define two vector variables 'x' and 'y', each of length Nr (i.e. " + \
          "the same length as the list of reactions). Then add these following " + \
          "constraints: \n" + \
          "-1 <= y <= 1\n" + \
          "x = dG0_prime + 3 * dG0_std * y\n" + \
          "Then use 'x' as the value of the standard Gibbs energy for all reactions."    
    
if __name__ == '__main__':
    pwd = os.path.realpath(os.path.curdir)
    _, directory = os.path.split(pwd)
    if directory == 'examples':
        sys.path.append(os.path.abspath('..'))
        os.chdir('..')
        
    python_main()
