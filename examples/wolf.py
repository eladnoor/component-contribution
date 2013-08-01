import sys
from subprocess import Popen, PIPE
import numpy as np

REACTION_FNAME = 'wolf_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = '../python/component_contribution.py'

def cmd_main():
    p1 = Popen(['cat', REACTION_FNAME], stdout=PIPE)
    p2 = Popen([PYTHON_BIN, PYTHON_SCRIPT_FNAME], stdin=p1.stdout,
               executable=PYTHON_BIN, stdout=PIPE)
    res = p2.communicate()[0]
    model_dG0 = np.array(res)
    print str(model_dG0)
    
def python_main():
    sys.path.append('../python')
    from training_data import TrainingData
    from component_contribution import ComponentContribution
    from kegg_model import KeggModel

    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    model = KeggModel.from_formulas(reaction_strings)    

    td = TrainingData()
    cc = ComponentContribution(td)
    model.add_thermo(cc)
    
    dG0_prime, dG0_std = model.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
    
    sys.stdout.write('[' + 
                     ', '.join([str(x) for x in model.dG0.flat]) + '; ' + 
                     ', '.join([str(x) for x in dG0_prime.flat]) + 
                     ']')    

if __name__ == '__main__':
    python_main()