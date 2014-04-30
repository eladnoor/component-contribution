import sys, logging, os

REACTION_FNAME = 'examples/wolf_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = 'python/component_contribution.py'


def python_main():
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    from python.training_data import TrainingData
    from python.component_contribution import ComponentContribution
    from python.kegg_model import KeggModel

    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    model = KeggModel.from_formulas(reaction_strings)    

    td = TrainingData()
    cc = ComponentContribution(td)
    model.add_thermo(cc)
    
    dG0_prime, dG0_std = model.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
    
    sys.stdout.write('[' + 
                     ', '.join(['%.1f' % x for x in model.dG0.flat]) + '; ' + 
                     ', '.join(['%.1f' % x for x in dG0_prime.flat]) + 
                     ']' + '\n')   
                     
    return cc.params

if __name__ == '__main__':
    pwd = os.path.realpath(os.path.curdir)
    _, directory = os.path.split(pwd)
    if directory == 'examples':
        sys.path.append(os.path.abspath('..'))
        os.chdir('..')
        
    params = python_main()
