import sys, logging, os

REACTION_FNAME = 'examples/wolf_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = 'python/component_contribution.py'

def python_main():
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    from python.component_contribution import ComponentContribution
    from python.kegg_model import KeggModel
    from python.kegg_reaction import KeggReaction

    reaction_strings = open(REACTION_FNAME, 'r').readlines()

    if True:
        cc = ComponentContribution()
        for s in reaction_strings:
            reaction = KeggReaction.parse_formula(s)
            ddG0 = reaction.get_transform_ddG0(pH=7.5, I=0.2, T=298.15)
            dG0_cc, s_cc = cc.get_dG0_r(reaction)
            dG0_prime_cc = dG0_cc + ddG0
            print "dG0 = %8.1f +- %5.1f, dG0' = %8.1f +- %5.1f" % \
                (dG0_cc, s_cc * 1.96, dG0_prime_cc, s_cc * 1.96)

    if True: # old piece of code left here for debugging purposes
        cc2 = ComponentContribution()
        model = KeggModel.from_formulas(reaction_strings)    
        model.add_thermo_old(cc2)
        
        dG0_prime, dG0_std = model.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
    
        for i in xrange(dG0_prime.shape[0]):
            dG0_cc = model.dG0[i, 0]
            dG0_prime_cc = dG0_prime[i, 0]
            s_cc = dG0_std[i, 0]
            print "dG0 = %8.1f +- %5.1f, dG0' = %8.1f +- %5.1f" % \
                (dG0_cc, s_cc * 1.96, dG0_prime_cc, s_cc * 1.96)

if __name__ == '__main__':
    pwd = os.path.realpath(os.path.curdir)
    _, directory = os.path.split(pwd)
    if directory == 'examples':
        sys.path.append(os.path.abspath('..'))
        os.chdir('..')
        
    python_main()
