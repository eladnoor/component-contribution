import sys, logging, os

REACTION_FNAME = 'examples/wolf_reactions.txt'
CC_CACHE_FNAME = 'cache/component_contribution.mat'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = 'python/component_contribution.py'

def python_main():
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    from python.component_contribution import ComponentContribution
    from python.kegg_reaction import KeggReaction

    if os.path.exists(CC_CACHE_FNAME):
        cc = ComponentContribution.from_matfile(CC_CACHE_FNAME)
    else:
        logging.info('Calculating the component-contributions from raw data')
        cc = ComponentContribution()
        cc.save_matfile(CC_CACHE_FNAME)

    reaction_strings = open(REACTION_FNAME, 'r').readlines()
    reactions = map(KeggReaction.parse_formula, reaction_strings)

    for i, r in enumerate(reactions):
        dG0, u_r, analysis = cc.get_dG0_r(r, True)
        ddG0 = r.get_transform_ddG0(pH=7, I=0.25, T=298.15)
        dG0_prime = dG0 + ddG0
        print '-'*50
        print r
        print "dG0' = %8.1f +- %5.1f" % (dG0_prime, u_r * 1.96)
        #for d in analysis:
        #    print 'Wrc = %3g, Wgc = %3g, N = %d : %s' % \
        #        (d['w_rc'], d['w_gc'], d['count'], d['reaction'])

if __name__ == '__main__':
    pwd = os.path.realpath(os.path.curdir)
    _, directory = os.path.split(pwd)
    if directory == 'examples':
        sys.path.append(os.path.abspath('..'))
        os.chdir('..')
        
    python_main()
