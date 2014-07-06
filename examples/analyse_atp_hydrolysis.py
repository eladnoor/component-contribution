import logging, os, csv
import numpy as np
logger = logging.getLogger('')
logger.setLevel(logging.INFO)
from python.component_contribution import ComponentContribution
from python.kegg_reaction import KeggReaction
from python.training_data import TrainingData

CC_CACHE_FNAME = 'cache/component_contribution.mat'

def main():

    td = TrainingData()

    if os.path.exists(CC_CACHE_FNAME):
        cc = ComponentContribution.from_matfile(CC_CACHE_FNAME, td)
    else:
        logging.info('Calculating the component-contributions from raw data')
        cc = ComponentContribution(td)
        cc.save_matfile(CC_CACHE_FNAME)

    G1 = np.matrix(cc.params['G1'])
    G2 = np.matrix(cc.params['G2'])
    G3 = np.matrix(cc.params['G3'])

    ############################################################################

    #reaction = KeggReaction.parse_formula('C00002 + C00001 <=> C00008 + C00009'); fname = 'atpase';
    reaction = KeggReaction.parse_formula('C00149 <=> C00122 + C00001'); fname = 'fumarase';

    x, g = cc._decompose_reaction(reaction)
    weights_rc = (x.T * G1).round(5)
    weights_gc = (x.T * G2 + g.T * G3).round(5)
    weights = weights_rc + weights_gc
    orders = sorted(range(weights.shape[1]),
                    key=lambda j:abs(weights[0, j]), reverse=True)

    output = csv.writer(open('res/%s_analysis.csv' % fname, 'w'))
    output.writerow(('Weight', 'dG\'0', 'dG0', 'reference', 'reaction'))
    for j in orders:
        if abs(weights[0, j]) < 1e-7:
            continue
                              
        output.writerow((weights_rc[0, j], td.dG0_prime[j], td.dG0[j],
                         td.reference[j], td.description[j]))

if __name__ == '__main__':
    main()