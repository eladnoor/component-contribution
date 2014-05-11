import sys, csv, os
import numpy as np
import matplotlib.pyplot as plt
from python.component_contribution import ComponentContribution
from python.kegg_model import KeggModel

cid2dG0 = {}
if not os.path.exists(sys.argv[1]):
    fp = open(sys.argv[1], 'w')
    cc = ComponentContribution()
    cc.train()
    csv_out = csv.writer(fp)
    csv_out.writerow(['cid', 'dG0_f'])
    for compound_id in cc.ccache.get_all_compound_ids():
        dG0_f = cc.get_major_ms_dG0_f(compound_id)
        csv_out.writerow([compound_id, '%8.2f' % dG0_f])
        cid2dG0[compound_id] = dG0_f
else:
    for row in csv.DictReader(open(sys.argv[1], 'r')):
        cid2dG0[row['cid']] = float(row['dG0_f'])

REACTION_FNAME = 'tests/report_gc_reactions.txt'
reaction_strings = open(REACTION_FNAME, 'r').readlines()
model = KeggModel.from_formulas(reaction_strings)    

# compare the dG0_r of the model to the one we get if multiplying
# the model stoichiometric matrix by the formation energies
model_dG0_f = np.matrix([cid2dG0[cid] for cid in model.cids]).T
model_dG0_r = model.S.T * model_dG0_f

cc2 = ComponentContribution()
model.add_thermo(cc2)
plt.plot(model.dG0, model_dG0_r, '.')
plt.show()

#dG0_prime, dG0_std = model.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)

