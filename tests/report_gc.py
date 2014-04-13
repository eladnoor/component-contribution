# -*- coding: utf-8 -*-
import sys
import numpy as np
import logging
from time import time

formatter = logging.Formatter('%(asctime)s %(filename)s %(lineno)s %(levelname)s  %(message)s')

stdout_handler = logging.StreamHandler(sys.stdout)
stderr_handler = logging.StreamHandler(sys.stderr)
file_handler = logging.FileHandler('/tmp/component_contribution_%f.log'%time())
stdout_handler.setFormatter(formatter)
stderr_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)

logger = logging.getLogger('')
logger.addHandler(stdout_handler)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

sys.path.append('../python')
import inchi2gv
from training_data import TrainingData
from component_contribution import ComponentContribution
from kegg_model import KeggModel

reaction_strings = open('../tests/report_gc_reactions.txt', 'r').readlines()
model = KeggModel.from_formulas(reaction_strings)    

td = TrainingData()
cc = ComponentContribution(td)
G = ComponentContribution.create_group_incidence_matrix(cc.train_cids)

groups_data = inchi2gv.init_groups_data()
decomposer = inchi2gv.InChIDecomposer(groups_data)
group_names = groups_data.GetGroupNames()

for i, group_name in enumerate(group_names):
    print '%s - %d examples\n' % (group_name, np.sum(G[:, i] != 0))
