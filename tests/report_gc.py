# -*- coding: utf-8 -*-
import sys
import numpy as np
import logging

logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
sys.path.append('../python')
from training_data import TrainingData
from component_contribution import ComponentContribution
from kegg_model import KeggModel
import inchi2gv

REACTION_FNAME = 'report_gc_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = '../python/component_contribution.py'

reaction_strings = open(REACTION_FNAME, 'r').readlines()
model = KeggModel.from_formulas(reaction_strings)    

td = TrainingData()
cc = ComponentContribution(td)
model.add_thermo(cc)

dG0_prime, dG0_std = model.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)

G = cc.params['G']

groups_data = inchi2gv.init_groups_data()
group_names = groups_data.GetGroupNames()
for i, group_name in enumerate(group_names):
    print '%s - %d examples' % (group_name, np.sum(G[:, i] != 0))
