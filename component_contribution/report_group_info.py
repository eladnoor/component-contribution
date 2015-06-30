# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 12:32:49 2015

@author: noore
"""

from component_contribution.component_contribution_trainer import ComponentContribution
import numpy as np
import csv, sys

cc = ComponentContribution.init()
group_names = cc.groups_data.GetGroupNames()
dG0_gc = cc.params['dG0_gc']
G = cc.params['G']

Ng = len(group_names)
group_names += map(cc.cids_joined.__getitem__, np.nonzero(G[:, Ng:])[0])

csv_writer = csv.writer(sys.stdout)
csv_writer.writerow(['Group', 'dG0_gc [kJ/mol]'])
csv_writer.writerows(zip(group_names, map(lambda x : '%.2f' % x, dG0_gc.flat)))
