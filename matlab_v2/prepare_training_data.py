# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:25:49 2015

@author: noore
"""
import sys
sys.path.append('..')

from component_contribution.component_contribution_trainer import ComponentContribution
from scipy.io import savemat
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        'Prepare all thermodynamic training data in a .mat file for running '
        'component contribution.')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                       help='the path to the .mat file that should be written '
                       'containing the training data')
    
    args = parser.parse_args()
    cc = ComponentContribution()

    d = {'w': cc.train_w,
         'b': cc.train_b,
         'G': cc.create_group_incidence_matrix(),
         'cids': cc.train_cids,
         'S': cc.train_S}
    savemat(args.outfile, d, oned_as='row')
