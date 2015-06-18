# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 09:53:00 2015

@author: noore
"""
from component_contribution.compound_cacher import CompoundCacher
from component_contribution.kegg_reaction import KeggReaction
from component_contribution import inchi2gv
from scipy.io import savemat, loadmat
import numpy as np
import os
import argparse

def decompose_reaction(ccache, decomposer, cids, G, reaction):
    """
        calculate the reaction stoichiometric vector and the group incidence
        vector (x and g)
    """
    Nc, Ng = G.shape
    
    x = np.matrix(np.zeros((Nc, 1)))
    x_prime = []
    G_prime = []

    for compound_id, coeff in reaction.iteritems():
        if compound_id in cids:
            i = cids.index(compound_id)
            x[i, 0] = coeff
        else:
            # Decompose the compound and calculate the 'formation energy'
            # using the group contributions.
            # Note that the length of the group contribution vector we get 
            # from CC is longer than the number of groups in "groups_data" 
            # since we artifically added fictive groups to represent all the 
            # non-decomposable compounds. Therefore, we truncate the 
            # dG0_gc vector since here we only use GC for compounds which
            # are not in cids_joined anyway.
            x_prime.append(coeff)
            comp = ccache.get_compound(compound_id)
            group_vec = decomposer.smiles_to_groupvec(comp.smiles_pH7)
            G_prime.append(group_vec.ToArray())

    if x_prime != []:
        g = np.matrix(x_prime) * np.vstack(G_prime)
    else:
        g = np.matrix(np.zeros((1, 1)))

    g.resize((Ng, 1))

    return x, g

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        'Parse a file with a list of reactions in KEGG format and prepare '
        'the necessary matrices for running CC.')
    parser.add_argument('--train_file', type=argparse.FileType('rb'),
                       help='path to the .mat file containing the CC parameters')
    parser.add_argument('--rxn_file', type=argparse.FileType('r'),
                       help='path to a text file containing the reactions')
    parser.add_argument('--out_file', type=argparse.FileType('wb'),
                       help='path to the output .mat file')
    
    args = parser.parse_args()

    ccache = CompoundCacher()
    groups_data = inchi2gv.init_groups_data()
    decomposer = inchi2gv.InChIDecomposer(groups_data)
    w, b, G, cids, S = map(loadmat(args.train_file).get, ['w', 'b', 'G', 'cids', 'S'])
    cids = list(cids.flat)
    Nc, Ng = G.shape
    
    model_X = []
    model_G = []
    for line in args.rxn_file.readlines():
        reaction = KeggReaction.parse_formula(line)
        try:
            x, g = decompose_reaction(ccache, decomposer, cids, G, reaction)
        except inchi2gv.GroupDecompositionError:
            x = np.zeros((Nc, 1))
            g = np.zeros((Ng, 1))
        model_X.append(list(x.flat))
        model_G.append(list(g.flat))

    mdict = {
                'X' : np.array(model_X).T,
                'G': np.array(model_G).T
            }
    savemat(args.out_file, mdict, appendmat=False, format='5', 
            long_field_names=False, do_compression=True, oned_as='row')