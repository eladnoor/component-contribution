#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 15:47:09 2017

@author: noore
"""

import logging
from numpy import matrix, array, load, zeros, logaddexp, sqrt, log
import os
import json
import gzip
import types
import re

COMPOUND_JSON_FNAME = 'cc_compounds.json.gz'
PREPROCESS_FNAME = 'cc_preprocess.npz'

R = 8.31e-3   # kJ/(K*mol)
DEFAULT_TEMP = 298.15  # K
DEFAULT_IONIC_STRENGTH = 0.1  # mM
DEFAULT_PH = 7.0
DEFAULT_PMG = 14.0
DEFAULT_PHASE = 'aqueous'
RT = R * DEFAULT_TEMP
RTlog10 = RT * log(10)

DEBYE_HUECKLE_A = 2.91482
DEBYE_HUECKLE_B = 1.6
MG_FORMATION_ENERGY = -455.3 # kJ/mol, formation energy of Mg2+

class Species(object):
    
    def __init__(self, d):
        self.dG0_f = d['dG0_f']
        self.phase = d['phase']
        self.nH = d['nH']
        self.nMg = d['nMg']
        self.z = d['z']
            
    def ddG_prime(self, pH, pMg, I):
        """
            Transform this individual estimate to difference conditions.
        """

        sqrt_I = sqrt(I)
        ddG_prime = 0

        # add the potential related to the pH
        if self.nH > 0:
            ddG_prime += self.nH * RTlog10 * pH

        # add the potential related to the ionic strength
        ddG_prime -= DEBYE_HUECKLE_A * (self.z ** 2 - self.nH) * sqrt_I / (1.0 + DEBYE_HUECKLE_B * sqrt_I)

        # add the potential related to the Mg ions
        if self.nMg > 0:
            ddG_prime += self.nMg * (RTlog10 * pMg - MG_FORMATION_ENERGY)

        return ddG_prime
        
    def dG0_prime(self, pH, pMg, I):
        """
            Transform this individual estimate to difference conditions.
        """
        dG0_f_prime = self.dG0_f + self.ddG_prime(pH, pMg, I)
        logging.info('nH = %2d, nMg = %2d, z = %2d, dG0_f = %6.1f -> dG\'0_f = %6.1f' %
                     (self.nH, self.nMg, self.z, self.dG0_f, dG0_f_prime))
        return dG0_f_prime
            
    
class Compound(object):
    
    def __init__(self, d):
        self.inchi = d.get('InChI', '')
        self.kegg_id = d['CID']
        self.compound_index = d.get('compound_index', -1)
        self.group_vector = d.get('group_vector', None)
        self.formula = d.get('formula', '')
        self.mass = d.get('mass', -1)
        self.num_electrons = d.get('num_electrons', 0)
        self.no_dg_explanation = d.get('error', '')
        if 'pmap' in d:
            self.source = d['pmap'].get('source', '')
            self.species_list = map(Species, d['pmap'].get('species', []))
            
        self.phase = DEFAULT_PHASE
        for sp in self.species_list:
            if sp.phase != DEFAULT_PHASE:
                if len(self.species_list) > 1:
                    raise ValueError('compound in non-aqueous phase must '
                                     'have only one species')
                self.phase = sp.phase
            
    def get_stoich_vector(self, Nc):
        x = matrix(zeros((Nc, 1)))
        i = self.compound_index
        if i is None:
            raise Exception('could not find index for ' + self.kegg_id)
        x[i, 0] = 1
        return x

    def get_group_incidence_vector(self, Ng):
        g = matrix(zeros((Ng, 1)))
        gv = self.group_vector
        if gv is None:
            raise Exception('could not find group vector for ' + self.kegg_id)
        for g_ind, g_count in gv:
            g[g_ind, 0] += g_count
        return g

    def dG0_prime(self, pH, pMg, I):
        """
            Get a detla-deltaG estimate for this group of species.
            I.e., this is the difference between the dG0 and the dG'0, which
            only depends on the pKa of the pseudoisomers, but not on their
            formation energies.

            Args:
                pH  - the pH to estimate at.
                pMg - the pMg to estimate at.
                I   - the ionic strength to estimate at.

            Returns:
                The estimated delta G in the given conditions or None.
        """
        if self.phase == DEFAULT_PHASE:
            # Compute per-species transforms, scaled down by R*T.
            dG0_prime_vec = array(map(lambda s: s.dG0_prime(pH, pMg, I),
                                      self.species_list))
            
            # Numerical issues: taking a sum of exp(v) for |v| quite large.
            # Use the fact that we take a log later to offset all values by a
            # constant (the minimum value).
            dG0_f_prime = -RT * logaddexp.reduce((-1.0 / RT) * dG0_prime_vec)
        else:
            dG0_f_prime = self.species_list[0].dG0_prime(pH, pMg, I)
        
        logging.info('KEGG_ID = %s, dG\'0_f = %.1f' % (self.kegg_id, dG0_f_prime))
        return dG0_f_prime
    
class Reaction(object):
    # load formation energies from the JSON file
    COMPOUND_DICT = {}
    for cd in json.load(gzip.open(COMPOUND_JSON_FNAME, 'r')):
        kegg_id = cd.get('CID', 'unknown')
        COMPOUND_DICT[kegg_id] = cd

    def __init__(self, kegg_id_to_coeff):
        self.kegg_id_to_coeff = kegg_id_to_coeff
        
        # Create the relevant "Compound" objects and store in a dictionary
        self.kegg_id_to_compound = {}
        for kegg_id in self.kegg_id_to_coeff.keys():
            compound = Compound(Reaction.COMPOUND_DICT[kegg_id])
            self.kegg_id_to_compound[kegg_id] = compound
    
    def kegg_ids(self):
        return self.kegg_id_to_coeff.keys()

    def get_coeff(self, kegg_id):
        return self.kegg_id_to_coeff.get(kegg_id, 0)

    def get_compound(self, kegg_id):
        return self.kegg_id_to_compound.get(kegg_id, None)
    
    def dG0_prime(self, pH, pMg, I):
        dG0_r_prime = 0
        for kegg_id in self.kegg_ids():
            coeff = self.get_coeff(kegg_id)
            compound = self.get_compound(kegg_id)
            dG0_r_prime += coeff * compound.dG0_prime(pH, pMg, I)
        return dG0_r_prime
    
    @staticmethod
    def parse_reaction_formula_side(s):
        """ 
            Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
            Ignores stoichiometry.
            
            Returns:
                The set of CIDs.
        """
        if s.strip() == "null":
            return {}
        
        compound_bag = {}
        for member in re.split('\s+\+\s+', s):
            tokens = member.split(None, 1)
            if len(tokens) == 0:
                continue
            if len(tokens) == 1:
                amount = 1
                key = member
            else:
                amount = float(tokens[0])
                key = tokens[1]
            compound_bag[key] = compound_bag.get(key, 0) + amount
        
        return compound_bag
    
    @staticmethod
    def parse_formula(formula, arrow='='):
        """ 
            Parses a two-sided formula such as: 2 C00001 => C00002 + C00003 
            
            Return:
                The set of substrates, products and the direction of the reaction
        """
        tokens = formula.split(arrow)
        if len(tokens) < 2:
            raise ValueError('Reaction does not contain the arrow sign (%s): %s'
                             % (arrow, formula))
        if len(tokens) > 2:
            raise ValueError('Reaction contains more than one arrow sign (%s): %s'
                             % (arrow, formula))
        
        left = tokens[0].strip()
        right = tokens[1].strip()
        
        sparse_reaction = {}
        for cid, count in Reaction.parse_reaction_formula_side(left).iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 
    
        for cid, count in Reaction.parse_reaction_formula_side(right).iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 
        
        return Reaction(sparse_reaction)
    
    
class Preprocessing(object):
    
    def __init__(self):
        # load pre-processing matrices (for the uncertainty estimation)
        relpath = os.path.dirname(os.path.realpath(__file__))
        cc_preprocess_fname = os.path.join(relpath, PREPROCESS_FNAME)
        cc_preprocess = load(cc_preprocess_fname)
    
        self.v_r = matrix(cc_preprocess['v_r'])
        self.v_g = matrix(cc_preprocess['v_g'])
        self.C1 = matrix(cc_preprocess['C1'])
        self.C2 = matrix(cc_preprocess['C2'])
        self.C3 = matrix(cc_preprocess['C3'])
        self.G1 = matrix(cc_preprocess['G1'])
        self.G2 = matrix(cc_preprocess['G2'])
        self.G3 = matrix(cc_preprocess['G3'])
        self.S = matrix(cc_preprocess['S'])
        self.kegg_ids = cc_preprocess['cids']
        self.Nc = self.C1.shape[0]
        self.Ng = self.C3.shape[0]
        assert self.C1.shape[0] == self.C1.shape[1]
        assert self.C1.shape[1] == self.C2.shape[0]
        assert self.C2.shape[1] == self.C3.shape[0]
        assert self.C3.shape[0] == self.C3.shape[1]
        assert self.C3.shape[0] == self.C3.shape[1]
    
    def reaction_to_vectors(self, reaction):
        # x is the stoichiometric vector of the reaction, only for the
        # compounds that appeared in the original training set for CC
        x = matrix(zeros((self.Nc, 1)))

        # g is the group incidence vector of all the other compounds
        g = matrix(zeros((self.Ng, 1)))
        
        for kegg_id in reaction.kegg_ids():
            coeff = reaction.get_coeff(kegg_id)
            compound = reaction.get_compound(kegg_id)
            x += coeff * compound.get_stoich_vector(self.Nc)
            g += coeff * compound.get_group_incidence_vector(self.Ng)

        return x, g

    def reactions_to_matrices(self, reactions):
        """
            Arguments:
                reaction - a KeggReaction object
            
            Returns:
                X        - the stoichiometric matrix of the reactions (only
                           for compounds that appear in the original training
                           set of CC)
                G        - the group incidence matrix (of all other compounds)
        """
        X = matrix(zeros((self.Nc, len(reactions))))
        G = matrix(zeros((self.Ng, len(reactions))))
        
        for i, reaction in enumerate(reactions):
            x, g = self.reaction_to_vectors(reaction)
            X[:, i] = x
            G[:, i] = g

        return X, G

    def dG0_prime(self, reactions, pH=DEFAULT_PH, pMg=DEFAULT_PMG,
                      I=DEFAULT_IONIC_STRENGTH):
        
        if type(reactions) != types.ListType:
            reactions = [reactions]

        dG0_r_prime = matrix(map(lambda r: r.dG0_prime(pH, pMg, I), reactions)).T

        X, G = self.reactions_to_matrices(reactions)
        U = X.T * self.C1 * X + \
            X.T * self.C2 * G + \
            G.T * self.C2.T * X + \
            G.T * self.C3 * G
        
        return dG0_r_prime, U
            
    @staticmethod
    def WriteCompoundAndCoeff(kegg_id, coeff):
        if coeff == 1:
            return kegg_id
        else:
            return "%g %s" % (coeff, kegg_id)

    @staticmethod
    def DictToReactionString(d):
        """String representation."""
        left = []
        right = []
        for kegg_id, coeff in sorted(d.iteritems()):
            _s = Preprocessing.WriteCompoundAndCoeff(kegg_id, -coeff)
            if coeff < 0:
                left.append(_s)
            elif coeff > 0:
                right.append(_s)
        return "%s %s %s" % (' + '.join(left), '<=>', ' + '.join(right))

    @staticmethod
    def Analyze(self, x, g):
        weights_rc = x.T * self.G1
        weights_gc = x.T * self.G2 + g.T * self.G3
        weights = weights_rc + weights_gc

        res = []
        for j in xrange(self.S.shape[1]):
            d = {self.kegg_ids[i]: self.S[i, j]
                 for i in xrange(self.Nc)
                 if self.S[i, j] != 0}
            r_string = self.DictToReactionString(d)
            res.append({'w': weights[0, j],
                        'w_rc': weights_rc[0, j].round(4),
                        'w_gc': weights_gc[0, j].round(4),
                        'reaction_string': r_string})
        res.sort(key=lambda d: abs(d['w']), reverse=True)
        return res

    def IsUsingGroupContributions(self, x, g):
        weights_gc = x.T * self.G2 + g.T * self.G3
        sum_w_gc = sum(abs(weights_gc).flat)
        logging.info('sum(w_gc) = %.2g' % sum_w_gc)
        return sum_w_gc > 1e-5


if __name__ == '__main__':
    p = Preprocessing()
