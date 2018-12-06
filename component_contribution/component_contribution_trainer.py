import gzip
import json
import logging
import os

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from . import Reaction, ccache, default_T, inchi2gv
from .linalg import LINALG
from .molecule import Molecule, OpenBabelError
from .training_data import FullTrainingData


logger = logging.getLogger(__name__)


class ComponentContribution(object):

    MSE_inf = 1e10
    CACHE_PARAM_FNAME = \
        resource_filename('component_contribution',
                          'cache/component_contribution.npz')

    def __init__(self, cache_file_name=CACHE_PARAM_FNAME, training_data=None):
        self.groups_data = inchi2gv.init_groups_data()
        self.decomposer = inchi2gv.InChIDecomposer(self.groups_data)
        self.group_names = self.groups_data.GetGroupNames()

        if cache_file_name and os.path.exists(cache_file_name):
            if training_data is not None:
                logger.warning('You provided both an existing cache file, and '
                               'a training dataset for ComponentContribution. '
                               'The training data will be ignored.')

            logger.info('Loading component-contributions params from: %s'
                        % cache_file_name)
            self.params = self.load_params(cache_file_name)
            self.cids = self.params['cids'].tolist()
            self.train_S = self.params['S']
            self.train_b = self.params['b']
            self.train_w = self.params['w']
            self.train_G = self.params['G']
            self.Nc = len(self.cids)
            self.Ng = len(self.group_names)
        else:
            logger.info('Calculating the component-contributions from raw '
                        'data')
            if training_data is None:
                training_data = FullTrainingData()

            self.params = None
            self.cids = training_data.cids
            self.train_S = training_data.stoichiometric_matrix
            self.train_b = training_data.dG0
            self.train_w = training_data.weight
            self.train_G = self.create_group_incidence_matrix()
            self.Nc = len(self.cids)
            self.Ng = len(self.group_names)
            self.params = self.train(self.train_S, self.train_G,
                                     self.train_b, self.train_w)
            self.params['cids'] = self.cids

            if cache_file_name:
                logger.info('Saving component-contributions params to: %s'
                            % cache_file_name)
                self.save_params(self.params, cache_file_name)

    @staticmethod
    def save_params(params, file_name):
        np.savez_compressed(file_name, **params)

    @staticmethod
    def load_params(file_name):
        params = np.load(file_name)
        return dict(params)

    def save_equilibrator_bin_data(self, npz_file_name):
        """
            write an NPZ file (Numpy binary files) containing only the
            preprocessed data needed for running Component Contribution
            estimations
        """
        preprocess_dict = {'cids': self.params['cids']}
        for k, v in self.params.items():
            if k.find('preprocess_') != -1:
                preprocess_dict[k.replace('preprocess_', '')] = v
        np.savez_compressed(npz_file_name, **preprocess_dict)

    def save_equilibrator_metadata(self, json_file_name):
        # write the JSON file containing the 'additiona' data on all the
        # compounds in eQuilibrator (i.e. formula, mass, pKa values, etc.)
        compound_json = []
        for i, compound_id in enumerate(ccache.all_compound_ids()):
            logger.debug("exporting " + compound_id)

            # skip compounds that cause a segmentation fault in openbabel
            if compound_id in ['KEGG:C09078', 'KEGG:C09093', 'KEGG:C09145',
                               'KEGG:C09246', 'KEGG:C10282', 'KEGG:C10286',
                               'KEGG:C10356', 'KEGG:C10359', 'KEGG:C10396',
                               'KEGG:C16818', 'KEGG:C16839', 'KEGG:C16857']:
                continue
            d = self.to_dict(compound_id)
            # override H2O as liquid phase only
            if compound_id in ['KEGG:C00001']:
                d['pmap']['species'][0]['phase'] = 'liquid'
            # override Sulfur as solid phase only
            if compound_id in ['KEGG:C00087']:
                d['pmap']['species'][0]['phase'] = 'solid'
            # add gas phase for O2 and N2
            if compound_id in ['KEGG:C00007', 'KEGG:C00697']:
                d['pmap']['species'].append({'phase': 'gas', 'dG0_f': 0})
            # add gas phase for H2
            if compound_id in ['KEGG:C00282']:
                d['pmap']['species'].append({'phase': 'gas',
                                             'dG0_f': 0,
                                             'nH': 2})
            # add gas phase for CO2
            if compound_id in ['KEGG:C00011']:
                d['pmap']['species'].append({'phase': 'gas',
                                             'dG0_f': -394.36})
            # add gas phase for CO
            if compound_id in ['KEGG:C00237']:
                d['pmap']['species'].append({'phase': 'gas',
                                             'dG0_f': -137.17})

            compound_json.append(d)

        json_bytes = json.dumps(compound_json, sort_keys=True, indent=4)
        with gzip.open(json_file_name, 'w') as new_json:
            new_json.write(json_bytes.encode('utf-8'))

    def get_major_ms_dG0_f(self, compound_id):
        """
            Returns the chemical formation energy of the major MS at pH 7.
            If the compound is part of the training set, returns the value
            that was calculated during training. Otherwise, we use pure
            group contribution (if possible) on the groups of the major MS.
        """
        if compound_id is None:
            raise ValueError('given compound ID is None')
        if self.params is None:
            self.train()

        if compound_id in self.cids:
            i = self.cids.index(compound_id)
            return self.params['dG0_cc'][i, 0]
        else:
            # Decompose the compound and calculate the 'formation energy'
            # using the group contributions.
            # Note that the length of the group contribution vector we get
            # from CC is longer than the number of groups in "groups_data"
            # since we artifically added fictive groups to represent all the
            # non-decomposable compounds. Therefore, we truncate the
            # dG0_gc vector since here we only use GC for compounds which
            # are not in cids anyway.
            comp = ccache.get_compound(compound_id)
            try:
                group_vec = self.decomposer.smiles_to_groupvec(comp.smiles)
                g = group_vec.as_array()
                dG0_gc = self.params['dG0_gc'][0:self.Ng]
                return g @ dG0_gc
            except inchi2gv.GroupDecompositionError:
                return np.nan

    def _decompose_reaction(self, reaction, raise_exception=True):
        # calculate the reaction stoichiometric vector and the group incidence
        # vector (x and g)
        x = np.zeros(self.Nc)
        total_gv = np.zeros(self.Ng)

        for compound_id, coeff in reaction.items():
            if compound_id in ['C00080', 'KEGG:C00080']:
                continue
            if compound_id in self.cids:
                i = self.cids.index(compound_id)
                x[i] = coeff
            else:
                # Decompose the compound and calculate the 'formation energy'
                # using the group contributions.
                # Note that the length of the group contribution vector we get
                # from CC is longer than the number of groups in "groups_data"
                # since we artifically added fictive groups to represent all
                # the non-decomposable compounds. Therefore, we truncate the
                # dG0_gc vector since here we only use GC for compounds which
                # are not in cids anyway.
                comp = ccache.get_compound(compound_id)
                try:
                    _gv = self.decomposer.smiles_to_groupvec(comp.smiles)
                    total_gv += _gv.as_array()
                except inchi2gv.GroupDecompositionError as exception:
                    if raise_exception:
                        logger.warning('Compound %s cannot be decomposed and '
                                       'is also not in the training set'
                                       % compound_id)
                        raise exception
                    else:
                        return np.zeros(self.Nc), np.zeros(self.Ng)

        return x, total_gv

    def get_dG0_r(self, reaction, raise_exception=False,
                  include_analysis=False):
        """
            Arguments:
            - reaction (a Reaction object)

            Returns:
            - dG0_cc
                        the CC estimation for this reaction's untransformed
                        dG0 (i.e. using the major MS at pH 7 for each of the
                        reactants)
            - sigma_cc
                        the standard error of the estimate. multiply by
                        1.96 to get the 95% confidence interval.
        """
        try:
            dG0_cc, U = self.get_dG0_r_multi([reaction], raise_exception=True)
            dG0_cc = dG0_cc[0]
            sigma_cc = np.sqrt(U[0, 0])
        except inchi2gv.GroupDecompositionError as exception:
            if raise_exception:
                raise exception
            if not include_analysis:
                return 0, np.sqrt(self.params['MSE_inf'])
            else:
                return 0, np.sqrt(self.params['MSE_inf']), []

        if not include_analysis:
            return dG0_cc, sigma_cc
        else:
            # Analyse the contribution of each training observation to this
            # reaction's dG0 estimate.
            G1 = self.params['preprocess_G1']
            G2 = self.params['preprocess_G2']
            G3 = self.params['preprocess_G3']
            S = self.params['preprocess_S']
            S_count = self.params['preprocess_S_count']
            cids = self.params['cids']

            x, g = self._decompose_reaction(reaction)

            # dG0_cc = (x*G1 + x*G2 + g*G3)*b
            weights_rc = (x @ G1).round(5)
            weights_gc = (x @ G2 + g @ G3[0:self.Ng, :]).round(5)
            weights = weights_rc + weights_gc

            orders = sorted(range(weights.shape[1]),
                            key=lambda j: abs(weights[0, j]), reverse=True)

            analysis = []
            for j in orders:
                if abs(weights[0, j]) < 1e-5:
                    continue
                r = Reaction({cids[i]: S[i, j] for i in range(S.shape[0])
                              if S[i, j] != 0})
                analysis.append({'index': j,
                                 'w_rc': weights_rc[0, j],
                                 'w_gc': weights_gc[0, j],
                                 'reaction': r,
                                 'count': int(S_count[0, j])})

            return dG0_cc, sigma_cc, analysis

    def get_dG0_r_multi(self, reactions, raise_exception=False):
        """
            Arguments:
            - reactions (a list of Reaction objects)
            - raise_exception (flag for non-decomposable compounds)

            Returns:
            - dG0
                      a 1D NumPy array containing the CC estimates for
                      the reactions' untransformed dG0
                      (i.e. using the major MS at pH 7 for each of the
                      reactants)
            - U
                      a 2D numpy array containing the covariance matrix
                      of the standard errors of the
                      estimates. one can use the eigenvectors of the matrix
                      to define a confidence high-dimensional space, or use
                      U as the covariance of a Gaussian used for sampling
                      (where dG0_cc is the mean of that Gaussian).
        """
        Nr = len(reactions)
        X = np.zeros((self.Nc, Nr))
        G = np.zeros((self.train_G.shape[1], Nr))
        for i, reaction in enumerate(reactions):
            x, g = self._decompose_reaction(reaction,
                                            raise_exception=raise_exception)
            X[:, i] = x
            G[:self.Ng, i] = g

        v_r = self.params['preprocess_v_r']
        v_g = self.params['preprocess_v_g']
        C1 = self.params['preprocess_C1']
        C2 = self.params['preprocess_C2']
        C3 = self.params['preprocess_C3']

        dG0 = X.T @ v_r + G.T @ v_g
        U = X.T @ C1 @ X + X.T @ C2 @ G + G.T @ C2.T @ X + G.T @ C3 @ G
        return dG0, U

    def get_dG0_r_prime(self, reaction, pH, ionic_strength, T,
                        raise_exception=False):
        """
            Arguments:
            - reaction (a Reaction objects)
            - pH
            - ionic_strength (in M)
            - T (temperature in Kalvin)
            - raise_exception (flag for non-decomposable compounds)

            Returns:
            - dG0_prime (CC estimates for the reaction's transformed dG0 in
                         kJ/mol)
            - sigma_cc (the standard error of the estimate. multiply by
                        1.96 to get the 95% confidence interval)
        """
        dG0_cc, sigma_cc = self.get_dG0_r(reaction, raise_exception)
        dG0_prime = dG0_cc + reaction.get_transform_ddG0(pH, ionic_strength, T)
        return dG0_prime, sigma_cc

    def get_dG0_r_prime_multi(self, reactions, pH, ionic_strength, T,
                              raise_exception=False):
        """
            Arguments:
            - reactions (a list of Reaction objects)
            - pH
            - ionic_strength (in M)
            - T (temperature in Kalvin)
            - raise_exception (flag for non-decomposable compounds)

            Returns:
            - dG0_prime
                      a 1D NumPy array containing the CC estimates for
                      the reactions' transformed dG0
            - U
                      a 2D numpy array containing the covariance matrix
                      of the standard errors of the
                      estimates. one can use the eigenvectors of the matrix
                      to define a confidence high-dimensional space, or use
                      U as the covariance of a Gaussian used for sampling
                      (where dG0_cc is the mean of that Gaussian).
        """
        dG0, U = self.get_dG0_r_multi(reactions, raise_exception)
        ddG0 = np.array([r.get_transform_ddG0(pH, ionic_strength, T)
                         for r in reactions])
        dG0_prime = dG0 + ddG0
        return dG0_prime, U

    def to_dict(self, compound_id):
        """
            adds the component-contribution estimation to the JSON
        """
        if compound_id is None:
            raise ValueError('given compound ID is None')
        if self.params is None:
            self.train()

        comp = ccache.get_compound(compound_id)
        d = {'compound_id': compound_id, 'inchi_key': comp.inchi_key}
        gv = None

        if compound_id in self.cids:
            i = self.cids.index(compound_id)
            gv = self.params['G'][i, :]
            major_ms_dG0_f = self.params['dG0_cc'][i]
            d['compound_index'] = i
        elif comp.smiles is not None:
            # decompose the compounds in the training_data and add to G
            try:
                group_def = self.decomposer.smiles_to_groupvec(comp.smiles)
                gv = group_def.as_array()
                # we need to truncate the dG0_gc matrix from all the group
                # dimensions that correspond to non-decomposable compounds
                # from the training set
                dG0_gc = self.params['dG0_gc'][0:self.Ng]
                major_ms_dG0_f = float(gv @ dG0_gc)
            except inchi2gv.GroupDecompositionError:
                d['error'] = ('We cannot estimate the formation energy of '
                              'this compound because its structure is too '
                              'small or too complex to decompose to groups')
                major_ms_dG0_f = np.nan
        else:
            d['error'] = ('We cannot estimate the formation energy of this '
                          'compound because it has no defined structure')
            major_ms_dG0_f = np.nan

        if gv is not None:
            sparse_gv = [(i, int(g)) for (i, g) in enumerate(gv.flat)
                         if g != 0]
            d['group_vector'] = sparse_gv

        if not np.isnan(major_ms_dG0_f):
            _species = list(comp.get_species(major_ms_dG0_f, default_T))
            d['pmap'] = {'source': 'Component Contribution (2013)',
                         'species': _species}

        d['num_electrons'] = comp.atom_bag.get('e-', 0)

        if comp.inchi is not None:
            d['InChI'] = comp.inchi
            try:
                mol = Molecule.FromInChI(str(comp.inchi))
                d['mass'] = mol.GetExactMass()
                d['formula'] = mol.GetFormula()
            except OpenBabelError:
                if compound_id == 'C00282':  # an exception for hydrogen
                    d['mass'] = 2.0157
                    d['formula'] = 'H2'
                else:
                    d['mass'] = 0
                    d['formula'] = ''

        return d

    def create_group_incidence_matrix(self):
        """
            Initialize G matrix, and then use the python script "inchi2gv.py"
            to decompose each of the compounds that has an InChI and save the
            decomposition as a row in the G matrix.
        """

        gv_data = []
        # decompose the compounds in the training_data and add to G
        for compound_id in self.cids:
            smiles = ccache.get_compound(compound_id).smiles
            try:
                gv_data.append(
                        list(self.decomposer.smiles_to_groupvec(smiles).flat))
            except inchi2gv.GroupDecompositionError:
                gv_data.append([0] * len(self.group_names))

        G = pd.DataFrame(index=self.cids,
                         columns=self.group_names,
                         dtype=float,
                         data=gv_data)

        for compound_id in G.index[(G == 0).all(1)]:
            # add a column for this compound, representing itself
            # as a new group
            G[compound_id] = 0.0

            # place a single '1' for this compound group decomposition
            G.at[compound_id, compound_id] = 1.0

        return G.values

    @staticmethod
    def train(S, G, b, w):
        """
            Estimate standard Gibbs energies of formation
        """
        assert type(S) == np.ndarray
        assert type(G) == np.ndarray
        assert type(b) == np.ndarray
        assert type(w) == np.ndarray

        m, n = S.shape
        assert G.shape[0] == m
        assert b.shape == (n,)
        assert w.shape == (n,)

        # Apply weighing
        W = np.diag(w.flat)
        GS = G.T @ S

        # Linear regression for the reactant layer (aka RC)
        inv_S, r_rc, P_R_rc, P_N_rc = LINALG._invert_project(S @ W)

        # Linear regression for the group layer (aka GC)
        inv_GS, r_gc, P_R_gc, P_N_gc = LINALG._invert_project(GS @ W)

        # calculate the group contributions
        dG0_gc = np.squeeze(inv_GS.T @ W @ b)

        # Calculate the contributions in the stoichiometric space
        dG0_rc = np.squeeze(inv_S.T @ W @ b)
        dG0_cc = np.squeeze(P_R_rc @ dG0_rc + P_N_rc @ G @ dG0_gc)

        # Calculate the residual error (unweighted squared error divided
        # by N - rank)
        e_rc = (S.T @ dG0_rc - b)
        MSE_rc = float((e_rc.T @ W @ e_rc) / (n - r_rc))
        # MSE_rc = (e_rc.T @ e_rc) / (n - r_rc)

        e_gc = (GS.T @ dG0_gc - b)
        MSE_gc = float((e_gc.T @ W @ e_gc) / (n - r_gc))
        # MSE_gc = (e_gc.T @ e_gc) / (n - r_gc)

        # Calculate the MSE of GC residuals for all reactions in ker(G).
        # This will help later to give an estimate of the uncertainty for such
        # reactions, which otherwise would have a 0 uncertainty in the GC
        # method.
        kerG_inds = list(np.where(np.all(GS == 0, 0))[0].flat)

        e_kerG = e_gc[kerG_inds]
        MSE_kerG = float((e_kerG.T @ e_kerG) / len(kerG_inds))

        MSE_inf = ComponentContribution.MSE_inf

        # Calculate the uncertainty covariance matrices
        inv_SWS, _, _, _ = LINALG._invert_project(S @ W @ S.T)
        inv_GSWGS, _, _, _ = LINALG._invert_project(GS @ W @ GS.T)

        V_rc = P_R_rc @ inv_SWS @ P_R_rc
        V_gc = P_N_rc @ G @ inv_GSWGS @ G.T @ P_N_rc
        V_inf = P_N_rc @ G @ P_N_gc @ G.T @ P_N_rc

        # Calculate the total of the contributions and covariances
        cov_dG0 = MSE_rc * V_rc + MSE_gc * V_gc + MSE_inf * V_inf

        # preprocessing matrices (for calculating the contribution of each
        # observation)
        G1 = P_R_rc @ inv_S.T @ W
        G2 = P_N_rc @ G @ inv_GS.T @ W
        G3 = inv_GS.T @ W

        S_uniq, P_col = LINALG._col_uniq(S)
        S_counter = np.sum(P_col, 0)
        preprocess_G1 = G1 @ P_col
        preprocess_G2 = G2 @ P_col
        preprocess_G3 = G3 @ P_col

        # preprocessing matrices (for quick calculation of uncertainty)
        preprocess_C1 = cov_dG0
        preprocess_C2 = MSE_gc * P_N_rc @ G @ inv_GSWGS + MSE_inf * G @ P_N_gc
        preprocess_C3 = MSE_gc * inv_GSWGS + MSE_inf * P_N_gc

        # Put all the calculated data in 'params' for the sake of debugging
        params  = {'b':              b,
                   'S':              S,
                   'w':              w,
                   'G':              G,
                   'dG0_rc':         dG0_rc,
                   'dG0_gc':         dG0_gc,
                   'dG0_cc':         dG0_cc,
                   'cov_dG0':        cov_dG0,
                   'V_rc':           V_rc,
                   'V_gc':           V_gc,
                   'V_inf':          V_inf,
                   'MSE_rc':         MSE_rc,
                   'MSE_gc':         MSE_gc,
                   'MSE_kerG':       MSE_kerG,
                   'MSE_inf':        MSE_inf,
                   'P_R_rc':         P_R_rc,
                   'P_R_gc':         P_R_gc,
                   'P_N_rc':         P_N_rc,
                   'P_N_gc':         P_N_gc,
                   'inv_S':          inv_S,
                   'inv_GS':         inv_GS,
                   'inv_SWS':        inv_SWS,
                   'inv_GSWGS':      inv_GSWGS,
                   'preprocess_v_r': dG0_cc,
                   'preprocess_v_g': dG0_gc,
                   'G1':             G1,
                   'G2':             G2,
                   'G3':             G3,
                   'preprocess_G1':  preprocess_G1,
                   'preprocess_G2':  preprocess_G2,
                   'preprocess_G3':  preprocess_G3,
                   'preprocess_S':   S_uniq,
                   'preprocess_S_count': S_counter,
                   'preprocess_C1':  preprocess_C1,
                   'preprocess_C2':  preprocess_C2,
                   'preprocess_C3':  preprocess_C3}
        return params


if __name__ == '__main__':
    logger.setLevel(logging.INFO)

    if os.path.exists(ComponentContribution.CACHE_PARAM_FNAME):
        logger.info('Found existing CC cache file. Deleting it and retraining '
                    'Compounent Contribution.')
        os.remove(ComponentContribution.CACHE_PARAM_FNAME)
    cc = ComponentContribution()
