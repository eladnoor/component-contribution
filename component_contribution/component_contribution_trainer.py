import os, logging
import numpy as np
from scipy.io import savemat, loadmat
from . import inchi2gv
from .training_data import TrainingData
from .kegg_reaction import KeggReaction
from .compound_cacher import CompoundCacher
from .thermodynamic_constants import default_T
from .molecule import Molecule, OpenBabelError
from .linalg import LINALG

base_path = os.path.split(os.path.realpath(__file__))[0]
CC_CACHE_FNAME = os.path.join(base_path, '../cache/component_contribution_python.mat')

class ComponentContribution(object):

    def __init__(self, training_data=None):
        if training_data is None:
            training_data = TrainingData()

        self.train_cids = list(training_data.cids)
        self.cids_joined = list(training_data.cids)

        self.train_S = training_data.S
        self.model_S_joined = np.matrix(self.train_S)
        self.train_S_joined = self.model_S_joined
        
        self.train_b = np.matrix(training_data.dG0).T
        self.train_w = np.matrix(training_data.weight).T
        self.train_G = None
        self.params = None

        self.ccache = CompoundCacher()
        self.groups_data = inchi2gv.init_groups_data()
        self.decomposer = inchi2gv.InChIDecomposer(self.groups_data)
        self.group_names = self.groups_data.GetGroupNames()
        
        self.Nc = len(self.cids_joined)
        self.Ng = len(self.group_names)

    @staticmethod
    def init():
        if os.path.exists(CC_CACHE_FNAME):
            logging.debug('Loading component-contributions from cache')
            return ComponentContribution.from_matfile(CC_CACHE_FNAME)
        else:
            logging.debug('Calculating the component-contributions from raw data')
            cc = ComponentContribution()
            cc.save_matfile(CC_CACHE_FNAME)
            return cc

    def save_matfile(self, file_name):
        if self.params is None:
            self.train()

        savemat(file_name, self.params, oned_as='row', do_compression=True)
    
    @staticmethod
    def from_matfile(file_name, training_data=None):
        cc = ComponentContribution(training_data=training_data)
        cc.params = loadmat(file_name)
        return cc
    
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
        
        if compound_id in self.cids_joined:
            i = self.cids_joined.index(compound_id)
            return self.params['dG0_cc'][i, 0]
        else:
            # Decompose the compound and calculate the 'formation energy'
            # using the group contributions.
            # Note that the length of the group contribution vector we get 
            # from CC is longer than the number of groups in "groups_data" 
            # since we artifically added fictive groups to represent all the 
            # non-decomposable compounds. Therefore, we truncate the 
            # dG0_gc vector since here we only use GC for compounds which
            # are not in cids_joined anyway.
            comp = self.ccache.get_compound(compound_id)
            try:
                group_vec = self.decomposer.smiles_to_groupvec(comp.smiles_pH7)
                g = np.matrix(group_vec.ToArray())
                dG0_gc = self.params['dG0_gc'][0:self.Ng, :]
                return float(np.dot(g, dG0_gc))
            except inchi2gv.GroupDecompositionError:
                return np.nan

    def _decompose_reaction(self, reaction):
        if self.params is None:
            self.train()
        
        cids = list(self.params['cids'])
        G = self.params['G']

        # calculate the reaction stoichiometric vector and the group incidence
        # vector (x and g)
        x = np.matrix(np.zeros((self.Nc, 1)))
        x_prime = []
        G_prime = []

        for compound_id, coeff in reaction.iteritems():
            if compound_id in self.cids_joined:
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
                comp = self.ccache.get_compound(compound_id)
                group_vec = self.decomposer.smiles_to_groupvec(comp.smiles_pH7)
                G_prime.append(group_vec.ToArray())

        if x_prime != []:
            g = np.matrix(x_prime) * np.vstack(G_prime)
        else:
            g = np.matrix(np.zeros((1, 1)))

        g.resize((G.shape[1], 1))

        return x, g

    def get_dG0_r(self, reaction, include_analysis=False):
        """
            Arguments:
                reaction - a KeggReaction object
            
            Returns:
                the CC estimation for this reaction's untransformed dG0 (i.e.
                using the major MS at pH 7 for each of the reactants)
        """
        try:
            x, g = self._decompose_reaction(reaction)
        except inchi2gv.GroupDecompositionError:
            if not include_analysis:
                return 0, 1e5
            else:
                return 0, 1e5, []

        v_r = np.matrix(self.params['preprocess_v_r'])
        v_g = np.matrix(self.params['preprocess_v_g'])
        C1  = np.matrix(self.params['preprocess_C1'])
        C2  = np.matrix(self.params['preprocess_C2'])
        C3  = np.matrix(self.params['preprocess_C3'])

        dG0_cc = float(x.T * v_r + g.T * v_g)
        s_cc_sqr = float(x.T * C1 * x + 2 * x.T * C2 * g + g.T * C3 * g)

        if not include_analysis:
            return dG0_cc, np.sqrt(s_cc_sqr)
        else:
            # Analyse the contribution of each training observation to this 
            # reaction's dG0 estimate.
            G1 = np.matrix(self.params['preprocess_G1'])
            G2 = np.matrix(self.params['preprocess_G2'])
            G3 = np.matrix(self.params['preprocess_G3'])
            S  = np.matrix(self.params['preprocess_S'])
            S_count = np.matrix(self.params['preprocess_S_count'])
            cids = self.params['cids']
            
            # dG0_cc = (x*G1 + x*G2 + g*G3)*b
            weights_rc = (x.T * G1).round(5)
            weights_gc = (x.T * G2 + g.T * G3).round(5)
            weights = weights_rc + weights_gc
    
            orders = sorted(range(weights.shape[1]),
                            key=lambda j:abs(weights[0, j]), reverse=True)
    
            analysis = []        
            for j in orders:
                if abs(weights[0, j]) < 1e-5:
                    continue
                r = KeggReaction({cids[i]:S[i,j] for i in xrange(S.shape[0])
                                  if S[i,j] != 0})
                analysis.append({'index': j,
                                 'w_rc': weights_rc[0, j],
                                 'w_gc': weights_gc[0, j],
                                 'reaction': r,
                                 'count': int(S_count[0, j])})

            return dG0_cc, np.sqrt(s_cc_sqr), analysis

    def get_dG0_r_multi(self, reactions):
        """
            Arguments:
                reaction - a KeggReaction object
            
            Returns:
                the CC estimation for this reaction's untransformed dG0 (i.e.
                using the major MS at pH 7 for each of the reactants)
        """
        X = []
        G = []
        for reaction in reactions:
            try:
                x, g = self._decompose_reaction(reaction)
            except inchi2gv.GroupDecompositionError:
                x = np.zeros((self.Nc, 1))
                g = np.zeros((self.params['G'].shape[1], 1))
            X.append(list(x.flat))
            G.append(list(g.flat))
        X = np.matrix(X).T
        G = np.matrix(G).T
        
        v_r = np.matrix(self.params['preprocess_v_r'])
        v_g = np.matrix(self.params['preprocess_v_g'])
        C1  = np.matrix(self.params['preprocess_C1'])
        C2  = np.matrix(self.params['preprocess_C2'])
        C3  = np.matrix(self.params['preprocess_C3'])

        dG0_cc = X.T * v_r + G.T * v_g
        U = X.T * C1 * X + X.T * C2 * G + G.T * C2.T * X + G.T * C3 * G
        return dG0_cc, U
        
    def get_compound_json(self, compound_id):
        """
            adds the component-contribution estimation to the JSON
        """
        if compound_id is None:
            raise ValueError('given compound ID is None')
        if self.params is None:
            self.train()

        d = {'CID': compound_id}
        comp = self.ccache.get_compound(compound_id)
        gv = None
        
        if compound_id in self.cids_joined:
            i = self.cids_joined.index(compound_id)
            gv = self.params['G'][i, :]
            major_ms_dG0_f = self.params['dG0_cc'][i, 0]
            d['compound_index'] = i
        elif comp.smiles_pH7 is not None:
            # decompose the compounds in the training_data and add to G
            try:
                group_def = self.decomposer.smiles_to_groupvec(comp.smiles_pH7)
                gv = np.matrix(group_def.ToArray())
                # we need to truncate the dG0_gc matrix from all the group
                # dimensions that correspond to non-decomposable compounds
                # from the training set
                dG0_gc = self.params['dG0_gc'][0:self.Ng, :]
                major_ms_dG0_f = float(np.dot(gv, dG0_gc))
            except inchi2gv.GroupDecompositionError:
                d['error'] = 'We cannot estimate the formation energy of this compound ' +\
                             'because its structure is too small or too complex to ' +\
                             'decompose to groups'
                major_ms_dG0_f = np.nan
        else:
            d['error'] = 'We cannot estimate the formation energy of this compound ' +\
                         'because it has no defined structure'
            major_ms_dG0_f = np.nan

        if gv is not None:
            sparse_gv = filter(lambda x: x[1] != 0, enumerate(gv.flat))
            d['group_vector'] = sparse_gv

        if not np.isnan(major_ms_dG0_f):
            d['pmap'] = {'source': 'Component Contribution (2013)',
                         'species': list(comp.get_species(major_ms_dG0_f, default_T))}

        d['num_electrons'] = comp.atom_bag.get('e-', 0)

        if comp.inchi is not None:
            d['InChI'] = comp.inchi
            try:
                mol = Molecule.FromInChI(str(comp.inchi))
                d['mass'] = mol.GetExactMass()
                d['formula'] = mol.GetFormula()
            except OpenBabelError:
                if compound_id == 'C00282': # an exception for hydrogen
                    d['mass'] = 2.0157
                    d['formula'] = 'H2'
                else:
                    d['mass'] = 0
                    d['formula'] = ''
            
        return d
    
    def estimate_kegg_model(self, model_S, model_cids):
    
        # standardize the CID list of the training data and the model
        # and create new (larger) matrices for each one
        cids_new = [cid for cid in model_cids if cid not in self.train_cids]

        self.cids_joined += cids_new
        self.Nc = len(self.cids_joined)
                
        self.model_S_joined = LINALG._zero_pad_S(
            model_S, model_cids, self.cids_joined)

        self.train_S_joined = LINALG._zero_pad_S(
            self.train_S, self.train_cids, self.cids_joined)

        self.train()
        
        dG0_cc = self.params['dG0_cc']
        cov_dG0 = self.params['cov_dG0']
        MSE_kerG = self.params['MSE_kerG']
        
        model_dG0 = self.model_S_joined.T * dG0_cc
        model_cov_dG0 = self.model_S_joined.T * cov_dG0 * self.model_S_joined 

        return model_dG0, model_cov_dG0, MSE_kerG
    
    def create_group_incidence_matrix(self):
        """
            Initialize G matrix, and then use the python script "inchi2gv.py" to
            decompose each of the compounds that has an InChI and save the
            decomposition as a row in the G matrix.
        """

        G = np.zeros((self.Nc, self.Ng))
        cpd_inds_without_gv = []
        
        # decompose the compounds in the training_data and add to G
        for i, compound_id in enumerate(self.cids_joined):
            smiles_pH7 = self.ccache.get_compound(compound_id).smiles_pH7
            try:
                group_def = self.decomposer.smiles_to_groupvec(smiles_pH7)
                for j in xrange(len(self.group_names)):
                    G[i, j] = group_def[j]
            except inchi2gv.GroupDecompositionError:
                # for compounds that have no InChI or are not decomposable
                # add a unique 1 in a new column
                cpd_inds_without_gv.append(i)

        N_non_decomposable = len(cpd_inds_without_gv)
        add_G = np.zeros((self.Nc, N_non_decomposable))
        for j, i in enumerate(cpd_inds_without_gv):
            add_G[i, j] = 1
        return np.matrix(np.hstack([G, add_G]))
    
    def train(self):
        """
            Estimate standard Gibbs energies of formation
        """
        self.train_G = self.create_group_incidence_matrix()

        S = self.train_S_joined
        G = self.train_G
        b = self.train_b
        w = self.train_w
        
        m, n = S.shape
        assert G.shape[0] == m
        assert b.shape == (n, 1)
        assert w.shape == (n, 1)

        # Apply weighing
        W = np.diag(w.flat)
        GS = G.T * S

        # Linear regression for the reactant layer (aka RC)
        inv_S, r_rc, P_R_rc, P_N_rc = LINALG._invert_project(S * W)

        # Linear regression for the group layer (aka GC)
        inv_GS, r_gc, P_R_gc, P_N_gc = LINALG._invert_project(GS * W)

        # calculate the group contributions
        dG0_gc = inv_GS.T * W * b

        # Calculate the contributions in the stoichiometric space
        dG0_rc = inv_S.T * W * b
        dG0_cc = P_R_rc * dG0_rc + P_N_rc * G * dG0_gc

        # Calculate the residual error (unweighted squared error divided by N - rank)
        e_rc = (S.T * dG0_rc - b)
        MSE_rc = float((e_rc.T * W * e_rc) / (n - r_rc))
        # MSE_rc = (e_rc.T * e_rc) / (n - r_rc)

        e_gc = (GS.T * dG0_gc - b)
        MSE_gc = float((e_gc.T * W * e_gc) / (n - r_gc))
        # MSE_gc = (e_gc.T * e_gc) / (n - r_gc)

        # Calculate the MSE of GC residuals for all reactions in ker(G).
        # This will help later to give an estimate of the uncertainty for such
        # reactions, which otherwise would have a 0 uncertainty in the GC method.
        kerG_inds = list(np.where(np.all(GS == 0, 0))[1].flat)
        
        e_kerG = e_gc[kerG_inds]
        MSE_kerG = float((e_kerG.T * e_kerG) / len(kerG_inds))

        MSE_inf = 1e10

        # Calculate the uncertainty covariance matrices
        # [inv_S_orig, ~, ~, ~] = invertProjection(S);
        # [inv_GS_orig, ~, ~, ~] = invertProjection(GS);
        inv_SWS, _, _, _ = LINALG._invert_project(S * W * S.T)
        inv_GSWGS, _, _, _ = LINALG._invert_project(GS * W * GS.T)


        #V_rc  = P_R_rc * (inv_S_orig.T * W * inv_S_orig) * P_R_rc
        #V_gc  = P_N_rc * G * (inv_GS_orig.T * W * inv_GS_orig) * G' * P_N_rc
        V_rc = P_R_rc * inv_SWS * P_R_rc
        V_gc  = P_N_rc * G * inv_GSWGS * G.T * P_N_rc
        # V_rc  = P_R_rc * (inv_S_orig.T * inv_S_orig) * P_R_rc
        # V_gc  = P_N_rc * G * (inv_GS_orig.T * inv_GS_orig) * G.T * P_N_rc
        V_inf = P_N_rc * G * P_N_gc * G.T * P_N_rc

        # Calculate the total of the contributions and covariances
        cov_dG0 = V_rc * MSE_rc + V_gc * MSE_gc + V_inf * MSE_inf

        # preprocessing matrices (for calculating the contribution of each 
        # observation)
        G1 = P_R_rc * inv_S.T * W
        G2 = P_N_rc * G * inv_GS.T * W
        G3 = inv_GS.T * W
        
        S_uniq, P_col = LINALG._col_uniq(S)
        S_counter = np.sum(P_col, 0)
        preprocess_G1 = G1 * P_col
        preprocess_G2 = G2 * P_col
        preprocess_G3 = G3 * P_col

        # preprocessing matrices (for quick calculation of uncertainty)
        preprocess_C1 = cov_dG0
        preprocess_C2 = MSE_gc * P_N_rc * G * inv_GSWGS + MSE_inf * G * P_N_gc
        preprocess_C3 = MSE_gc * inv_GSWGS + MSE_inf * P_N_gc

        # Put all the calculated data in 'params' for the sake of debugging
        self.params = {'b':              self.train_b,
                       'train_S':        self.train_S_joined,
                       'model_S':        self.model_S_joined,
                       'train_cids':     self.train_cids,
                       'cids':           self.cids_joined,
                       'w':              self.train_w,
                       'G':              self.train_G,
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

