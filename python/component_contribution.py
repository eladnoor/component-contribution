import sys
import numpy as np
from scipy.io import savemat
from training_data import TrainingData
from kegg_model import KeggModel
from compound_cacher import CompoundCacher
import inchi2gv

class ComponentContribution(object):
    
    def __init__(self, training_data):
        """
            Initialize G matrix, and then use the python script "inchi2gv.py" to decompose each of the 
            compounds that has an InChI and save the decomposition as a row in the G matrix.
        """
        self.train_cids = training_data.cids
        self.train_S = training_data.S
        self.train_b = np.matrix(training_data.dG0).T
        self.train_w = np.matrix(training_data.weight).T
        self.train_G = None
        self.train_S_joined = None
        self.model_S_joined = None
        self.params = None

    def savemat(self, fname):
        if self.params is None:
            raise Exception('One cannot call savemat() before calling train()')

        savemat(fname, self.params, oned_as='row')
    
    def estimate_kegg_model(self, model_S, model_cids):
        # standardize the CID list of the training data and the model
        # and create new (larger) matrices for each one
        cids_new = [cid for cid in model_cids if cid not in self.train_cids]
        self.cids_joined = self.train_cids + cids_new
        self.model_S_joined = ComponentContribution._zero_pad_S(
            model_S, model_cids, self.cids_joined)
        self.train_S_joined = ComponentContribution._zero_pad_S(
            self.train_S, self.train_cids, self.cids_joined)

        self.train_G = ComponentContribution.create_group_incidence_matrix(self.cids_joined)
        self.train()
        
        dG0_cc = self.params['dG0_cc']
        cov_dG0 = self.params['cov_dG0']
        MSE_kerG = self.params['MSE_kerG']
        
        model_dG0 = self.model_S_joined.T * dG0_cc
        model_cov_dG0 = self.model_S_joined.T * cov_dG0 * self.model_S_joined 

        return model_dG0, model_cov_dG0, MSE_kerG
    
    @staticmethod
    def _zero_pad_S(S, cids_orig, cids_joined):
        """
            takes a stoichiometric matrix with a given list of IDs 'cids' and adds
            0-rows so that the list of IDs will be 'cids_joined'
        """
        if not set(cids_orig).issubset(cids_joined):
            raise Exception('The full list is missing some IDs in "cids"')
        full_S = np.zeros((len(cids_joined), S.shape[1]))
        for i, cid in enumerate(cids_orig):
            full_S[cids_joined.index(cid), :] = S[i, :]
        
        return np.matrix(full_S)
        
    @staticmethod
    def create_group_incidence_matrix(cids):
        """
            create the group incidence matrix
        """
        ccache = CompoundCacher()
        groups_data = inchi2gv.init_groups_data()
        decomposer = inchi2gv.InChIDecomposer(groups_data)
        group_names = groups_data.GetGroupNames()
        
        Nc = len(cids)
        Ng = len(group_names)
        G = np.zeros((Nc, Ng))
        cpd_inds_without_gv = []
        
        # decompose the compounds in the training_data and add to G
        for i, cid in enumerate(cids):
            inchi_pH7 = ccache.get_kegg_compound(cid).inchi_pH7
            try:
                group_def = decomposer.inchi_to_groupvec(inchi_pH7)
                for j in xrange(len(group_names)):
                    G[i, j] = group_def[j]
            except inchi2gv.GroupDecompositionError:
                # for compounds that have no InChI or are not decomposable
                # add a unique 1 in a new column
                cpd_inds_without_gv.append(i)

        add_G = np.zeros((Nc, len(cpd_inds_without_gv)))
        for j, i in enumerate(cpd_inds_without_gv):
            add_G[i, j] = 1
        return np.matrix(np.hstack([G, add_G]))
    
    @staticmethod
    def _invert_project(A, eps=1e-10, method='octave'):
        n, m = A.shape
        if method == 'octave':
            from oct2py import Oct2Py
            oc = Oct2Py()
            U, S, V = oc.svd(A)
            s = np.diag(S)
            U = np.matrix(U)
            V = np.matrix(V)
            r = sum(abs(s) > eps)
            inv_S = np.matrix(np.diag([1.0/s[i] for i in xrange(r)]))
            inv_A = V[:, :r] * inv_S * U[:, :r].T
            P_R   = U[:, :r] * U[:, :r].T
            P_N   = U[:, r:] * U[:, r:].T
        elif method == 'numpy':
            # numpy.linalg.svd returns U, s, V_H such that
            # A = U * s * V_H
            # however, matlab and octave return U, S, V such that
            # V needs to be transposed when multiplied:
            # A = U * S * V.T
            U, s, V_H = np.linalg.svd(A, full_matrices=True)
            V = V_H.T
            r = sum(abs(s) > eps)
            inv_S = np.matrix(np.diag([1.0/s[i] for i in xrange(r)]))
            inv_A = V[:, :r] * inv_S * U[:, :r].T
            P_R   = U[:, :r] * U[:, :r].T
            P_N   = np.eye(n) - P_R
        elif method == 'nosvd':
            inv_A = A.T * np.linalg.inv(A * A.T + np.eye(n)*4e-6).T
            # then the solution for (A.T * x = b) will be given by (x = inv_A.T * b)
            P_R = A * inv_A
            P_N = np.eye(n) - P_R
            r = sum(np.abs(np.linalg.eig(P_R)[0]) > 0.5)
        else:
            raise ArgumentError('method argument must be "octave", "numpy" or "nosvd"')

        return inv_A, r, P_R, P_N
        
    def train(self):
        """
            Estimate standard Gibbs energies of formation
        """
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
        inv_S, r_rc, P_R_rc, P_N_rc = ComponentContribution._invert_project(S * W)

        # Linear regression for the group layer (aka GC)
        inv_GS, r_gc, P_R_gc, P_N_gc = ComponentContribution._invert_project(GS * W)

        # Calculate the contributions in the stoichiometric space
        dG0_rc = inv_S.T * W * b
        dG0_gc = G * inv_GS.T * W * b
        dG0_cc = P_R_rc * dG0_rc + P_N_rc * dG0_gc

        # Calculate the residual error (unweighted squared error divided by N - rank)
        e_rc = (S.T * dG0_rc - b)
        MSE_rc = float((e_rc.T * W * e_rc) / (n - r_rc))
        # MSE_rc = (e_rc.T * e_rc) / (n - r_rc)

        e_gc = (S.T * dG0_gc - b)
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
        inv_SWS, _, _, _ = ComponentContribution._invert_project(S * W * S.T)
        inv_GSWGS, _, _, _ = ComponentContribution._invert_project(GS * W * GS.T)


        #V_rc  = P_R_rc * (inv_S_orig.T * W * inv_S_orig) * P_R_rc
        #V_gc  = P_N_rc * G * (inv_GS_orig.T * W * inv_GS_orig) * G' * P_N_rc
        V_rc = P_R_rc * inv_SWS * P_R_rc
        V_gc  = P_N_rc * G * inv_GSWGS * G.T * P_N_rc
        # V_rc  = P_R_rc * (inv_S_orig.T * inv_S_orig) * P_R_rc
        # V_gc  = P_N_rc * G * (inv_GS_orig.T * inv_GS_orig) * G.T * P_N_rc
        V_inf = P_N_rc * G * P_N_gc * G.T * P_N_rc

        # Calculate the total of the contributions and covariances
        cov_dG0 = V_rc * MSE_rc + V_gc * MSE_gc + V_inf * MSE_inf

        # Put all the calculated data in 'params' for the sake of debugging
        self.params = {'b':           self.train_b,
                       'train_S':     self.train_S_joined,
                       'model_S':     self.model_S_joined,
                       'train_cids':  self.train_cids,
                       'cids':        self.cids_joined,
                       'w':           self.train_w,
                       'G':           self.train_G,
                       'dG0_rc':      dG0_rc,
                       'dG0_gc':      dG0_gc,
                       'dG0_cc':      dG0_cc,
                       'cov_dG0':     cov_dG0,
                       'V_rc':        V_rc,
                       'V_gc':        V_gc,
                       'V_inf':       V_inf,
                       'MSE_rc':      MSE_rc,
                       'MSE_gc':      MSE_gc,
                       'MSE_kerG':    MSE_kerG,
                       'MSE_inf':     MSE_inf,
                       'P_R_rc':      P_R_rc,
                       'P_R_gc':      P_R_gc,
                       'inv_S':       inv_S,
                       'inv_GS':      inv_GS,
                       'inv_SWS':     inv_SWS,
                       'inv_GSWGS':   inv_GSWGS}

if __name__ == '__main__':
    from kegg_model import KeggModel
    
    reaction_strings = sys.stdin.readlines()
    td = TrainingData()
    cc = ComponentContribution(td)
    model = KeggModel.from_formulas(reaction_strings)
    model.add_thermo(cc)
    
    dG0_prime, dG0_std = model.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
    
    sys.stdout.write('[' + 
                     ', '.join([str(x) for x in model.dG0.flat]) + '; ' + 
                     ', '.join([str(x) for x in dG0_prime.flat]) + 
                     ']')    
