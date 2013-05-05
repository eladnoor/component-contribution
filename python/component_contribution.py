import sys
import numpy as np
from inchi2gv import GroupsData, InChI2GroupVector, GROUP_CSV, GroupDecompositionError
from training_data import TrainingData
from kegg_model import KeggModel
from compound_cacher import CompoundCacher

class ComponentContribution(object):
    
    def __init__(self, training_data):
        """
            Initialize G matrix, and then use the python script "inchi2gv.py" to decompose each of the 
            compounds that has an InChI and save the decomposition as a row in the G matrix.
        """
        self.ccache = CompoundCacher.getInstance()
        
        self.groups_data = GroupsData.FromGroupsFile(GROUP_CSV, transformed=False)
        self.inchi2gv = InChI2GroupVector(self.groups_data)
        self.group_names = self.groups_data.GetGroupNames()
        
        self.training_data = training_data

    
    def estimate_kegg_model(self, kegg_model):
        # standardize the CID list of the training data and the model
        # and create new (larger) matrices for each one
        cids = sorted(set(kegg_model.cids + self.training_data.cids))
        model_S = ComponentContribution._zero_pad_S(
            kegg_model.S, kegg_model.cids, cids)
        train_S = ComponentContribution._zero_pad_S(
            self.training_data.S, self.training_data.cids, cids)

        G = self.create_group_incidence_matrix(cids)
        
        dG0_f, cov_dG0, params = self.train(train_S, G)
        
        model_dG0 = model_S.T * dG0_f
        model_cov_dG0 = model_S.T * cov_dG0 * model_S 

        return model_dG0, model_cov_dG0
    
    @staticmethod
    def _zero_pad_S(S, cids, all_cids):
        """
            takes a stoichiometric matrix with a given list of IDs 'cids' and adds
            0-rows so that the list of IDs will be 'all_cids'
        """
        if not set(cids).issubset(all_cids):
            raise Exception('The full list is missing some IDs in "cids"')
        full_S = np.zeros((len(all_cids), S.shape[1]))
        for i, cid in enumerate(cids):
            full_S[all_cids.index(cid), :] = S[i, :]
        
        return np.matrix(full_S)
        
    def create_group_incidence_matrix(self, cids):
        """
            create the group incidence matrix
        """
        Nc = len(cids)
        Ng = len(self.group_names)
        G = np.zeros((Nc, Ng))
        cpd_inds_without_gv = []
        
        # decompose the compounds in the training_data and add to G
        for i, cid in enumerate(cids):
            inchi = self.ccache.get_kegg_compound(cid).inchi
            try:
                group_def = self.inchi2gv.EstimateInChI(inchi)
                for j in xrange(len(self.group_names)):
                    G[i, j] = group_def[0, j]
            except GroupDecompositionError as e:
                # for compounds that have no InChI or are not decomposable
                # add a unique 1 in a new column
                cpd_inds_without_gv.append(i)

        add_G = np.zeros((Nc, len(cpd_inds_without_gv)))
        for j, i in enumerate(cpd_inds_without_gv):
            add_G[i, j] = 1
        return np.matrix(np.hstack([G, add_G]))
    
    @staticmethod
    def _invert_project(A, eps=1e-10):
        U, s, V = np.linalg.svd(A, full_matrices=True)
        r = (s > eps).nonzero()[0].shape[0] # the rank of A
        inv_S = np.matrix(np.diag([1.0/s[i] for i in xrange(r)]))
        inv_A = (V[:r, :].T) * inv_S * (U[:, :r].T)
        P_R = U[:,:r] * (U[:,:r].T) # a projection matrix onto the row-space of A
        P_N = U[:,r:] * (U[:,r:].T) # a projection matrix onto the null-space of A
        return inv_A, r, P_R, P_N
        
    def train(self, S, G):
        """
            Estimate standard Gibbs energies of formation
        """
        b = np.matrix(self.training_data.dG0).T
        w = np.matrix(self.training_data.weight).T
        
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

        # Put all the calculated data in 'params' for the sake of debugging
        params = {}
        params['contributions'] = [dG0_rc, dG0_gc]
        params['covariances']   = [V_rc, V_gc, V_inf]
        params['MSEs']          = [MSE_rc, MSE_gc, MSE_inf]
        params['projections']   = [P_R_rc,
                                   P_R_gc * G.T * P_N_rc,
                                   P_N_gc * G.T * P_N_rc,
                                   P_R_gc,
                                   P_N_rc,
                                   P_N_gc]
        params['inverses']      = [inv_S, inv_GS, inv_SWS, inv_GSWGS]

        # Calculate the total of the contributions and covariances
        cov_dG0 = V_rc * MSE_rc + V_gc * MSE_gc + V_inf * MSE_inf
        
        return dG0_cc, cov_dG0, params

if __name__ == '__main__':
    from kegg_model import KeggModel
    
    reaction_strings = sys.stdin.readlines()
    td = TrainingData()
    cc = ComponentContribution(td)
    model = KeggModel.parse_kegg_model(reaction_strings)
    model_dG0, model_cov_dG0 = cc.estimate_kegg_model(model)
    sys.stdout.write('[[' + ', '.join([str(x) for x in model_dG0.flat]) + ']]')