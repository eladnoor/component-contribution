import numpy as np
from inchi2gv import GroupsData, InChI2GroupVector, GROUP_CSV
from training_data import TrainingData
from kegg_model import KeggModel

class ComponentContribution(object):
    
    def __init__(self, kegg_model, training_data):
        """
            Initialize G matrix, and then use the python script "inchi2gv.py" to decompose each of the 
            compounds that has an InChI and save the decomposition as a row in the G matrix.
        """
        self.groups_data = GroupsData.FromGroupsFile(GROUP_CSV, transformed=False)
        self.inchi2gv = InChI2GroupVector(self.groups_data)
        self.group_names = self.groups_data.GetGroupNames()

        self.kegg_model = kegg_model
        self.training_data = training_data
        ComponentContribution._standardize_models(self.kegg_model, self.training_data)
        
        # create the group incidence matrix
        self.G = np.zeros((len(self.training_data.cids), len(self.group_names)))
        self.has_gv = np.ones((len(self.training_data.cids), 1))
        
        # decompose the compounds in the training_data and add to G
        for i, cid in enumerate(self.training_data.cids):
            # for compounds that have no InChI, add a unique 1 in a new column
            inchi = self.training_data.cid2compound[cid].inchi
            group_def = self.inchi2gv.EstimateInChI(inchi)
            #if length(group_def) == length(Groups) % decompoisition successful, place it in G
            #    G(i, 1:length(group_def)) = group_def;
            #    training_data.has_gv(i) = true;
            #elseif isempty(group_def) % decomposition failed, add a unique 1 in a new column
            #    G(:, end+1) = 0; 
            #    G(i, end) = 1;
            #    training_data.has_gv(i) = false;
            #else
            #    fprintf('InChI = %s\n', inchi);
            #    fprintf('*************\n%s\n', getGroupVectorFromInchi(inchi, false));
            #    error('ERROR: while trying to decompose compound C%05d', training_data.cids(i));

    @staticmethod 
    def _standardize_models(kegg_model, training_data):
        """
            map between the model and the training data compounds
            and zero-pad the S matrices so that their rows (compounds) will be
            consistent
        """
        all_cids = sorted(set(kegg_model.cids + training_data.cids))
        kegg_model.S = ComponentContribution._zero_pad_S(
                                kegg_model.S, kegg_model.cids, all_cids)
        kegg_model.cids = all_cids
        
        training_data.S = ComponentContribution._zero_pad_S(
                                training_data.S, training_data.cids, all_cids)
        training_data.cids = all_cids
    
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
        
    #def create_group_incidence_matrix(self):
        
if __name__ == '__main__':
    td = TrainingData()
    model = KeggModel.load_kegg_model('../examples/wolf_reactions.txt')
    cc = ComponentContribution(model, td)