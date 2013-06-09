import re
import numpy as np
from kegg_reaction import KeggReaction
from compound_cacher import CompoundCacher

class KeggModel(object):
    
    def __del__(self):
        self.ccache.dump()
    
    def __init__(self, S, cids):
        self.S = S
        self.cids = cids
        assert len(self.cids) == self.S.shape[0]
        self.ccache = CompoundCacher.getInstance()
    
    @staticmethod
    def from_file(fname, arrow='<=>', has_reaction_ids=False):
        """
        reads a file containing reactions in KEGG format
        
        Arguments:
           fname            - the filename to read
           arrow            - the string used as the 'arrow' in each reaction (default: '<=>')
           has_reaction_ids - a boolean flag indicating if there is a column of
                              reaction IDs (separated from the reaction with
                              whitespaces)
        
        Return a KeggModel
        """
        fd = open(fname, 'r')
        reaction_strings = fd.readlines()
        fd.close()
        return KeggModel.from_formulas(reaction_strings, arrow, has_reaction_ids)
        
    @staticmethod
    def from_formulas(reaction_strings, arrow='<=>', has_reaction_ids=False):
        """
        parses a list of reactions in KEGG format
        
        Arguments:
           reaction_strings - a list of reactions in KEGG format
           arrow            - the string used as the 'arrow' in each reaction (default: '<=>')
           has_reaction_ids - a boolean flag indicating if there is a column of
                              reaction IDs (separated from the reaction with
                              whitespaces)
        
        Return values:
           S     - a stoichiometric matrix
           cids  - the KEGG compound IDs in the same order as the rows of S
        """
        
        cids = set()
        reactions = []
        for line in reaction_strings:
            if has_reaction_ids:
                tokens = re.split('(\w+)\s+(.*)', line, maxsplit=1)
                line = tokens[1]
            reaction = KeggReaction.parse_formula(line, arrow)
            if not reaction.is_balanced():
                raise ValueError('Model contains unbalanced reactions')
            cids = cids.union(reaction.keys())
            reactions.append(reaction)
        
        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids'.
        cids = sorted(cids)
        S = np.zeros((len(cids), len(reactions)))
        for i, reaction in enumerate(reactions):
            for cid, coeff in reaction.iteritems():
                S[cids.index(cid), i] = coeff
                
        return KeggModel(S, cids)
            
    def get_transform_ddG0(self, pH, I, T):
        """
        needed in order to calculate the transformed Gibbs energies of the 
        model reactions.
        
        Returns:
            an array (whose length is self.S.shape[1]) with the differences
            between DrG0_prime and DrG0. Therefore, one must add this array
            to the chemical Gibbs energies of reaction (DrG0) to get the 
            transformed values
        """
        ddG0_compounds = np.zeros((self.S.shape[0], 1))
        for i, cid in enumerate(self.cids):
            comp = self.ccache.get_kegg_compound(cid)
            ddG0_compounds[i, 0] = comp.transform(pH, I, T)
        
        ddG0_forward = np.dot(self.S.T, ddG0_compounds)
        return ddG0_forward
        
    def check_S_balance(self):
        elements, Ematrix = self.ccache.get_kegg_ematrix(self.cids)
        conserved = Ematrix.T * self.S
        rxnFil = np.any(conserved[:,range(self.S.shape[1])],axis=0)
        "set all stoi coefficients in unbalanced reactions to 0"
        self.S[:,np.nonzero(rxnFil)[1]] = 0
        return self
