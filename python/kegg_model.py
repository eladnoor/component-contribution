import re, csv, logging
import numpy as np
from kegg_reaction import KeggReaction
from compound_cacher import CompoundCacher

class KeggModel(object):
    
    def __del__(self):
        self.ccache.dump()
    
    def __init__(self, S, cids, rids=None):
        self.S = S
        self.cids = cids
        self.rids = rids
        assert len(self.cids) == self.S.shape[0]
        if self.rids is not None:
            assert len(self.rids) == self.S.shape[1]
        self.ccache = CompoundCacher.getInstance()
    
    @staticmethod
    def from_file(fname, arrow='<=>', format='kegg', has_reaction_ids=False):
        """
        reads a file containing reactions in KEGG format
        
        Arguments:
           fname            - the filename to read
           arrow            - the string used as the 'arrow' in each reaction (default: '<=>')
           format           - the text file format provided ('kegg', 'tsv' or 'csv')
           has_reaction_ids - a boolean flag indicating if there is a column of
                              reaction IDs (separated from the reaction with
                              whitespaces)
        
        Return a KeggModel
        """
        fd = open(fname, 'r')
        if format == 'kegg':
            model = KeggModel.from_formulas(fd.readlines(), arrow, has_reaction_ids)
        elif format == 'tsv':
            model = KeggModel.from_csv(fd, has_reaction_ids=has_reaction_ids, delimiter='\t')
        elif format == 'csv':
            model = KeggModel.from_csv(fd, has_reaction_ids=has_reaction_ids, delimiter=None)
        fd.close()
        return model
    
    @staticmethod
    def from_csv(fd, has_reaction_ids=True, delimiter=None):
        csv_reader = csv.reader(fd, delimiter=delimiter)
        if has_reaction_ids:
            rids = csv_reader.next()
            rids = rids[1:]
        else:
            rids = None
        S = []
        cids = []
        for i, row in enumerate(csv_reader):
            cids.append(int(row[0]))
            S.append([float(x) for x in row[1:]])
        S = np.array(S)
        return KeggModel(S, cids, rids)
    
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
        if has_reaction_ids:
            rids = []
        else:
            rids = None
        reactions = []
        for line in reaction_strings:
            if has_reaction_ids:
                tokens = re.split('(\w+)\s+(.*)', line, maxsplit=1)
                rids.append(tokens[0])
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
                
        return KeggModel(S, cids, rids)

    def add_thermo(self, cc):
        dG0_cc, cov_dG0 = cc.estimate_kegg_model(self.S, self.cids)
        self.dG0 = dG0_cc
        self.cov_dG0 = cov_dG0
        
    def get_transformed_dG0(self, pH, I, T):
        """
            returns the estimated dG0_prime and the standard deviation of
            each estimate (i.e. a measure for the uncertainty).
        """
        dG0_prime = self.dG0 + self._get_transform_ddG0(pH=pH, I=I, T=T)
        dG0_std = np.matrix(np.sqrt(np.diag(self.cov_dG0))).T
        
        print np.nonzero(dG0_std < 1e-5)[0]
        
        return dG0_prime, dG0_std

    def _get_transform_ddG0(self, pH, I, T):
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
        unbalanced_ind = np.nonzero(rxnFil)[1]
        if unbalanced_ind != []:
            logging.warning('There are (%d) unbalanced reactions in S. ' 
                            'Setting their coefficients to 0.' % 
                            len(unbalanced_ind.flat))
            if self.rids is not None:
                logging.warning('These are the unbalanced reactions: ' +
                                ', '.join([self.rids[i] for i in unbalanced_ind.flat]))
                    
            self.S[:, unbalanced_ind] = 0
        return self
