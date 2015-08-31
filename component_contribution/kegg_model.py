import re, csv, logging
import numpy as np
from kegg_reaction import KeggReaction
from compound_cacher import CompoundCacher
from kegg_errors import KeggParseException

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
        self.ccache = CompoundCacher()

        # remove H+ from the stoichiometric matrix if it exists
        if 'C00080' in self.cids:
            i = self.cids.index('C00080')
            self.S = np.vstack((self.S[:i,:], self.S[i+1:,:]))
            self.cids.pop(i)

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
            cids.append(row[0])
            S.append([float(x) for x in row[1:]])
        S = np.array(S)

        return KeggModel(S, cids, rids)
    
    @staticmethod
    def from_kegg_reactions(kegg_reactions, has_reaction_ids=False):
        if has_reaction_ids:
            rids = [r.rid for r in kegg_reactions]
        else:
            rids = None

        cids = set()
        for reaction in kegg_reactions:
            cids = cids.union(reaction.keys())
        
        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids'.
        cids = sorted(cids)
        S = np.matrix(np.zeros((len(cids), len(kegg_reactions))))
        for i, reaction in enumerate(kegg_reactions):
            S[:, i] = np.matrix(reaction.dense(cids))
        
        logging.debug('Successfully loaded %d reactions (involving %d unique compounds)' %
                      (S.shape[1], S.shape[0]))
        return KeggModel(S, cids, rids)
    
    @staticmethod
    def from_formulas(reaction_strings, arrow='<=>', has_reaction_ids=False,
                      raise_exception=False):
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
        try:
            reactions = []
            not_balanced_count = 0
            for line in reaction_strings:
                rid = None
                if has_reaction_ids:
                    tokens = re.findall('(\w+)\s+(.*)', line.strip())[0]
                    rid = tokens[0]
                    line = tokens[1]
                try:
                    reaction = KeggReaction.parse_formula(line, arrow, rid)
                except KeggParseException as e:
                    logging.warning(str(e))
                    reaction = KeggReaction({})
                if not reaction.is_balanced(fix_water=True, raise_exception=raise_exception):
                    not_balanced_count += 1
                    logging.warning('Model contains an unbalanced reaction: ' + line)
                    reaction = KeggReaction({})
                reactions.append(reaction)
                logging.debug('Adding reaction: ' + reaction.write_formula())
            
            if not_balanced_count > 0:
                warning_str = '%d out of the %d reactions are not chemically balanced' % \
                              (not_balanced_count, len(reaction_strings))
                logging.debug(warning_str)
            return KeggModel.from_kegg_reactions(reactions, has_reaction_ids)
        
        except ValueError as e:
            if raise_exception:
                raise e
            else:
                logging.debug(str(e))
                return None

    def add_thermo(self, cc):
        # check that all CIDs in the reaction are already cached by CC
        Nc, Nr = self.S.shape
        reactions = []
        for j in xrange(Nr):
            sparse = {self.cids[i]:self.S[i,j] for i in xrange(Nc)
                      if self.S[i,j] != 0}
            reaction = KeggReaction(sparse)
            reactions.append(reaction)
            
        self.dG0, self.cov_dG0 = cc.get_dG0_r_multi(reactions)
        
    def get_transformed_dG0(self, pH, I, T):
        """
            returns the estimated dG0_prime and the standard deviation of
            each estimate (i.e. a measure for the uncertainty).
        """
        dG0_prime = self.dG0 + self._get_transform_ddG0(pH=pH, I=I, T=T)
        dG0_std = np.matrix(np.sqrt(np.diag(self.cov_dG0))).T
        U, s, V = np.linalg.svd(self.cov_dG0, full_matrices=True)
        sqrt_Sigma = np.matrix(U) * np.matrix(np.diag(s**0.5)) * np.matrix(V)
        return dG0_prime, dG0_std, sqrt_Sigma

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
        ddG0_compounds = np.matrix(np.zeros((self.S.shape[0], 1)))
        for i, cid in enumerate(self.cids):
            comp = self.ccache.get_compound(cid)
            ddG0_compounds[i, 0] = comp.transform_pH7(pH, I, T)
        
        ddG0_forward = np.dot(self.S.T, ddG0_compounds)
        return ddG0_forward
        
    def check_S_balance(self, fix_water=False):
        elements, Ematrix = self.ccache.get_element_matrix(self.cids)
        conserved = Ematrix.T * self.S

        if fix_water:
            # This part only looks for imbalanced oxygen and uses extra
            # H2O molecules (on either side of the reaction equation) to
            # balance them. Keep in mind that also the e- balance is affected
            # by the water (and hydrogen is not counted at all).
            if 'C00001' not in self.cids:
                self.S = np.vstack([self.S, np.zeros((1, self.S.shape[1]))])
                self.cids.append('C00001')
                elements, Ematrix = self.ccache.get_element_matrix(self.cids)
            
            i_h2o = self.cids.index('C00001')
            add_water = -conserved[elements.index('O'), :]
            self.S[i_h2o, :] += add_water
            conserved += Ematrix[i_h2o, :].T * add_water

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

    def write_reaction_by_index(self, r):
        sparse = dict([(cid, self.S[i, r]) for i, cid in enumerate(self.cids)
                       if self.S[i, r] != 0])
        if self.rids is not None:
            reaction = KeggReaction(sparse, rid=self.rids[r])
        else:
            reaction = KeggReaction(sparse)
        return reaction.write_formula()
        
    def get_unidirectional_S(self):
        S_plus = np.copy(self.S)
        S_minus = np.copy(self.S)
        S_plus[self.S < 0] = 0
        S_minus[self.S > 0] = 0
        return S_minus, S_plus
        