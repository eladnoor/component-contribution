import re
import numpy as np
from kegg_reaction import parse_kegg_reaction
from compound_cacher import CompoundCacher

class KeggModel(object):
    
    def __init__(self, S, cids):
        self.S = S
        self.cids = cids
    
    @staticmethod
    def load_kegg_model(fname, arrow='<=>', has_reaction_ids=False):
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
        return KeggModel.parse_kegg_model(reaction_strings, arrow, has_reaction_ids)

    @staticmethod
    def parse_kegg_model(reaction_strings, arrow, has_reaction_ids):
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
            reaction = parse_kegg_reaction(line, arrow)
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
        
