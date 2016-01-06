import re
import numpy as np
import logging
from compound_cacher import CompoundCacher
from kegg_errors import KeggParseException

class KeggReaction(object):

    def __init__(self, sparse, arrow='<=>', rid=None):
        for cid, coeff in sparse.iteritems():
            if not (isinstance(coeff, float) or isinstance(coeff, int)):
                raise ValueError('All values in KeggReaction must be integers or floats')
        self.sparse = dict(filter(lambda (k,v):v, sparse.items()))
        self.arrow = arrow
        self.rid = rid
        self.ccache = CompoundCacher()

    def keys(self):
        return self.sparse.keys()
        
    def iteritems(self):
        return self.sparse.iteritems()

    def __str__(self):
        return self.write_formula()

    def reverse(self):
        """
            reverse the direction of the reaction by negating all stoichiometric
            coefficients
        """
        self.sparse = dict( (k, -v) for (k, v) in self.sparse.iteritems() )

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
                try:
                    amount = float(tokens[0])
                except ValueError:
                    raise KeggParseException(
                        "Non-specific reaction: %s" % s)
                key = tokens[1]
                
            try:
                compound_bag[key] = compound_bag.get(key, 0) + amount
            except ValueError:
                raise KeggParseException(
                    "Non-specific reaction: %s" % s)
        
        return compound_bag

    @staticmethod
    def parse_formula(formula, arrow='<=>', rid=None):
        """ 
            Parses a two-sided formula such as: 2 C00001 => C00002 + C00003 
            
            Return:
                The set of substrates, products and the direction of the reaction
        """
        tokens = formula.split(arrow)
        if len(tokens) < 2:
            raise KeggParseException('Reaction does not contain the arrow sign (%s): %s'
                                     % (arrow, formula))
        if len(tokens) > 2:
            raise KeggParseException('Reaction contains more than one arrow sign (%s): %s'
                                     % (arrow, formula))
        
        left = tokens[0].strip()
        right = tokens[1].strip()
        
        sparse_reaction = {}
        for cid, count in KeggReaction.parse_reaction_formula_side(left).iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 

        for cid, count in KeggReaction.parse_reaction_formula_side(right).iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 

        return KeggReaction(sparse_reaction, arrow, rid=rid)

    @staticmethod
    def write_compound_and_coeff(compound_id, coeff):
        if coeff == 1:
            return compound_id
        else:
            return "%g %s" % (coeff, compound_id)

    def write_formula(self):
        """String representation."""
        left = []
        right = []
        for cid, coeff in sorted(self.sparse.iteritems()):
            if coeff < 0:
                left.append(KeggReaction.write_compound_and_coeff(cid, -coeff))
            elif coeff > 0:
                right.append(KeggReaction.write_compound_and_coeff(cid, coeff))
        return "%s %s %s" % (' + '.join(left), self.arrow, ' + '.join(right))

    def _get_reaction_atom_bag(self, raise_exception=False):
        """
            Use for checking if all elements are conserved.
            
            Returns:
                An atom_bag of the differences between the sides of the reaction.
                E.g. if there is one extra C on the left-hand side, the result will
                be {'C': -1}.
        """
        try:
            cids = list(self.keys())
            coeffs = map(self.sparse.__getitem__, cids)
            coeffs = np.matrix(coeffs)
    
            cached_cids = set(map(str, self.ccache.compound_id2inchi.keys()))
            if not cached_cids.issuperset(cids):
                missing_cids = set(cids).difference(cached_cids)
                warning_str = 'The following compound IDs are not in the cache, ' + \
                              'make sure they appear in kegg_additions.tsv and ' + \
                              'then run compound_cacher.py: ' + \
                              ', '.join(sorted(missing_cids))
                raise ValueError(warning_str)
        
            elements, Ematrix = self.ccache.get_element_matrix(cids)
            conserved = coeffs * Ematrix
    
            if np.any(np.isnan(conserved), 1):
                warning_str = 'cannot test reaction balancing because of unspecific ' + \
                              'compound formulas: %s' % self.write_formula()
                raise ValueError(warning_str)
            
            atom_bag = {}        
            if np.any(conserved != 0, 1):
                logging.debug('unbalanced reaction: %s' % self.write_formula())
                for j, c in enumerate(conserved.flat):
                    if c != 0:
                        logging.debug('there are %d more %s atoms on the right-hand side' %
                                      (c, elements[j]))
                        atom_bag[str(elements[j])] = c
            return atom_bag
            
        except ValueError as e:
            if raise_exception:
                raise e
            else:
                logging.debug(str(e))
                return None

    def is_balanced(self, fix_water=False, raise_exception=False):
        reaction_atom_bag = self._get_reaction_atom_bag(raise_exception)

        if reaction_atom_bag is None: # this means some compound formulas are missing
            return False

        if fix_water and 'O' in reaction_atom_bag:
            self.sparse.setdefault('C00001', 0)
            self.sparse['C00001'] += -reaction_atom_bag['O']
            if self.sparse['C00001'] == 0:
                del self.sparse['C00001']
            reaction_atom_bag = self._get_reaction_atom_bag()

        return len(reaction_atom_bag) == 0

    def is_empty(self):
        return len(self.sparse) == 0
            
    def dense(self, cids):
        s = np.matrix(np.zeros((len(cids), 1)))
        for cid, coeff in self.iteritems():
            s[cids.index(cid), 0] = coeff
        return s

    def get_transform_ddG0(self, pH, I, T):
        """
        needed in order to calculate the transformed Gibbs energies of
        reactions.
        
        Returns:
            The difference between DrG0_prime and DrG0 for this reaction.
            Therefore, this value must be added to the chemical Gibbs
            energy of reaction (DrG0) to get the transformed value.
        """
        ddG0_forward = 0
        for compound_id, coeff in self.iteritems():
            comp = self.ccache.get_compound(compound_id)
            ddG0_forward += coeff * comp.transform_pH7(pH, I, T)
        return ddG0_forward
        
if __name__ == '__main__':
    reaction = KeggReaction.parse_formula('C00149 <=> C00036')
    print reaction._get_reaction_atom_bag()
