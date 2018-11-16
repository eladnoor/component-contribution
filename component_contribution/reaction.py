import re
import numpy as np
import logging
from .compound_cache import ccache
from .exceptions import ParseException

class Reaction(object):

    def __init__(self, sparse, arrow='<=>', rid=None):
        for cid, coeff in sparse.items():
            if type(coeff) not in [float, int]:
                raise ValueError('All values in Reaction must be integers or floats')
        
        self.sparse = dict([(k, v) for (k, v) in sparse.items() if v != 0])
        self.arrow = arrow
        self.rid = rid

    def keys(self):
        return self.sparse.keys()
        
    def items(self):
        return self.sparse.items()

    def __str__(self):
        return self.write_formula()

    def reverse(self):
        """
            reverse the direction of the reaction by negating all stoichiometric
            coefficients
        """
        self.sparse = dict( (k, -v) for (k, v) in self.sparse.items() )

    @staticmethod
    def parse_reaction_formula_side(s):
        """ 
            Parses the side formula, e.g. '2 KEGG:C00001 + KEGG:C00002 + 3 KEGG:C00003'
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
                    raise ParseException(
                        "Non-specific reaction: %s" % s)
                key = tokens[1]
                
            try:
                compound_bag[key] = compound_bag.get(key, 0) + amount
            except ValueError:
                raise ParseException(
                    "Non-specific reaction: %s" % s)
        
        return compound_bag

    @staticmethod
    def parse_formula(formula, arrow='<=>', rid=None):
        """ 
            Parses a two-sided formula such as: 2 KEGG:C00001 => KEGG:C00002 + KEGG:C00003 
            
            Return:
                The set of substrates, products and the direction of the reaction
        """
        tokens = formula.split(arrow)
        if len(tokens) < 2:
            raise ParseException('Reaction does not contain the arrow sign (%s): %s'
                                     % (arrow, formula))
        if len(tokens) > 2:
            raise ParseException('Reaction contains more than one arrow sign (%s): %s'
                                     % (arrow, formula))
        
        left = tokens[0].strip()
        right = tokens[1].strip()
        
        sparse_reaction = {}
        for cid, count in Reaction.parse_reaction_formula_side(left).items():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 

        for cid, count in Reaction.parse_reaction_formula_side(right).items():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 

        return Reaction(sparse_reaction, arrow, rid=rid)

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
        for cid, coeff in sorted(self.sparse.items()):
            if coeff < 0:
                left.append(Reaction.write_compound_and_coeff(cid, -coeff))
            elif coeff > 0:
                right.append(Reaction.write_compound_and_coeff(cid, coeff))
        return "%s %s %s" % (' + '.join(left), self.arrow, ' + '.join(right))

    def _get_reaction_atom_bag(self, raise_exception=False):
        """
            Use for checking if all elements are conserved.
            
            Returns:
                An atom_bag of the differences between the sides of the reaction.
                E.g. if there is one extra C on the left-hand side, the result will
                be {'C': -1}.
        """
        cids = list(self.keys())
        coeffs = np.array(list(map(self.sparse.__getitem__, cids)))

        element_df = ccache.get_element_data_frame(cids)
        if element_df.shape[1] == 0 or np.any( (element_df == 0).all(axis=1) ):
            warning_str = 'cannot generate the reaction atom bag because ' + \
                          'compounds have unspecific formulas: ' + \
                          '%s' % self.write_formula()
            if raise_exception:
                raise ValueError(warning_str)
            else:
                logging.warning(warning_str)
                return None

        conserved = coeffs @ element_df.values

        atom_bag = {}        
        if np.any(conserved != 0):
            logging.debug('unbalanced reaction: %s' % self.write_formula())
            for j, c in enumerate(conserved.flat):
                if c != 0:
                    logging.debug('there are %d more %s atoms on the right-hand side' %
                                  (c, element_df.columns[j]))
                    atom_bag[str(element_df.columns[j])] = c
        return atom_bag
            
    def is_balanced(self, fix_protons=True, fix_water=False,
                    raise_exception=False):
        reaction_atom_bag = self._get_reaction_atom_bag(
                raise_exception=raise_exception)

        if reaction_atom_bag is None:
            # this means some compound formulas are missing
            return False

        if fix_water and 'O' in reaction_atom_bag:
            self.sparse.setdefault('KEGG:C00001', 0)
            self.sparse['KEGG:C00001'] += -reaction_atom_bag['O']
            if self.sparse['KEGG:C00001'] == 0:
                del self.sparse['KEGG:C00001']
            reaction_atom_bag = self._get_reaction_atom_bag()

        if fix_protons and 'H' in reaction_atom_bag:
            self.sparse.setdefault('KEGG:C00080', 0)
            self.sparse['KEGG:C00080'] += -reaction_atom_bag['H']
            if self.sparse['KEGG:C00080'] == 0:
                del self.sparse['KEGG:C00080']
            reaction_atom_bag = self._get_reaction_atom_bag()

        return len(reaction_atom_bag) == 0

    def is_empty(self):
        return len(self.sparse) == 0
            
    def dense(self, cids):
        s = np.zeros((len(cids), 1))
        for cid, coeff in self.items():
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
        for compound_id, coeff in self.items():
            if compound_id in ['C00080', 'KEGG:C00080']:
                continue # H+ is ignored in the Legendre transform
            comp = ccache.get_compound(compound_id)
            ddG0_forward += coeff * comp.transform_p_h_7(pH, I, T)
        return ddG0_forward
        
if __name__ == '__main__':
    reaction = Reaction.parse_formula('KEGG:C00149 <=> KEGG:C00036')
    print(reaction)
    print('Is balanced: ', reaction.is_balanced())
    print(reaction._get_reaction_atom_bag())

    reaction = Reaction.parse_formula('KEGG:C00149 + KEGG:C00003 <=> KEGG:C00036 + KEGG:C00004')
    print(reaction)
    print('Is balanced: ', reaction.is_balanced())
    print(reaction._get_reaction_atom_bag())

    reaction = Reaction.parse_formula('KEGG:C00149 + KEGG:C00138 <=> KEGG:C00139 + KEGG:C00036')
    print(reaction)
    try:
        reaction.is_balanced()
    except ValueError as e:
        print(str(e))
        
