import re
import kegg_errors

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
        if len(tokens) == 1:
            amount = 1
            key = member
        else:
            try:
                amount = float(tokens[0])
            except ValueError:
                raise kegg_errors.KeggParseException(
                    "Non-specific reaction: %s" % s)
            key = tokens[1]
            
        if key[0] != 'C':
            raise kegg_errors.KeggNonCompoundException(
                "Compound ID doesn't start with C: %s" % key)
        try:
            cid = int(key[1:])
            compound_bag[cid] = compound_bag.get(cid, 0) + amount
        except ValueError:
            raise kegg_errors.KeggParseException(
                "Non-specific reaction: %s" % s)
    
    return compound_bag

def parse_kegg_reaction(formula, arrow='<=>'):
    """ 
        Parses a two-sided formula such as: 2 C00001 => C00002 + C00003 
        
        Return:
            The set of substrates, products and the direction of the reaction
    """
    tokens = formula.split(arrow)
    if len(tokens) < 2:
        raise ValueError('Reaction does not contain the arrow sign: ' + formula)
    if len(tokens) > 2:
        raise ValueError('Reaction contains more than one arrow sign: ' + formula)
    
    left = tokens[0].strip()
    right = tokens[1].strip()
    
    sparse_reaction = {}
    for cid, count in parse_reaction_formula_side(left).iteritems():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 

    for cid, count in parse_reaction_formula_side(right).iteritems():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 

    return sparse_reaction