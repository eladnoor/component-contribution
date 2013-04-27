#!/usr/bin/python

class KeggParseException(Exception):
    pass
        
class KeggNonCompoundException(Exception):
    pass

class KeggReactionNotBalancedException(Exception):
    def __init__(self, msg, atom_bag=None):
        Exception.__init__(self, msg)
        self.atom_bag = atom_bag or {}
    
class KeggMissingModuleException(Exception):
    pass