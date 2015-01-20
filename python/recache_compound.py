# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 21:00:31 2014

@author: eladn
"""
import sys
from python.compound_cacher import CompoundCacher

compound_id = sys.argv[1]
CompoundCacher.RebuildCompoundJSON()
ccache = CompoundCacher()
ccache.remove(compound_id)
comp = ccache.get_compound(compound_id)
ccache.dump()
