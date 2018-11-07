# The MIT License (MIT)
#
# Copyright (c) 2013 The Weizmann Institute of Science.
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from __future__ import absolute_import

import logging
import json
from six import string_types
import pandas as pd
from collections import defaultdict

from component_contribution.compound import Compound
from component_contribution.singleton import Singleton

import os
base_path = os.path.split(os.path.realpath(__file__))[0]
DEFAULT_CACHE_FNAME = os.path.join(base_path, '../cache/compounds.csv') 

class CompoundCache(Singleton):
    """
    CompoundCache is a singleton that handles caching of Compound objects for
    the component-contribution package.  The Compounds are retrieved by their ID
    (e.g., KEGG COMPOUND ID, ChEBI Id or HMDB in most cases) or InChI Key.  The
    first time a Compound is requested, it is obtained from the relevant
    database and a Compound object is created (this takes a while because it
    usually involves internet communication and then invoking the ChemAxon
    plugin for calculating the pKa values for that structure). Any further
    request for the same Compound ID will draw the object from the cache. When
    the method dump() is called, all cached data is written to a file that will
    be loaded in future python sessions.
    """

    COLUMNS = [
        "inchi", "name", "atom_bag", "p_kas", "smiles", "major_ms",
        "number_of_protons", "charges"]
    
    SERIALIZE_COLUMNS = ["atom_bag", "p_kas", "number_of_protons", "charges"]

    def __init__(self, cache_file_name=None):
        # a lookup table for compound IDs
        self._compound_id_to_inchi_key = dict()
        
        # a lookup table for InChIKeys
        self._inchi_key_to_compound_ids = defaultdict(set)
        
        # an internal cache for Compound objects that have already been created
        self.compound_dict = {} 

        if not os.path.exists(DEFAULT_CACHE_FNAME):
            self._data = pd.DataFrame(columns=self.COLUMNS)
            self._data.index.name = 'inchi_key'
            self.to_data_frame.to_csv(DEFAULT_CACHE_FNAME)
        else:
            self._data = pd.read_csv(DEFAULT_CACHE_FNAME, index_col=0)
            
            for col in self.SERIALIZE_COLUMNS:
                self._data[col] = self._data[col].apply(json.loads)
            
            self._read_cross_refs(self._data)
            self._data.drop('cross_references', axis=1, inplace=True)
            self._data['inchi'].fillna('', inplace=True)
            self._data['smiles'].fillna('', inplace=True)

    def _read_cross_refs(self, data):
        for inchi_key, row in data.iterrows():
            for compound_id in row.cross_references.split(";"):
                self._compound_id_to_inchi_key[compound_id] = inchi_key
                self._inchi_key_to_compound_ids[inchi_key].add(compound_id)

    def all_compound_ids(self):
        return sorted(self._compound_id_to_inchi_key.keys())

    def to_data_frame(self):
        """
        Writes the cache to a file.

        Parameters
        ----------
        file_name : str
            The file name.

        """
        df = pd.DataFrame(self._data)

        df['cross_references'] = \
            df.index.map(self._inchi_key_to_compound_ids.get)
        df['cross_references'] = df['cross_references'].apply(
            ';'.join)

        for col in self.SERIALIZE_COLUMNS:
            df[col].apply(json.dumps, inplace=True)

        return df

    def get_compound(self, compound_id, compute_pkas=True):
        if compound_id in self._compound_id_to_inchi_key:  # compound exists
            logging.debug('Cache hit for %s' % compound_id)
            return self.get(self._compound_id_to_inchi_key[compound_id])
        # compound_id is an InChI Key and exists
        elif compound_id in self._data.index:
            logging.debug('Cache hit for InChiKey %s' % compound_id)
            return self.get(compound_id)
        else:  # compound does not exist.
            logging.debug('Cache miss, calculating pKas for %s' % compound_id)
            cpd = Compound.get(compound_id, compute_pkas)
            if cpd.inchi_key in self._data.index:
                cpd = self.get(cpd.inchi_key)
                self._compound_id_to_inchi_key[compound_id] = cpd.inchi_key
            else:
                self.add(cpd)

            return cpd

    def remove(self, inchi_key):
        if inchi_key not in self._data.index:
            raise KeyError(inchi_key)

        del self.compound_dict[inchi_key]
        del self._inchi_key_to_compound_ids
        self._data.drop(inchi_key)

        for compound_id in self._compound_id_to_inchi_key:
            if self._compound_id_to_inchi_key[compound_id] == inchi_key:
                del self._compound_id_to_inchi_key[compound_id]

    def get(self, inchi_key):
        if inchi_key not in self._data.index:
            raise KeyError(inchi_key)

        if inchi_key in self.compound_dict:
            cpd = self.compound_dict[inchi_key]
        else:
            data = self._data.loc[inchi_key]
            cpd = Compound(
                inchi_key, data.inchi, data.atom_bag, data.p_kas, data.smiles,
                data.major_ms, data.number_of_protons, data.charges)
            self.compound_dict[inchi_key] = cpd

        return cpd

    def add(self, cpd):
        if cpd.inchi_key not in self._inchi_key_to_compound_ids:
            self._inchi_key_to_compound_ids[cpd.inchi_key] = set()

        self._inchi_key_to_compound_ids[cpd.inchi_key].add(cpd.compound_id)
        self._compound_id_to_inchi_key[cpd.compound_id] = cpd.inchi_key
        self.compound_dict[cpd.inchi_key] = cpd
        self._data.loc[cpd.inchi_key] = [
            cpd.inchi, cpd.name, cpd.atom_bag, cpd.p_kas, cpd.smiles,
            cpd.major_microspecies, cpd.number_of_protons, cpd.charges]

    def get_element_data_frame(self, compound_ids):
        if isinstance(compound_ids, string_types):
            compound_ids = [compound_ids]

        # gather the "atom bags" of all compounds in a list 'atom_bag_list'
        atom_bags = [self.get_compound(cid).atom_bag or {}
                     for cid in compound_ids]

        elements = sorted(set([e for bag in atom_bags for e in bag.keys()]))

        # create the elemental matrix, where each row is a compound and each
        # column is an element (or e-)
        element_df = pd.DataFrame(index=compound_ids, columns=elements,
                                  dtype=int).fillna(0)

        for cid, atom_bag in zip(compound_ids, atom_bags):
            for elem in elements:
                element_df[elem][cid] = atom_bag.get(elem, 0)
        
        return element_df

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    ccache = CompoundCache()

    # find out all the KEGG compound IDs that are used in the training data
    # and add them to the cache
    from component_contribution.training_data import TrainingData
    thermo_df, _ = TrainingData.get_all_thermo_params()
    cids = set(['C00001', 'C00080'])
    for rxn in thermo_df['reaction']:
        cids = cids.union(rxn.keys())
    cids = sorted(cids)
    
    for cid in cids:
        ccache.get_compound('KEGG:%s' % cid)
    
    ccache.dump()
    
