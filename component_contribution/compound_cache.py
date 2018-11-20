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


import atexit
import json
import logging
import sqlite3
from collections import defaultdict
from contextlib import ExitStack

import pandas as pd
import pandas.io.sql as pd_sql
from importlib_resources import path

import component_contribution.cache
from component_contribution.compound import Compound


logger = logging.getLogger(__name__)
file_manager = ExitStack()
atexit.register(file_manager.close)

DEFAULT_CACHE_FNAME = str(file_manager.enter_context(
    path(component_contribution.cache, "compounds.sqlite")))


class CompoundCache(object):
    """
    CompoundCache is a singleton that handles caching of Compound objects for
    the component-contribution package.  The Compounds are retrieved by their
    ID (e.g., KEGG COMPOUND ID, ChEBI Id or HMDB in most cases) or InChI Key.
    The first time a Compound is requested, it is obtained from the relevant
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

    def __init__(self, cache_file_name=DEFAULT_CACHE_FNAME):
        # connect to the SQL database
        self.con = sqlite3.connect(cache_file_name)

        # a lookup table for compound IDs
        self._compound_id_to_inchi_key = dict()

        # a lookup table for InChIKeys
        self._inchi_key_to_compound_ids = defaultdict(set)

        # an internal cache for Compound objects that have already been created
        self.compound_dict = {}

        # a flag telling the Cache that it needs to rewrite itself in the
        # filesystem. updated any time a new compound is cached.
        self.requires_update = False

        if pd_sql.has_table('compound_cache', self.con):
            df = pd.read_sql('SELECT * FROM compound_cache', self.con)
            df.set_index('inchi_key', inplace=True, drop=True)
            self.from_data_frame(df)
        else:
            self._data = pd.DataFrame(columns=self.COLUMNS)
            self._data.index.name = 'inchi_key'
            self.requires_update = True
            self.dump()

    def dump(self, cache_file_name=DEFAULT_CACHE_FNAME):
        if self.requires_update:
            df = self.to_data_frame()
            pd_sql.to_sql(df, 'compound_cache', self.con, if_exists='replace')
            self.requires_update = False
            self.con.commit()

    def _read_cross_refs(self, df):
        for inchi_key, row in df.iterrows():
            for compound_id in row.cross_references.split(";"):
                self._compound_id_to_inchi_key[compound_id] = inchi_key
                self._inchi_key_to_compound_ids[inchi_key].add(compound_id)

    def all_compound_ids(self):
        return sorted(self._compound_id_to_inchi_key.keys())

    def to_data_frame(self):
        """
            Creates a DataFrame that can be written to the SQLite database.
            Some of the columns need to be serialized before the export.

            Returns:
                pandas.DataFrame
        """
        df = self._data.copy()

        df['cross_references'] = \
            df.index.map(self._inchi_key_to_compound_ids.get)
        df['cross_references'] = df['cross_references'].apply(
            ';'.join)

        for col in self.SERIALIZE_COLUMNS:
            df[col] = df[col].apply(json.dumps)

        return df

    def from_data_frame(self, df):
        """
            Reads all the cached data from a DataFrame

            Arguments:
                df - a DataFrame created earlier by to_data_frame()
        """
        self._data = df.copy()

        for col in self.SERIALIZE_COLUMNS:
            self._data[col] = self._data[col].apply(json.loads)

        self._read_cross_refs(df)
        self._data.drop('cross_references', axis=1, inplace=True)
        self._data['major_ms'] = self._data['major_ms'].apply(int)
        self._data['inchi'].fillna('', inplace=True)
        self._data['smiles'].fillna('', inplace=True)

    def exists(self, compound_id):
        return ((compound_id in self._compound_id_to_inchi_key) or
                (compound_id in self._data.index))

    def add_cross_link(self, inchi_key, new_compound_id):
        logging.debug('Adding a cross-link from %s to %s' %
                      (new_compound_id, inchi_key))
        cpd = self.get(inchi_key)
        self._compound_id_to_inchi_key[new_compound_id] = cpd.inchi_key
        self.requires_update = True

    def get_compound(self, compound_id, compute_pkas=True):
        # compound exists
        if compound_id in self._compound_id_to_inchi_key:
            logging.debug('Cache hit for %s' % compound_id)
            return self.get(self._compound_id_to_inchi_key[compound_id])
        # compound_id is an InChI Key and exists
        elif compound_id in self._data.index:
            logging.debug('Cache hit for InChiKey %s' % compound_id)
            return self.get(compound_id)
        # compound does not exist.
        else:
            logging.debug('Cache miss, calculating pKas for %s' % compound_id)
            cpd = Compound.get(compound_id, compute_pkas)
            if cpd.inchi_key in self._data.index:
                self.add_cross_link(cpd.inchi_key, compound_id)
            else:
                logging.debug('Adding the new InChiKey to the cache: %s'
                              % cpd.inchi_key)
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
            data = self._data.loc[inchi_key, :]
            cpd = Compound(
                inchi_key, data.inchi, data.atom_bag, data.p_kas, data.smiles,
                int(data.major_ms), data.number_of_protons, data.charges)
            self.compound_dict[inchi_key] = cpd

        return cpd

    def add(self, cpd):
        if cpd.inchi_key not in self._inchi_key_to_compound_ids:
            self._inchi_key_to_compound_ids[cpd.inchi_key] = set()

        self._inchi_key_to_compound_ids[cpd.inchi_key].add(cpd.compound_id)
        self._compound_id_to_inchi_key[cpd.compound_id] = cpd.inchi_key
        self.compound_dict[cpd.inchi_key] = cpd
        self._data.loc[cpd.inchi_key, :] = [
            cpd.inchi, cpd.name, cpd.atom_bag, cpd.p_kas, cpd.smiles,
            int(cpd.major_microspecies), cpd.number_of_protons, cpd.charges]
        self.requires_update = True

    def get_element_data_frame(self, compound_ids):
        if isinstance(compound_ids, str):
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


# this is the only place where one should use the constructore.
# we wish to only have one instance of the cache (i.e. use it as a singleton)

# legacy code for transitioning from CSV to SQLite
if False:
    con = sqlite3.connect(DEFAULT_CACHE_FNAME)
    df = pd.read_csv('component_contribution/cache/compounds.csv',
                     index_col=0)
    df.index.name = 'inchi_key'
    df.to_sql('compound_cache', con, if_exists='replace')
    con.commit()
    con.close()

    con = sqlite3.connect(DEFAULT_CACHE_FNAME)
    df = pd.read_sql('SELECT * FROM compound_cache', con)
    df.set_index('inchi_key', inplace=True)
    print(df)

    con.close()

ccache = CompoundCache()


if __name__ == '__main__':
    from .training_data import FullTrainingData
    logger.setLevel(logging.DEBUG)

    # Cache all the compounds that are part of the training data.
    # Calling the constructor here,
    # already invokes queries for all the relevant compounds (since their
    # structure is needed in order to balance the reactions and the pKas
    # are needed in order to perform the reverse Legendre transform).
    td = FullTrainingData()

    ccache.dump()
