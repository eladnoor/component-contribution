# The MIT License (MIT)
#
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

import pybel
from requests import get, exceptions

class DatabaseInterface():

    def __init__(self):
        self._registry = {}

    def register(self, database_name, mol_retrieving_function):
        self._registry[database_name] = mol_retrieving_function

    def get_molecule(self, database_name, accession):
        if database_name in self._registry:
            # TODO: We might want to catch `requests.HTTPError` here.
            return self._registry[database_name](accession)
        else:
            raise KeyError("Unknown database: '{}'.".format(database_name))


def get_kegg_molecule(accession):
    try:
        response = get(
            "http://rest.kegg.jp/get/cpd:{}/mol".format(accession))
        response.raise_for_status()
        molstring = str(response.text)
    except exceptions.HTTPError:
        return None
        
    return pybel.readstring("mol", molstring)


def get_hmdb_molecule(accession):
    # TODO: The identifier correction should be moved outside where needed.
    # Ensure that the numeric part of the HMDB identifier is correctly padded.
    response = get(
        "http://www.hmdb.ca/structures/metabolites/HMDB{:0>7}.mol".format(
            accession[4:]))
    response.raise_for_status()
    return pybel.readstring("mol", str(response.text))


def get_inchi_molecule(accession):
    return pybel.readstring("inchi", accession)


def get_chebi_molecule(accession):
    raise NotImplementedError("Coming soon!")


databases = DatabaseInterface()
databases.register("KEGG", get_kegg_molecule)
databases.register("ChEBI", get_chebi_molecule)
databases.register("HMDB", get_hmdb_molecule)
databases.register("InChI", get_inchi_molecule)
