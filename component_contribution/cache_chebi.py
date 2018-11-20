# The MIT License (MIT)
#
# Copyright (c) 2013 The Weizmann Institute of Science.
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
# Copyright (c) 2018 Institute for Molecular Systems Biology, ETH Zurich
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

import logging
from ftplib import FTP

import pandas as pd
import pybel

from .compound import Compound
from .compound_cache import ccache


logger = logging.getLogger(__name__)

CHEBI_FTP = ('ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited'
             '/chebiId_inchi.tsv')

ftp = FTP('ftp.ebi.ac.uk')
ftp.login()
ftp.cwd('pub/databases/chebi/Flat_file_tab_delimited')

table = []
ftp.retrlines('RETR chebiId_inchi.tsv',
              lambda s: table.append(s.split('\t', 1)))
df = pd.DataFrame(data=table[1:], columns=table[0])

for i, row in df.iterrows():
    compound_id = 'ChEBI:' + row['CHEBI_ID']
    print('%05d) Adding %s to cache' % (i, compound_id))
    if ccache.exists(compound_id):
        print('- already in cache, skipping')
    else:
        molecule = pybel.readstring("inchi", row['InChI'])
        cpd = Compound.from_molecule(compound_id, molecule,
                                     compute_pkas=False)
        if not ccache.exists(cpd.inchi_key):
            print('- new compound, calculating pKas and adding to cache')
            cpd = Compound.from_molecule(compound_id, molecule,
                                         compute_pkas=True)
        else:
            print('- this InChI is already in the cache, adding cross-link')

        ccache.add(cpd)

    if (i+1) % 10 == 0:
        print('Dumping cache')
        ccache.dump()
