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

import bioservices.kegg
import pandas as pd

kegg = bioservices.kegg.KEGG()
cid2name = kegg.list('cpd')
cid2name = filter(lambda x: len(x) == 2, map(lambda l : l.split('\t'), cid2name.split('\n')))
cid_df = pd.DataFrame(cid2name, columns=['cpd', 'names'])
cid_df['cpd'] = cid_df['cpd'].apply(lambda x: x[4:])
cid_df['name'] = cid_df['names'].apply(lambda s: s.split(';')[0])
cid_df.set_index('cpd', inplace=True)
cid_df['inchi'] = None

for cid in cid_df.index[0:10]:
    ChEBI = re.findall('ChEBI: ([\d\s]+)\n', kegg.get(cid))
    if len(ChEBI) == 0:
        print 'Cannot find a ChEBI for %s' % cid
    elif len(ChEBI) > 1:
        print 'Error parsing compound %s' % cid
    else:
        cid2chebi.at[cid, 'ChEBI'] = ChEBI[0]

cid2chebi.to_csv(settings.KEGG2CHEBI_FNAME)
