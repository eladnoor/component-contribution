# -*- coding: utf-8 -*-
"""
Created on Tue May 31 10:57:02 2016

@author: noore
"""

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
