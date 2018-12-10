# Script for generating the group decompositions of all KEGG compounds
# in a Matlab-friendly format

import sys
sys.path.append('src')

import argparse
import gzip
import json
import pandas as pd
from component_contribution.inchi2gv import (
    init_groups_data, InChIDecomposer, GroupDecompositionError)


def MakeOpts():
    parser = argparse.ArgumentParser(
            description='Decompose all KEGG compounds into group vectors')
    parser.add_argument('-v', '--verbose', action="store_true",
                        help="print warnings to stderr")
    parser.add_argument('output', type=str,
                        help="path to output CSV file")
    return parser

if __name__ == '__main__':
    parser = MakeOpts()
    args = parser.parse_args()

    groups_data = init_groups_data()
    decomposer = InChIDecomposer(groups_data)
    group_names = groups_data.GetGroupNames()
    
    data = []
    with gzip.open('src/component_contribution/data/kegg_compounds.json.gz') as fp:
        for d in json.load(fp):
            try:
                gv = decomposer.inchi_to_groupvec(d['inchi'])
                error_msg = None
            except GroupDecompositionError as e:
                gv = [None] * len(group_names)
                error_msg = str(e)
                if args.verbose:
                    sys.stderr.write(d['compound_id'] + ': ' + error_msg + '\n')

            data.append([d['compound_id'], d['inchi']] + gv + [error_msg])
    
    cols = ['kegg_id', 'inchi'] + group_names + ['error']
    df = pd.DataFrame(data=data, columns=cols).set_index('kegg_id')
    sys.stderr.write('Writing results to: ' + args.output)
    df.to_csv(args.output)
