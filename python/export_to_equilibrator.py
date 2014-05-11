import json, gzip, logging, os
from python.component_contribution import ComponentContribution

if not os.path.isdir('res'):
    os.mkdir('res')

OUTPUT_JSON = 'res/cc_compounds.json.gz'

if __name__ == '__main__':
    cc = ComponentContribution()
    cc.train()
    compound_json = []

    for i, compound_id in enumerate(cc.ccache.get_all_compound_ids()):
        logging.info("exporting " + compound_id)

        # skip compounds that cause a segmentation fault in openbabel
        if compound_id in ['C09078', 'C09093', 'C09145', 'C09246',
                           'C10282', 'C10286', 'C10356', 'C10359',
                           'C10396', 'C16818', 'C16839', 'C16857']: 
            continue
        compound_json.append(cc.get_compound_json(compound_id))
    new_json = gzip.open(OUTPUT_JSON, 'w')
    json.dump(compound_json, new_json, sort_keys=True, indent=4)
    new_json.close()