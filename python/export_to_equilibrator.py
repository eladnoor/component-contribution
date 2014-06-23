import json, gzip, logging, os
import numpy as np
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
        sub_json = cc.get_compound_json(compound_id)
        if compound_id in ['C00001']:
            sub_json['pmap']['species'][0]['phase'] = 'liquid'
        if compound_id in ['C00087']:
            sub_json['pmap']['species'][0]['phase'] = 'solid'
        if compound_id in ['C00007', 'C00697']:
            sub_json['pmap']['species'].append({'phase':'gas', 'dG0_f':0})
        if compound_id in ['C00282']:
            sub_json['pmap']['species'].append({'phase':'gas', 'dG0_f':0, 'nH':2})
        if compound_id in ['C00011']:
            sub_json['pmap']['species'].append({'phase':'gas', 'dG0_f':-394.36})
        if compound_id in ['C00237']:
            sub_json['pmap']['species'].append({'phase':'gas', 'dG0_f':-137.17})
            
        compound_json.append(sub_json)
    new_json = gzip.open(OUTPUT_JSON, 'w')
    json.dump(compound_json, new_json, sort_keys=True, indent=4)
    new_json.close()
    
    v_r = cc.params['preprocess_v_r']
    v_g = cc.params['preprocess_v_g']
    C1  = cc.params['preprocess_C1']
    C2  = cc.params['preprocess_C2']
    C3  = cc.params['preprocess_C3']

    np.savez_compressed('res/cc_preprocess', v_r=v_r, v_g=v_g, C1=C1, C2=C2, C3=C3)
