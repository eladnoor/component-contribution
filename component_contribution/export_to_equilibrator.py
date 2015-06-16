import json, gzip, logging, os
import numpy as np
from .component_contribution_trainer import ComponentContribution

if not os.path.isdir('res'):
    os.mkdir('res')

OUTPUT_JSON = '../equilibrator/data/cc_compounds.json.gz'
OUTPUT_NPZ = '../equilibrator/data/cc_preprocess'

if __name__ == '__main__':
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    cc = ComponentContribution.init()

    # write the JSON file containing the 'additiona' data on all the compounds
    # in eQuilibrator (i.e. formula, mass, pKa values, etc.) 
    compound_json = []
    for i, compound_id in enumerate(cc.ccache.get_all_compound_ids()):
        logging.debug("exporting " + compound_id)

        # skip compounds that cause a segmentation fault in openbabel
        if compound_id in ['C09078', 'C09093', 'C09145', 'C09246',
                           'C10282', 'C10286', 'C10356', 'C10359',
                           'C10396', 'C16818', 'C16839', 'C16857']: 
            continue
        sub_json = cc.get_compound_json(compound_id)
        if compound_id in ['C00001']: # override H2O as liquid phase only
            sub_json['pmap']['species'][0]['phase'] = 'liquid'
        if compound_id in ['C00087']: # override Sulfur as solid phase only
            sub_json['pmap']['species'][0]['phase'] = 'solid'
        if compound_id in ['C00007', 'C00697']: # add gas phase for O2 and N2
            sub_json['pmap']['species'].append({'phase':'gas', 'dG0_f':0})
        if compound_id in ['C00282']: # add gas phase for H2
            sub_json['pmap']['species'].append({'phase':'gas', 'dG0_f':0, 'nH':2})
        if compound_id in ['C00011']: # add gas phase for CO2
            sub_json['pmap']['species'].append({'phase':'gas', 'dG0_f':-394.36})
        if compound_id in ['C00237']: # add gas phase for CO
            sub_json['pmap']['species'].append({'phase':'gas', 'dG0_f':-137.17})
            
        compound_json.append(sub_json)
    new_json = gzip.open(OUTPUT_JSON, 'w')
    json.dump(compound_json, new_json, sort_keys=True, indent=4)
    new_json.close()
    
    # write an NPZ file (Numpy binary files) of preprocessed data needed for 
    # running Component Contribution estimations quickly
    
    v_r = cc.params['preprocess_v_r']
    v_g = cc.params['preprocess_v_g']
    C1  = cc.params['preprocess_C1']
    C2  = cc.params['preprocess_C2']
    C3  = cc.params['preprocess_C3']
    G1  = cc.params['preprocess_G1']
    G2  = cc.params['preprocess_G2']
    G3  = cc.params['preprocess_G3']
    S   = cc.params['preprocess_S']
    S_count = cc.params['preprocess_S_count']
    cids = cc.params['cids']

    np.savez_compressed(OUTPUT_NPZ,
                        v_r=v_r, v_g=v_g, C1=C1, C2=C2, C3=C3,
                        G1=G1, G2=G2, G3=G3,
                        S=S, S_count=S_count, cids=cids)
                        
                        
    logging.info("Now run 'clear_database' and then 'load_database' in eQuilibrator")
