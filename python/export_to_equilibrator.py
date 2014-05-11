import json, gzip, logging
from python.component_contribution import ComponentContribution
from python.thermodynamic_constants import default_T

INPUT_JSON = 'data/equilibrator_compounds.json.gz'
OUTPUT_JSON = 'data/cc_compounds.json.gz'

def update_json_data(d, cc):
    """
        adds the component-contribution estimation to the JSON
    """
    compound_id = d['CID']
    
    # decompose the compound and calculate the 'formation energy'
    # using the group contributions
    major_ms_dG0_f = cc.get_major_ms_dG0_f(compound_id)
    comp = cc.ccache.get_compound(compound_id)

    pmap = {'priority': 3,
            'source': 'New Component Contribution (2014)',
            'species': list(comp.get_species(major_ms_dG0_f, default_T))}
    d['pmaps'] = d.setdefault('pmaps', []) + [pmap]
        
    return d

if __name__ == '__main__':
    cc = ComponentContribution()
    cc.train()
    compound_json = []

    for i, d in enumerate(json.load(gzip.open(INPUT_JSON,'r'))):
        compound_id = d['CID']
        logging.info("exporting " + compound_id)
        update_json_data(d, cc)
        compound_json.append(d)
    new_json = gzip.open(OUTPUT_JSON, 'w')
    json.dump(compound_json, new_json, sort_keys=True, indent=4)
    new_json.close()