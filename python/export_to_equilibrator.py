import json, gzip, numpy, logging
from python.compound_cacher import CompoundCacher
from python.training_data import TrainingData
from python.component_contribution import ComponentContribution
from python.inchi2gv import init_groups_data, InChIDecomposer, GroupDecompositionError
from python.thermodynamic_constants import default_T


def update_json_data(d, dG0_rc, decomposer):
    """
        adds the component-contribution estimation to the JSON
    """
    compound_id = d['CID']
    comp = ccache.get_compound(compound_id)
    
    # decompose the compound and calculate the 'formation energy'
    # using the group contributions
    try:
        group_vec = decomposer.smiles_to_groupvec(comp.smiles_pH7)
        g = group_vec.ToArray()
        major_ms_dG0_f = float(numpy.dot(g, dG0_rc[:g.shape[1],0]))

        pmap = {'priority': 3,
                'source': 'New Component Contribution (2014)',
                'species': list(comp.get_species(major_ms_dG0_f, default_T))}
        d['pmaps'] = d.setdefault('pmaps', []) + [pmap]
    except GroupDecompositionError:
        pass
        
    return d

if __name__ == '__main__':
    td = TrainingData()
    cc = ComponentContribution(td)
    cc.train_without_model()
    dG0_rc = cc.params['dG0_rc']

    ccache = CompoundCacher()
    groups_data = init_groups_data()
    decomposer = InChIDecomposer(groups_data)
    group_names = groups_data.GetGroupNames()

    compound_json = []

    for i, d in enumerate(json.load(gzip.open('data/equilibrator_compounds.json.gz','r'))):
        compound_id = d['CID']
        logging.info("exporting " + compound_id)
        update_json_data(d, dG0_rc, decomposer)
        compound_json.append(d)
        if i == 10:
            break
    new_json = gzip.open('data/cc_compounds.json.gz', 'w')
    json.dump(compound_json, new_json, sort_keys=True, indent=4)
    new_json.close()