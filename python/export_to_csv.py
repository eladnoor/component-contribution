import logging, os, csv
from python.component_contribution import ComponentContribution
from python.thermodynamic_constants import default_T

if not os.path.isdir('res'):
    os.mkdir('res')

OUTPUT_CSV = 'res/cc_compounds.csv'

if __name__ == '__main__':
    cc = ComponentContribution.init()
    output_csv = csv.writer(open(OUTPUT_CSV, 'w'))
    output_csv.writerow(['Compound ID', 'nH', 'charge', 'dG0_f'])

    for i, compound_id in enumerate(cc.ccache.get_all_compound_ids()):
        logging.info("exporting " + compound_id)

        # skip compounds that cause a segmentation fault in openbabel
        if compound_id in ['C09078', 'C09093', 'C09145', 'C09246',
                           'C10282', 'C10286', 'C10356', 'C10359',
                           'C10396', 'C16818', 'C16839', 'C16857']: 
            continue

        comp = cc.ccache.get_compound(compound_id)
        major_ms_dG0_f = cc.get_major_ms_dG0_f(compound_id)
        for d in comp.get_species(major_ms_dG0_f, default_T):
            output_csv.writerow([compound_id, d['nH'], d['z'], d['dG0_f']])