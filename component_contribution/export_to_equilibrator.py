import logging, os
from component_contribution.component_contribution_trainer import ComponentContribution

if not os.path.isdir('res'):
    os.mkdir('res')

OUTPUT_JSON = 'res/cc_compounds.json.gz'
OUTPUT_NPZ = 'res/cc_preprocess'

if __name__ == '__main__':
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    cc = ComponentContribution.init()
    
    logging.info("Wrote compound data to: " + OUTPUT_JSON)
    cc.save_compound_data(OUTPUT_JSON)  
    
    logging.info("Wrote preprocessing data to: " + OUTPUT_NPZ)
    cc.save_preprocessing_data(OUTPUT_NPZ)
    
    logging.info("Now run 'clear_database' and then 'load_database' in eQuilibrator")
