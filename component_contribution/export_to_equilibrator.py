import logging
import os

from pkg_resources import resource_filename

from component_contribution.component_contribution_trainer import (
    ComponentContribution)


logger = logging.getLogger(__name__)


if not os.path.isdir('res'):
    os.mkdir('res')

CACHE_EQUILIBRATOR_JSON_FNAME = \
    resource_filename('component_contribution',
                      'cache/cc_equilibrator.json.gz')

CACHE_EQUILIBRATOR_NPZ_FNAME = \
    resource_filename('component_contribution',
                      'cache/cc_equilibrator.npz')

if __name__ == '__main__':
    logger.setLevel(logging.INFO)
    cc = ComponentContribution()

    logging.info("Writing metadata to: " + CACHE_EQUILIBRATOR_JSON_FNAME)
    cc.save_equilibrator_metadata(CACHE_EQUILIBRATOR_JSON_FNAME)

    logging.info("Writing binary data to: " + CACHE_EQUILIBRATOR_NPZ_FNAME)
    cc.save_equilibrator_bin_data(CACHE_EQUILIBRATOR_NPZ_FNAME)

    logging.info("Now run 'clear_database' and then 'load_database' in "
                 "eQuilibrator")
