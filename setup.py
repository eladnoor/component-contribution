#!/usr/bin/env python


import versioneer
from setuptools import setup

mydata_files = ['TECRDB.csv',
                'formation_energies_transformed.csv',
                'redox.csv']
data_files = [('/train/tecrdb', ['component_contribution/data/TECRDB.csv']),
              ('/train/formation', ['component_contribution/data/formation_energies_transformed.csv']),
              ('/train/redox', ['component_contribution/data/redox.csv']),
              ('/train/toy', ['component_contribution/data/toy_training_data.csv'])]

# Most arguments are set in the `setup.cfg`.
setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    data_files=data_files
)
