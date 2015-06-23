from setuptools import setup

setup(
    name = 'ComponentContribution',
    version = '1.0.1',
    author = 'Elad Noor',
    author_email='noor@imsb.biol.ethz.ch',
    description = 'Standard reaction Gibbs energy estimation for biochemical reactions',
    license = 'MIT',
    packages=['component_contribution'],
    url='https://github.com/eladnoor/component-contribution',
    install_requires=['scipy>=0.14.0',
                      'numpy>=1.6.2',
                      'oct2py==3.1.0'],
    data_files=[('data', ['data/TECRDB.tsv', 
                          'data/redox.tsv',
                          'data/formation_energies_transformed.tsv',
                          'data/equilibrator_compounds.json.gz',
                          'data/kegg_additions.tsv',
                          'data/kegg_compounds.json.gz']),
                ('cache', ['cache/compounds.json.gz',
                           'cache/component_contribution_python.mat'])
               ],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',
    
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
    
        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',
    
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
    ],
)

