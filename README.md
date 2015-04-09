component-contribution
======================

Standard reaction Gibbs energy estimation for biochemical reactions

Requirements (for the python version):
* python == 2.7
* numpy >= 1.6.2
* scipy >= 0.14
* oct2py >= 3.1 (dependent on Octave)
* Open-Babel >= 2.3.1 ('babel' must be in PATH, including python bindings)
* ChemAxon's Marvin >= 5.11 ('cxcalc' must be in PATH)

Requirements (for the MATLAB version):
* Matlab >= 7
* python == 2.7
* numpy >= 1.6.2
* scipy >= 0.14
* Open-Babel >= 2.3.1 ('babel' must be in PATH, including python bindings)
* ChemAxon's Marvin >= 5.11 ('cxcalc' must be in PATH)

For more information on the method behind component-contribution, please view our open access paper:

Noor E, Haraldsdóttir HS, Milo R, Fleming RMT (2013)
Consistent Estimation of Gibbs Energy Using Component Contributions,
PLoS Comput Biol 9:e1003098, DOI: 10.1371/journal.pcbi.1003098

Please, quote this paper if you publish work that uses component-contribution.

Description of files in /data:
------------------------------
* compounds.tsv - table of all KEGG compounds, their name and InChI (used only in Matlab code)
* fixed_mapping.tsv - table mapping some KEGG compound IDs to BiGG IDs (overriding the InChI-based mapping)
* formation_energies_transformed.tsv - table of biochemical formation energies (used for training CC)
* kegg_additions.tsv - table of compounds that are missing from KEGG together with their InChI
* kegg_compounds.json.gz - JSON of all KEGG compounds including their InChI and names
* redox.tsv - table of reduction potentials (used for training CC)
* TECRDB.tsv - table of K'eq values from the NIST database (http://xpdb.nist.gov/enzyme_thermodynamics/)
