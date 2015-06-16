component-contribution
======================

Standard reaction Gibbs energy estimation for biochemical reactions

### Requirements (for the python version):
* python == 2.7
* numpy >= 1.6.2
* scipy >= 0.14
* oct2py >= 3.1 (dependent on Octave)
* Open-Babel >= 2.3.1
* ChemAxon's Marvin >= 5.11

### Requirements (for the MATLAB version):
* Matlab >= 7
* python == 2.7
* numpy >= 1.6.2
* scipy >= 0.14
* Open-Babel >= 2.3.1 ('babel' must be in PATH, including python bindings)
* ChemAxon's Marvin >= 5.11 ('cxcalc' must be in PATH)

For more information on the method behind component-contribution, please view our open access paper:

Noor E, Haraldsdóttir HS, Milo R, Fleming RMT (2013)
[Consistent Estimation of Gibbs Energy Using Component Contributions](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003098),
PLoS Comput Biol 9:e1003098, DOI: 10.1371/journal.pcbi.1003098

Please, quote this paper if you publish work that uses component-contribution.

## Description of files in /data
* compounds.tsv - table of all KEGG compounds, their name and InChI (used only in Matlab code)
* fixed_mapping.tsv - table mapping some KEGG compound IDs to BiGG IDs (overriding the InChI-based mapping)
* formation_energies_transformed.tsv - table of biochemical formation energies (used for training CC)
* kegg_additions.tsv - table of compounds that are missing from KEGG together with their InChI
* kegg_compounds.json.gz - JSON of all KEGG compounds including their InChI and names
* redox.tsv - table of reduction potentials (used for training CC)
* TECRDB.tsv - table of K'eq values from the NIST database (http://xpdb.nist.gov/enzyme_thermodynamics/)

## Installing on Windows
1. install 32-bit python2.7
  * recommended version [Winpython](http://winpython.github.io)
  * you must install the 32-bit version since OpenBabel is not compiled for 64-bit windows
  * make sure to register the python interpreter before installing openbabel python bindings
  * add `python.exe` to PATH: [instructions](http://docs.python.org/2/using/windows.html#excursus-setting-environment-variables)
2. install OpenBabel (version 2.3.2)
  * [instructions](http://openbabel.org/wiki/Category:Installation)
  * [installer](http://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/OpenBabel2.3.2a_Windows_Installer.exe/download)
3. install OpenBabel python (version 1.8) bindings
  * [instructions](http://open-babel.readthedocs.org/en/latest/UseTheLibrary/PythonInstall.html#windows)
  * [installer](http://sourceforge.net/projects/openbabel/files/openbabel-python/1.8/openbabel-python-1.8.py27.exe/download)
4. optional: install Octave
  * Not required for the Matlab version of component-contribution
  * [instructions](http://blink1073.github.io/oct2py/source/installation.html)
  * there is no need to use "pip" to install oct2py as it is already part of winpython
  * [installer](http://sourceforge.net/projects/octave/files/Octave%20Windows%20binaries/)
  * add `octave.exe` to PATH
5. optional: install "Marvin Suite" by ChemAxon
  * Marvin is only required for adding structures of new compounds that are not in the KEGG database
  * [instructions](http://www.chemaxon.com/download/marvin-suite/)
  * add `cxcalc.bat` to PATH
  * you will need to get a license to use ChemAxon (it is free for academic use)

## Installing on Ubuntu
1. as root: `apt-get install openbabel octave`
2. as root: `pip install -U numpy scipy oct2py`
3. optional: install "Marvin Suite" by ChemAxon
  * Marvin is only required for adding structures of new compounds that are not in the KEGG database
  * [instructions](http://www.chemaxon.com/download/marvin-suite/)
  * add `cxcalc` to PATH, e.g. using a symbolic link from `/usr/bin/cxcalc`
  * you will need to get a license to use ChemAxon (it is free for academic use)
