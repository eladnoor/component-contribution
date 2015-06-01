import sys
from distutils.version import StrictVersion

err_num = 0
inchi_ATP = 'InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1/f/h18-19,21,23H,11H2'
smiles_ATP_pH7 = 'NC1=C2N=CN([C@@H]3O[C@H](COP([O-])(=O)OP([O-])(=O)OP(O)([O-])=O)[C@@H](O)[C@H]3O)C2=NC=N1'
ATP_group_dict = {43:1, 50:1, 72:1, 74:3, 103:1, 119:2, 127:1, 147:2, 148:2, 149:1, 157:2, 160:1, 162:1}

# Test that openbabel can be imported
try:
    import openbabel
except ImportError:
    sys.stderr.write('OpenBabel is not installed, or was installed without '
                     'python bindings. Please go to http://openbabel.org/wiki/Python '
                     'and follow the installation instructions.\n')
    err_num += 1
               
# Test that numpy can be imported and its version is rather new
try:
    import numpy
    if StrictVersion(numpy.__version__) < StrictVersion('1.6.2'):
        sys.stderr.write('WARNING: your NumPy version is lower than 1.6.2 '
                         'and might not work properly. Please upgrade to '
                         'a newer version.\n')
except ImportError:
    sys.stderr.write('NumPy is not installed. Please go to http://www.numpy.org '
                     'and follow the installation instructions.\n')
    err_num += 1

try:
    from component_contribution.compound_cacher import CompoundCacher
    ccache = CompoundCacher()
    atp_comp = ccache.get_compound('C00002')
    assert(smiles_ATP_pH7 == atp_comp.smiles_pH7)
except AssertionError:
    sys.stderr.write('Internal Error: the SMILES string for ATP is wrong.\n')
    err_num += 1
except Exception as e:
    sys.stderr.write('Error using Compound Cacher: ' + str(e))
    err_num += 1

# Test inchi2gv.py
try:
    from component_contribution import inchi2gv    
    groups_data = inchi2gv.init_groups_data()
    decomposer = inchi2gv.InChIDecomposer(groups_data)
    groupvec1 = decomposer.inchi_to_groupvec(inchi_ATP)
    groupvec2 = decomposer.smiles_to_groupvec(smiles_ATP_pH7)
    print groupvec1
    print groupvec2
    for group_ind, group_count in enumerate(groupvec1.Flatten()):
        assert(ATP_group_dict.get(group_ind, 0) == group_count)
        
except ImportError:
    sys.stderr.write('Cannot import the python script inchi2gv. Make sure the file '
                     'inchi2gv.py is located in the folder '
                     '"component-contribution/python/".\n')
    err_num += 1
except inchi2gv.GroupDecompositionError as e:
    sys.stderr.write('Internal Error: cannot decompose the compound ATP.\n')
    err_num += 1
except AssertionError:
    sys.stderr.write('Internal Error: the decomposition of ATP is not correct.\n')
    err_num += 1

if err_num == 0:
    sys.stderr.write('Success\n')
sys.exit(err_num)
