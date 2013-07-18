import sys, logging
sys.path.append('../python')
from compound import Compound
from inchi2gv import init_groups_data, InChI2GroupVector, GroupDecompositionError
from compound_cacher import CompoundCacher
from molecule import Molecule

#logger = logging.getLogger('')
#logger.setLevel(logging.DEBUG)
ccache = CompoundCacher.getInstance('../cache/compounds.json')
groups_data = init_groups_data()
group_list = groups_data.GetGroupNames()
inchi2gv_converter = InChI2GroupVector(groups_data)

patterns = ['c~[O;+0]', 'c~[O;+1]', 'c~[n;+1]~c', 'c~[n;+0]~c', 'c~[n;-1]~c']

for cid in [255, 1007]:
    comp = ccache.get_kegg_compound(cid)
    print "-"*50, '\nC%05d' % cid
    inchi = comp.inchi
    mol = Molecule.FromInChI(inchi)
    print mol.ToSmiles()
    
    print mol.FindSmarts("c~[n;+1]~c")
    
    try:
        groupvec = inchi2gv_converter.InChI2GroupVector(inchi)
        sys.stdout.write(str(groupvec) + '\n')
    except GroupDecompositionError as e:
        sys.stderr.write(str(e) + '\n')
        sys.stderr.write(e.GetDebugTable())
    