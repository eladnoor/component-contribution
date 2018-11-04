import unittest
from component_contribution import CompoundCache, inchi2gv, \
    Compound, GroupDecompositionError

import logging
logging.getLogger().setLevel(logging.INFO)

class TestGroupDecomposition(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestGroupDecomposition, self).__init__(*args, **kwargs)
        self.ccache = CompoundCache()
        groups_data = inchi2gv.init_groups_data()
        self.group_names = groups_data.GetGroupNames()
        self.decomposer = inchi2gv.InChIDecomposer(groups_data)
        
    def test_atp_major_ms(self):
        inchi_ATP = 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1'
        smiles_ATP_neutral = 'Nc1c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O)c2ncn1'

        atp_cpd = self.ccache.get_compound('KEGG:C00002')
        self.assertIsInstance(atp_cpd, Compound)
        self.assertEqual(inchi_ATP, atp_cpd.inchi)
        self.assertDictEqual(self.decompose_inchi(inchi_ATP),
                             self.decompose_smiles(smiles_ATP_neutral))
        
    def decompose_cid(self, cid):
        cpd = self.ccache.get_compound('KEGG:' + cid)
        self.assertIsInstance(cpd, Compound)
        groups = self.decompose_inchi(cpd.inchi)
        del cpd
        return groups
    
    def decompose_inchi(self, inchi):
        try:
            gr = self.decomposer.inchi_to_groupvec(inchi)
            return dict(filter(lambda x: x[1], zip(self.group_names, gr)))
        except GroupDecompositionError:   
            return None

    def decompose_smiles(self, smiles):
        try:
            gr = self.decomposer.smiles_to_groupvec(smiles)
            return dict(filter(lambda x: x[1], zip(self.group_names, gr)))
        except GroupDecompositionError:   
            return None

    def test_atp(self):
        # test the decomposition of ATP into groups
        atp_groups = self.decompose_cid('C00002')
        
        # test the decomposition of ATP into groups
        self.assertEqual(atp_groups['-C- [H2 Z0 Mg0]'], 1)
        self.assertEqual(atp_groups['-OPO2-OPO2- [H2 Z0 Mg0]'], 1)
        self.assertEqual(atp_groups['-OPO3 [H2 Z0 Mg0]'], 1)
        self.assertEqual(atp_groups['2-ring =c< [H0 Z0 Mg0]'], 2)
        self.assertEqual(atp_groups['ring -C< [H1 Z0 Mg0]'], 2)
        self.assertEqual(atp_groups['ring -n= [H0 Z0 Mg0]'], 3)
        self.assertEqual(atp_groups['ring >c-N [H2 Z0 Mg0]'], 1)
        self.assertEqual(atp_groups['ring >C-O [H2 Z0 Mg0]'], 2)
        self.assertEqual(atp_groups['ring =n< [H0 Z0 Mg0]'], 1)
        
    def test_phosphate(self):
        # test the decomposition of Pi into groups
        self.assertIsNone(self.decompose_cid('C00009'))

    def test_acetate(self):
        self.assertDictEqual(self.decompose_cid('C00033'),
            {'-C [H3 Z0 Mg0]': 1, '-COO [H1 Z0 Mg0]': 1, 'Origin [H0 Z0 Mg0]': 1})


