# The MIT License (MIT)
#
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from component_contribution import (
    ChemAxonNotFoundError, Compound, GroupDecompositionError, ccache, inchi2gv)


class TestGroupDecomposition:

    def __init__(self, *args, **kwargs):
        super(TestGroupDecomposition, self).__init__(*args, **kwargs)

        groups_data = inchi2gv.init_groups_data()
        self.group_names = groups_data.GetGroupNames()
        self.decomposer = inchi2gv.InChIDecomposer(groups_data)

    def test_atp_major_ms(self):
        inchi_ATP = 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1'
        smiles_ATP_neutral = 'Nc1c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O)c2ncn1'
        self.assertDictEqual(self.decompose_inchi(inchi_ATP),
                             self.decompose_smiles(smiles_ATP_neutral))
        try:
            atp_cpd = ccache.get_compound('KEGG:C00002')
            self.assertIsInstance(atp_cpd, Compound)
            self.assertEqual(inchi_ATP, atp_cpd.inchi)
        except ChemAxonNotFoundError:
            self.skipTest('cannot test without ChemAxon')

    def decompose_cid(self, cid):
        cpd = ccache.get_compound('KEGG:' + cid)
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
        try:
            atp_groups_res = self.decompose_cid('C00002')
        except ChemAxonNotFoundError:
            smiles = 'Nc1c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O)c2ncn1'
            atp_groups_res = self.decompose_smiles(smiles)

        atp_groups_ref = {'-C- [H2 Z0 Mg0]': 1,
                          '-OPO2-OPO2- [H2 Z0 Mg0]': 1,
                          '-OPO3 [H2 Z0 Mg0]': 1,
                          '2-ring =c< [H0 Z0 Mg0]': 2,
                          'Origin [H0 Z0 Mg0]': 1,
                          'ring -C< [H1 Z0 Mg0]': 2,
                          'ring -O- [H0 Z0 Mg0]': 1,
                          'ring -n= [H0 Z0 Mg0]': 3,
                          'ring =c- [H1 Z0 Mg0]': 2,
                          'ring =c< [H0 Z0 Mg0]': 1,
                          'ring =n< [H0 Z0 Mg0]': 1,
                          'ring >C-O [H2 Z0 Mg0]': 2,
                          'ring >c-N [H2 Z0 Mg0]': 1}

        self.assertDictEqual(atp_groups_res, atp_groups_ref)

    def test_phosphate(self):
        # test the decomposition of Pi into groups
        try:
            pi_groups_res = self.decompose_cid('C00009')
        except ChemAxonNotFoundError:
            pi_groups_res = self.decompose_smiles('OP(=O)(O)O')
        self.assertIsNone(pi_groups_res)

    def test_acetate(self):
        try:
            ace_groups_res = self.decompose_cid('C00033')
        except ChemAxonNotFoundError:
            ace_groups_res = self.decompose_smiles('CC(=O)O')

        ace_groups_ref = {'-C [H3 Z0 Mg0]': 1,
                          '-COO [H1 Z0 Mg0]': 1,
                          'Origin [H0 Z0 Mg0]': 1}
        self.assertDictEqual(ace_groups_res, ace_groups_ref)
