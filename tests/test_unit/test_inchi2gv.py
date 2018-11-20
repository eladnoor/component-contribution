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


import pytest

from component_contribution import (
    ChemAxonNotFoundError, Compound, GroupDecompositionError, ccache, inchi2gv)


try:
    # TODO: Should introduce a utility that runs cxcalc --help to test.
    ccache.get_compound('KEGG:C00002')
except ChemAxonNotFoundError:
    pytest.skip("Cannot test without ChemAxon.", allow_module_level=True)


@pytest.fixture(scope="module")
def groups_data():
    return inchi2gv.init_groups_data()


@pytest.fixture(scope="module")
def group_names(groups_data):
    return groups_data.GetGroupNames()


@pytest.fixture(scope="module")
def decomposer(groups_data):
    return inchi2gv.InChIDecomposer(groups_data)


# TODO: Should this be a utility function in the package?
def decompose_inchi(inchi, group_names, decomposer):
    try:
        gr = decomposer.inchi_to_groupvec(inchi)
        # TODO: Why apply a filter? Should this be a dict-comprehension?
        return dict(filter(lambda x: x[1], zip(group_names, gr)))
    except GroupDecompositionError:
        return None


def decompose_cid(cid, group_names, decomposer):
    cpd = ccache.get_compound('KEGG:' + cid)
    assert isinstance(cpd, Compound)
    groups = decompose_inchi(cpd.inchi, group_names, decomposer)
    del cpd
    return groups


def decompose_smiles(smiles, group_names, decomposer):
    try:
        gr = decomposer.smiles_to_groupvec(smiles)
        return dict(filter(lambda x: x[1], zip(group_names, gr)))
    except GroupDecompositionError:
        return None


def test_atp_major_ms(group_names, decomposer):
    atp_inchi = 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1'
    atp_smiles_neutral = 'Nc1c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O)c2ncn1'
    assert decompose_inchi(atp_inchi, group_names, decomposer) == \
        decompose_smiles(atp_smiles_neutral, group_names, decomposer)

    atp_cpd = ccache.get_compound('KEGG:C00002')
    assert isinstance(atp_cpd, Compound)
    assert atp_cpd.inchi == atp_inchi


@pytest.mark.parametrize("compound_id, compound_smiles, exp_groups", [
    # ATP
    ('C00002',
     'Nc1c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O)c2ncn1',
     {
         '-C- [H2 Z0 Mg0]': 1,
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
         'ring >c-N [H2 Z0 Mg0]': 1
     }),
    # Phosphate
    ('C00009',
     'OP(=O)(O)O',
     None),
    # Acetate
    ('C00033',
     'CC(=O)O',
     {
         '-C [H3 Z0 Mg0]': 1,
         '-COO [H1 Z0 Mg0]': 1,
         'Origin [H0 Z0 Mg0]': 1
     }),

])
def test_decomposition(compound_id, compound_smiles, exp_groups, group_names,
                       decomposer):
    """Test the decomposition of a compound into groups."""
    groups_res = decompose_cid(compound_id, group_names, decomposer)
    assert groups_res == exp_groups
    groups_res = decompose_smiles(compound_smiles, group_names, decomposer)
    assert groups_res == exp_groups
