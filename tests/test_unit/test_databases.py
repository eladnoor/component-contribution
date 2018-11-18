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

import component_contribution.databases as db


def test_get_molecule(mocker):
    interface = db.DatabaseInterface()
    response = mocker.Mock(return_value=7)
    interface.register("mock", response)
    assert interface.get_molecule("mock", "foo") == 7
    response.assert_called_once_with("foo")
    with pytest.raises(KeyError):
        interface.get_molecule("candyland", "foo")


@pytest.mark.parametrize("compound_id, url", [
    ("C14818", "http://rest.kegg.jp/get/cpd:C14818/mol"),
    ("C00237", "http://rest.kegg.jp/get/cpd:C00237/mol"),
])
def test_get_kegg_molecule(mocker, compound_id, url):
    mocked_get = mocker.patch("component_contribution.databases.get")
    mocked_readstring = mocker.patch(
        "component_contribution.databases.readstring")
    db.get_kegg_molecule(compound_id)
    mocked_get.assert_called_once_with(url)
    assert mocked_readstring.call_count == 1


@pytest.mark.parametrize("compound_id, url", [
    ("HMDB0000122",
     "http://www.hmdb.ca/structures/metabolites/HMDB0000122.mol"),
    ("HMDB0000108",
     "http://www.hmdb.ca/structures/metabolites/HMDB0000108.mol"),
    ("HMDB0108",
     "http://www.hmdb.ca/structures/metabolites/HMDB0000108.mol"),
])
def test_get_hmdb_molecule(mocker, compound_id, url):
    mocked_get = mocker.patch("component_contribution.databases.get")
    mocked_readstring = mocker.patch(
        "component_contribution.databases.readstring")
    db.get_hmdb_molecule(compound_id)
    mocked_get.assert_called_once_with(url)
    assert mocked_readstring.call_count == 1


@pytest.mark.parametrize("compound_id", [
    "foo"
    "bar"
])
def test_get_inchi_molecule(mocker, compound_id):
    mocked_readstring = mocker.patch(
        "component_contribution.databases.readstring")
    db.get_inchi_molecule(compound_id)
    mocked_readstring.assert_called_once_with("inchi", compound_id)


@pytest.mark.parametrize("compound_id, inchi", [
    ("CHEBI:15846", "InChI=1S/C21H27N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,(H5-,22,23,24,25,33,34,35,36,37)/p+1/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1"),
    ("CHEBI:15422", "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"),
])
def test_get_chebi_molecule(mocker, compound_id, inchi):
    mocked_readstring = mocker.patch(
        "component_contribution.databases.readstring")
    db.get_chebi_molecule(compound_id)
    mocked_readstring.assert_called_once_with("inchi", inchi)
