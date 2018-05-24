# -*- encoding: utf-8 -*-

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


def test_get_chebi_molecule():
    with pytest.raises(NotImplementedError):
        db.get_chebi_molecule(None)
