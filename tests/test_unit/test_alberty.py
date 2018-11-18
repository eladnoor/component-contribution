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

from component_contribution import ChemAxonNotFoundError, ccache


try:
    ccache.get_compound('KEGG:C00002')
except ChemAxonNotFoundError:
    pytest.skip("Cannot test without ChemAxon.", allow_module_level=True)


@pytest.fixture(scope="module")
def atp_comp():
    return ccache.get_compound('KEGG:C00002')


def test_atp(atp_comp):
    assert set(atp_comp.p_kas) == {12.6, 7.42, 4.93, 3.29, 1.55, 0.9}
    assert atp_comp.major_microspecies == 2


def test_transform(atp_comp):
    temperature = 300.0
    expected_delta_delta_g = [0, 72.3, 114.9, 143.2, 162.1, 171.0, 176.2]
    assert len(atp_comp.number_of_protons) == len(expected_delta_delta_g)
    delta_delta_g = [atp_comp._ddG(0, i, temperature) for i in range(7)]
    assert delta_delta_g == pytest.approx(expected_delta_delta_g)
