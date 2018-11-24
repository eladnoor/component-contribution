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
    expected_transforms = [(5.0, 0.0, 300.0, 371.6),
                           (7.0, 0.0, 300.0, 521.6),
                           (9.0, 0.0, 300.0, 662.5),
                           (7.0, 0.5, 300.0, 520.6)]
    for p_h, ionic_strength, temperature, expected_t in expected_transforms:
        assert atp_comp.transform(p_h, ionic_strength, temperature) == \
            pytest.approx(expected_t, rel=1e-3)


def test_transform(atp_comp):
    temperature = 300.0
    expected_delta_delta_g = [0, 72.3, 114.9, 143.2, 162.1, 171.0, 176.2]
    expected_delta_delta_g = [114.9 - x for x in expected_delta_delta_g]
    assert len(atp_comp.species) == len(expected_delta_delta_g)
    delta_delta_g = [ms.ddG_over_T * temperature for ms in atp_comp.species]
    assert delta_delta_g == pytest.approx(expected_delta_delta_g, rel=1e-3)
