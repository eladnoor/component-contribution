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

from component_contribution import ComponentContribution, Reaction
from component_contribution.training_data import ToyTrainingData


@pytest.fixture(scope="module")
def training_data():
    return ToyTrainingData()


@pytest.fixture(scope="module")
def comp_contribution(training_data):
    return ComponentContribution(cache_file_name=None,
                                 training_data=training_data)


@pytest.fixture(scope="module")
def reaction():
    formula = 'KEGG:C00002 + KEGG:C00001 <=> KEGG:C00008 + KEGG:C00009'
    return Reaction.parse_formula(formula)


def test_train(comp_contribution):
    assert comp_contribution.params['b'].shape == (48,)
    assert comp_contribution.params['dG0_rc'].shape == (97,)
    assert comp_contribution.params['dG0_gc'].shape == (172,)
    assert comp_contribution.params['dG0_cc'].shape == (97,)


def test_delta_g_zero_calculation(comp_contribution):
    formula = 'KEGG:C00002 + KEGG:C00001 <=> KEGG:C00008 + KEGG:C00009'
    rxn = Reaction.parse_formula(formula)
    delta_g_zero, sigma = comp_contribution.get_dG0_r(rxn)
    assert delta_g_zero == pytest.approx(12.6, rel=1e-2)
    assert sigma == pytest.approx(1.76, rel=1e-2)


@pytest.mark.parametrize(
    "ph_value, ionic_strength, temperature, exp_delta_g_zero_prime, exp_sigma",
    [
        (6.0, 0.25, 298.15, -24.2, 1.76),
        (7.0, 0.25, 298.15, -26.6, 1.76),
        (8.0, 0.25, 298.15, -31.5, 1.76),
    ])
def test_delta_g_zero_prime_calculation(
        ph_value, ionic_strength, temperature, exp_delta_g_zero_prime,
        exp_sigma, reaction, comp_contribution):
    delta_g_zero_prime, sigma = comp_contribution.get_dG0_r_prime(
        reaction, pH=ph_value, ionic_strength=ionic_strength, T=temperature)
    assert delta_g_zero_prime == pytest.approx(exp_delta_g_zero_prime,
                                               rel=1e-2)
    assert sigma == pytest.approx(exp_sigma, rel=1e-2)
