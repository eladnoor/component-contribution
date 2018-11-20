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

from component_contribution.chemaxon import (
    verify_cxcalc, get_dissociation_constants, get_formula_and_charge)


if not verify_cxcalc():
    pytest.skip("Cannot test without ChemAxon.", allow_module_level=True)


@pytest.mark.parametrize("inchi, exp_formula, exp_charge, exp_diss_table", [
    # Orthophosphate
    ('InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-3',
     'O4P',
     -3,
     [1.80, 6.95, 12.90]),
    # D-Erythrulose
    ('InChI=1S/C4H8O4/c5-1-3(7)4(8)2-6/h3,5-7H,1-2H2/t3-/m1/s1',
     'C4H8O4',
     0,
     [-8.53, -3.86, -3.33, -3.01, 12.54, 13.9, 15.45]),
])
def test_pka(inchi, exp_formula, exp_charge, exp_diss_table):
    formula, charge = get_formula_and_charge(inchi)
    assert formula == exp_formula
    assert charge == exp_charge
    diss_table, major_ms = get_dissociation_constants(inchi)
    assert len(diss_table) == len(exp_diss_table)
    assert sorted(diss_table) == pytest.approx(sorted(exp_diss_table))
