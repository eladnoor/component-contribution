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

from component_contribution.chemaxon import (
    ChemAxonNotFoundError, get_dissociation_constants, get_formula_and_charge)


class TestChemaxon:

    def __init__(self, *args, **kwargs):
        super(TestChemaxon, self).__init__(*args, **kwargs)

    def test_pka(self):
        try:
            compound_list = [('InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-3', 'O4P', -3, [1.80, 6.95, 12.90]), # orthophosphate
                             ('InChI=1S/C4H8O4/c5-1-3(7)4(8)2-6/h3,5-7H,1-2H2/t3-/m1/s1', 'C4H8O4', 0, [-8.53, -3.86, -3.33, -3.01, 12.54, 13.9, 15.45])] # D-Erythrulose
            for inchi, formula, charge, diss_table in compound_list:
                _formula, _charge = get_formula_and_charge(inchi)
                _diss_table, major_ms = get_dissociation_constants(inchi)
                self.assertEqual(formula, _formula)
                self.assertEqual(charge, _charge)
                self.assertEqual(len(diss_table), len(_diss_table))
                for pka1, pka2 in zip(sorted(diss_table), sorted(_diss_table)):
                    self.assertAlmostEqual(pka1, pka2, 2)

        except ChemAxonNotFoundError:
            self.skipTest('cannot test without ChemAxon')
