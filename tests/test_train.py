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

import unittest
from component_contribution.training_data import ToyTrainingData
from component_contribution import ComponentContribution, Reaction

class TestTraining(unittest.TestCase):
    
    td = ToyTrainingData()
    cc = ComponentContribution(cache_file_name=None, training_data=td)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    def test_train(self):
        self.assertEqual(self.cc.params['b'].shape, (48,))
        self.assertEqual(self.cc.params['dG0_rc'].shape, (97,))
        self.assertEqual(self.cc.params['dG0_gc'].shape, (172,))
        self.assertEqual(self.cc.params['dG0_cc'].shape, (97,))
        
    def test_dG0_calculation(self):
        formula = 'KEGG:C00002 + KEGG:C00001 <=> KEGG:C00008 + KEGG:C00009'
        rxn = Reaction.parse_formula(formula)
        dG0, sigma = self.cc.get_dG0_r(rxn)
        self.assertAlmostEqual(dG0, 12.6, 1)
        self.assertAlmostEqual(sigma, 1.8, 1)
        
    def test_dG0_prime_calculation(self):
        formula = 'KEGG:C00002 + KEGG:C00001 <=> KEGG:C00008 + KEGG:C00009'
        rxn = Reaction.parse_formula(formula)
        dG0_prime, sigma = self.cc.get_dG0_r_prime(rxn, pH=7.0, I=0.25, T=298.15)
        self.assertAlmostEqual(dG0_prime, -26.6, 1)
        self.assertAlmostEqual(sigma, 1.8, 1)
        
        dG0_prime_ph8, _ = self.cc.get_dG0_r_prime(rxn, pH=8.0, I=0.25, T=298.15)
        self.assertAlmostEqual(dG0_prime_ph8, -31.5, 1)

        dG0_prime_ph6, _ = self.cc.get_dG0_r_prime(rxn, pH=6.0, I=0.25, T=298.15)
        self.assertAlmostEqual(dG0_prime_ph6, -24.2, 1)
