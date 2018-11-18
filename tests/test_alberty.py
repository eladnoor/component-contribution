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


from component_contribution import ChemAxonNotFoundError, ccache


class TestAlbertyTransform:

    def __init__(self, *args, **kwargs):
        super(TestAlbertyTransform, self).__init__(*args, **kwargs)
        try:
            self.atp_comp = ccache.get_compound('KEGG:C00002')
            self.missing_chemaxon = False
        except ChemAxonNotFoundError:
            self.missing_chemaxon = True

    def test_atp(self):
        if self.missing_chemaxon:
            self.skipTest('cannot test without ChemAxon')
        self.assertSetEqual(set([12.6, 7.42, 4.93, 3.29, 1.55, 0.9]),
                            set(self.atp_comp.p_kas))
        self.assertEqual(self.atp_comp.major_microspecies, 2)

    def test_transform(self):
        if self.missing_chemaxon:
            self.skipTest('cannot test without ChemAxon')
        T = 300.0
        ddG_ref_list = [0, 72.3, 114.9, 143.2, 162.1, 171.0, 176.2]

        self.assertEqual(len(self.atp_comp.number_of_protons),
                         len(ddG_ref_list))

        ddG_res_list = [self.atp_comp._ddG(0, i, T) for i in range(7)]
        for i in range(7):
            self.assertAlmostEqual(ddG_ref_list[i], ddG_res_list[i], 1)
