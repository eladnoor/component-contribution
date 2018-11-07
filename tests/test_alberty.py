import unittest
from component_contribution import CompoundCache, ChemAxonNotFoundError

class TestAlbertyTransform(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestAlbertyTransform, self).__init__(*args, **kwargs)
        self.ccache = CompoundCache()
        try:
            self.atp_comp = self.ccache.get_compound('KEGG:C00002')
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
