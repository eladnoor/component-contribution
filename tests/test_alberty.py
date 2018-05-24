import unittest
from component_contribution import CompoundCache

class TestAlbertyTransform(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestAlbertyTransform, self).__init__(*args, **kwargs)
        self.ccache = CompoundCache()
        self.atp_comp = self.ccache.get_compound('KEGG:C00002')

    def test_atp(self):
        self.assertSetEqual(set([12.6, 7.42, 4.93, 3.29, 1.55, 0.9]),
                            set(self.atp_comp.p_kas))
        self.assertEqual(self.atp_comp.major_microspecies, 2)
        
    def test_transform(self):
        pH = 7.0
        I = 0.2
        T = 300.0
        
        species_transform_values = [(0, 406.3),
                                    (1, 478.6),
                                    (2, 521.2),
                                    (3, 549.5),
                                    (4, 568.4),
                                    (5, 577.3),
                                    (6, 582.5)]
        
        for i, t in species_transform_values:
            self.assertAlmostEqual(self.atp_comp.transform(i, pH, I, T), t, 1)

