import unittest

class TestChemaxon(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestChemaxon, self).__init__(*args, **kwargs)
        
    def test_pka(self):
        from component_contribution.chemaxon import get_formula_and_charge, \
            get_dissociation_constants

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

if __name__ == "__main__":
    unittest.main()