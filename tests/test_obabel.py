import unittest

class TestOpenBabel(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(TestOpenBabel, self).__init__(*args, **kwargs)
        
    def test_openbabel(self):
        import pybel

        smiles = "C[C@H](O)C(=O)CC(=O)C(=O)[O-]"
        inchi = 'InChI=1S/C6H8O5/c1-3(7)4(8)2-5(9)6(10)11/h3,7H,2H2,1H3,(H,10,11)/p-1/t3-/m0/s1'
        
        mol = pybel.readstring('smiles', smiles)
        res = mol.write('inchi').strip()
        
        self.assertEqual(res, inchi)
