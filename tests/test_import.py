import unittest
from distutils.version import StrictVersion

class TestDependenceis(unittest.TestCase):

    def test_pybel(self):
        # Test that openbabel can be imported
        try:
            import pybel
        except ImportError:
            self.fail('''OpenBabel is not installed, or was installed without 
                         python bindings. Please go to http://openbabel.org/wiki/Python
                         and follow the installation instructions.''')
            
                       
    def test_numpy(self):
        # Test that numpy can be imported and its version is rather new
        try:
            import numpy
            if StrictVersion(numpy.__version__) < StrictVersion('1.6.2'):
                self.fail('''WARNING: your NumPy version is lower than 1.6.2
                             and might not work properly. Please upgrade to
                             a newer version.''')
        except ImportError:
            self.fail('''NumPy is not installed. Please go to http://www.numpy.org
                         and follow the installation instructions.''')
