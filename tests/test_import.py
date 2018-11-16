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
