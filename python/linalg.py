# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 11:21:14 2015

@author: noore
"""

try:
   import oct2py
   SVD_METHOD = 'octave'
except ImportError:
   raise Exception('It seems that you do not have Octave and/or Oct2Py installed. '
                   'Calculating the Component-Contributions requires these packages.')

import numpy as np

class LINALG(object):
    oc = oct2py.Oct2Py()

    @staticmethod    
    def _svd_octave(A):
        U, S, V = LINALG.oc.svd(A)
        S = np.matrix(S)
        U = np.matrix(U)
        V = np.matrix(V)
        return U, S, V
        
    @staticmethod
    def _svd_dgesdd(A):

        # numpy.linalg.svd returns U, s, V such that
        # A = U * s * V

        # however, matlab and octave return U, S, V such that
        # V needs to be conjugate transposed when multiplied:
        # A = U * S * V.H

        # we would like to stick to the latter standard, so we return
        # the transposed V here (assuming it is real)

        U, s, V = np.linalg.svd(A, full_matrices=True)
        S = np.matrix(np.diag(s))
        U = np.matrix(U)
        V = np.matrix(V)
        return U, S, V.T

    @staticmethod
    def _svd_dgesvd(A):
        """
            This LAPACK function is the one that matches the matlab
            and octave implementations, but so far I couldn't find 
            a good python wrapper for it, except for oct2py (which
            is difficult to install on some systems and requires
            the entire Octave framework).
        """
        pass


    @staticmethod
    def svd(A):
        return LINALG._svd_octave(A)