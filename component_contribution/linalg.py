# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 11:21:14 2015

@author: noore
"""

import numpy as np
import scipy

class LINALG(object):

    @staticmethod
    def svd(A):

        # numpy.linalg.svd returns U, s, V such that
        # A = U * s * V

        # however, matlab and octave return U, S, V such that
        # V needs to be conjugate transposed when multiplied:
        # A = U * S * V.H

        # we would like to stick to the latter standard, so we return
        # the transposed V here (assuming it is real)

        U, s, V = scipy.linalg.svd(A, full_matrices=True)
        S = np.matrix(np.zeros(A.shape))
        np.fill_diagonal(S, s)
        U = np.matrix(U)
        V = np.matrix(V)
        return U, S, V.T

    @staticmethod
    def _zero_pad_S(S, cids_orig, cids_joined):
        """
            takes a stoichiometric matrix with a given list of IDs 'cids' and adds
            0-rows so that the list of IDs will be 'cids_joined'
        """
        if not set(cids_orig).issubset(cids_joined):
            raise Exception('The full list is missing some IDs in "cids"')

        full_S = np.zeros((len(cids_joined), S.shape[1]))
        for i, cid in enumerate(cids_orig):
            S_row = S[i, :]
            full_S[cids_joined.index(cid), :] = S_row

        return np.matrix(full_S)

    @staticmethod
    def _invert_project(A, eps=1e-10):
        n, m = A.shape
        U, S, V = LINALG.svd(A)
        inv_A = V * np.linalg.pinv(S) * U.T

        r = (S > eps).sum()
        P_R   = U[:, :r] * U[:, :r].T
        P_N   = U[:, r:] * U[:, r:].T

        return inv_A, r, P_R, P_N

    @staticmethod
    def _row_uniq(A):
        """
            A procedure usually performed before linear regression (i.e. solving Ax = y).
            If the matrix A contains repeating rows, it is advisable to combine
            all of them to one row, and the observed value corresponding to that
            row will be the average of the original observations.

            Input:
                A - a 2D NumPy array

            Returns:
                A_unique, P_row

                where A_unique has the same number of columns as A, but with
                unique rows.
                P_row is a matrix that can be used to map the original rows
                to the ones in A_unique (all values in P_row are 0 or 1).
        """
        # convert the rows of A into tuples so we can compare them
        A_tuples = [tuple(A[i,:].flat) for i in xrange(A.shape[0])]
        A_unique = list(sorted(set(A_tuples), reverse=True))

        # create the projection matrix that maps the rows in A to rows in
        # A_unique
        P_col = np.matrix(np.zeros((len(A_unique), len(A_tuples))))

        for j, tup in enumerate(A_tuples):
            # find the indices of the unique row in A_unique which correspond
            # to this original row in A (represented as 'tup')
            i = A_unique.index(tup)
            P_col[i, j] = 1

        return np.matrix(A_unique), P_col

    @staticmethod
    def _col_uniq(A):
        A_unique, P_col = LINALG._row_uniq(A.T)
        return A_unique.T, P_col.T
