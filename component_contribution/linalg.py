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

from numpy import array, zeros, fill_diagonal
from scipy.linalg import svd
from numpy.linalg import pinv

class LINALG(object):

    @staticmethod
    def _zero_pad_S(S, cids_orig, cids_joined):
        """
            takes a stoichiometric matrix with a given list of IDs 'cids' and adds
            0-rows so that the list of IDs will be 'cids_joined'
        """
        if not set(cids_orig).issubset(cids_joined):
            raise Exception('The full list is missing some IDs in "cids"')

        full_S = zeros((len(cids_joined), S.shape[1]))
        for i, cid in enumerate(cids_orig):
            S_row = S[i, :]
            full_S[cids_joined.index(cid), :] = S_row

        return full_S

    @staticmethod
    def _invert_project(A, eps=1e-10):
        n, m = A.shape
        U, s, Vh = svd(A, full_matrices=True)
        S = zeros(A.shape)
        fill_diagonal(S, s)
        inv_A = Vh.T @ pinv(S) @ U.T

        r = (S > eps).sum()
        P_R   = U[:, :r] @ U[:, :r].T
        P_N   = U[:, r:] @ U[:, r:].T

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
        A_tuples = [tuple(A[i,:].flat) for i in range(A.shape[0])]
        A_unique = list(sorted(set(A_tuples), reverse=True))

        # create the projection matrix that maps the rows in A to rows in
        # A_unique
        P_col = zeros((len(A_unique), len(A_tuples)))

        for j, tup in enumerate(A_tuples):
            # find the indices of the unique row in A_unique which correspond
            # to this original row in A (represented as 'tup')
            i = A_unique.index(tup)
            P_col[i, j] = 1

        return array(A_unique, ndmin=2), P_col

    @staticmethod
    def _col_uniq(A):
        A_unique, P_col = LINALG._row_uniq(A.T)
        return A_unique.T, P_col.T

