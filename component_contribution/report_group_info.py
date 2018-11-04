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


from component_contribution.component_contribution_trainer import ComponentContribution
import numpy as np
import csv, sys

cc = ComponentContribution.init()
group_names = cc.groups_data.GetGroupNames()
dG0_gc = cc.params['dG0_gc']
G = cc.params['G']

Ng = len(group_names)
group_names += map(cc.cids_joined.__getitem__, np.nonzero(G[:, Ng:])[0])

csv_writer = csv.writer(sys.stdout)
csv_writer.writerow(['Group', 'dG0_gc [kJ/mol]'])
csv_writer.writerows(zip(group_names, map(lambda x : '%.2f' % x, dG0_gc.flat)))
