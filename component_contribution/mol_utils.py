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

from __future__ import absolute_import

from collections import defaultdict

import openbabel
import pybel
from six import iteritems


PERIODIC_TABLE = _obElements = openbabel.OBElementTable()


def atom_bag_and_charge(molecule):
    """
    Compute the atom bag and the formal charge of a molecule.

    The formal charge is calculated by summing the formal charge of each atom
    in the molecule.

    Parameters
    ----------
    molecule : pybel.Molecule
        A molecule object.

    Returns
    -------
    dict
        A dictionary of atom counts.
    int
        The formal charge of the molecule.

    """

    assert isinstance(molecule, pybel.Molecule)

    atom_bag = defaultdict(int)
    formal_charge = 0
    
    # Make all hydrogens explicit so we can properly count them
    molecule = molecule.clone
    molecule.addh()
    
    for atom in molecule.atoms:
        assert isinstance(atom, pybel.Atom)
        symbol = PERIODIC_TABLE.GetSymbol(atom.atomicnum)
        atom_bag[symbol] += 1
        formal_charge += atom.formalcharge

    n_protons = sum(c * PERIODIC_TABLE.GetAtomicNum(elem)
                    for (elem, c) in iteritems(atom_bag))

    atom_bag['e-'] = n_protons - formal_charge

    return atom_bag, formal_charge


def remove_atoms(molecule, indices):
    """
    Removes atoms from a molecule. Modifications are applied in place.

    Parameters
    ----------

    molecule : pybel.Molecule
        A molecule object.
    indices : list
        A list of indices to remove.

    Returns
    -------
    molecule : pybel.Molecule
        The molecule object modified.
    """
    assert isinstance(molecule, pybel.Molecule)
    molecule.OBMol.BeginModify()
    for i in sorted(indices, reverse=True):
        molecule.OBMol.DeleteAtom(molecule.OBMol.GetAtom(i+1))
    molecule.OBMol.EndModify()
    return molecule
