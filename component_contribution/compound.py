# The MIT License (MIT)
#
# Copyright (c) 2013 The Weizmann Institute of Science.
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

import numpy as np
import pybel
from scipy.special import logsumexp
from pkg_resources import resource_stream
import pandas as pd
from collections import defaultdict

from component_contribution.databases import databases
from component_contribution.mol_utils import atom_bag_and_charge
from component_contribution.thermodynamic_constants import R, debye_huckel
from component_contribution import chemaxon

MIN_PH = 0.0
MAX_PH = 14.0


COMPOUND_EXCEPTIONS = {
    # We add an exception for H+ (and put nH = 0) in order to
    # eliminate its effect of the Legendre transform.
    'KEGG:C00080': ({'H': 1}, [], None, 0, [0], [0]),
    # ChemAxon gets confused with the structure of sulfur
    #  (returns a protonated form, [SH-], at pH 7).
    'KEGG:C00087': ({'S': 1, 'e-': 16}, [], 'S', 0, [0], [0]),
    # ChemAxon gets confused with the structure of carbon
    # monoxide (protonated form, [CH]#[O+], at pH 7).
    'KEGG:C00237': ({'C': 1, 'O': 1, 'e-': 14}, [], '[C-]#[O+]',
                    0, [0], [0]),
    'KEGG:C00282': ({'H': 2, 'e-': 2}, [], None, 0, [2], [0]),
    # When given the structure of carbonic acid,
    # ChemAxon returns the pKas for CO2(tot), i.e.
    # it assumes the non-hydrated CO2 species is one
    # of the pseudoisomers, and the lower pKa value
    # is 6.05 instead of 3.78. Here, we introduce
    # a new "KEGG" compound that will represent
    # pure bicarbonate (without CO2(sp)) and
    # therefore plug in the pKa values from
    # Alberty's book.
    'KEGG:C01353': ({'C': 1, 'H': 1, 'O': 3, 'e-': 32}, [10.33, 3.43],
                    'OC(=O)[O-]', 1, [0, 1, 2], [-2, -1, 0]),
    # Metal Cations get multiple pKa values from ChemAxon, which is
    # obviously a bug. We override the important ones here:
    # Ca2+
    'KEGG:C00076': ({'Ca': 1, 'e-': 18}, [], '[Ca++]', 0, [0], [2]),
    # K+
    'KEGG:C00238': ({'K': 1, 'e-': 18}, [], '[K+]', 0, [0], [1]),
    # Mg2+
    'KEGG:C00305': ({'Mg': 1, 'e-': 10}, [], '[Mg++]', 0, [0], [2]),
    # Fe2+
    'KEGG:C14818': ({'Fe': 1, 'e-': 24}, [], '[Fe++]', 0, [0], [2]),
    # Fe3+
    'KEGG:C14819': ({'Fe': 1, 'e-': 23}, [], '[Fe+++]', 0, [0], [3]),
    # ferredoxin(red)
    'KEGG:C00138': ({'Fe': 1, 'e-': 26}, [], None, 0, [0], [0]),
    # ferredoxin(ox)
    'KEGG:C00139': ({'Fe': 1, 'e-': 25}, [], None, 0, [0], [1])
}

COMPOUND_ADDITIONS = pd.read_csv(
    resource_stream('component_contribution', '/data/compound_additions.csv')).set_index('cid')


class Compound(object):

    def __init__(self, inchi_key, inchi, atom_bag, p_kas, smiles,
                 major_microspecies, number_of_protons, charges,
                 compound_id=None):
        assert type(atom_bag) in [dict, defaultdict]
        assert type(p_kas) == list
        if not type(major_microspecies) in [int, np.int64]:
            raise AssertionError('major_ms is not an integer, but a %s' %
                                 type(major_microspecies))
        assert type(number_of_protons) == list
        assert type(charges) == list
        
        self.inchi_key = inchi_key
        self.compound_id = compound_id
        self.name = compound_id
        self.inchi = inchi
        self.atom_bag = atom_bag
        self.p_kas = p_kas
        self.smiles = smiles
        self.major_microspecies = major_microspecies
        self.number_of_protons = number_of_protons
        self.charges = charges

    @classmethod
    def get(cls, compound_id, compute_pkas=True):
        if compound_id in COMPOUND_ADDITIONS.index:
            inchi = COMPOUND_ADDITIONS.at[compound_id, 'inchi']
            molecule = pybel.readstring("inchi", inchi)
        else:
            try:
                database, accession = compound_id.split(":", 1)
            except ValueError:
                # assume by default that this is a KEGG compound ID
                database = 'KEGG'
                accession = compound_id
                
            molecule = databases.get_molecule(database, accession)
        return cls.from_molecule(compound_id, molecule, compute_pkas)

    @classmethod
    def from_molecule(cls, compound_id, molecule, compute_pkas=True):
        if molecule is not None:
            inchi = molecule.write("inchi").strip()
            inchi_key = molecule.write("inchikey").strip()
        else:
            inchi = ''
            inchi_key = ''
            
        if not inchi_key:
            # probably a compound without an explicit chemical structure
            # or formula. we thus use the compound_id instead of the InChIKey
            # TODO: find a way to map these compounds between different databases
                        
            return cls(inchi_key=compound_id, inchi='', atom_bag={}, p_kas=[],
                       smiles=None, major_microspecies=0,
                       number_of_protons=[], charges=[],
                       compound_id=compound_id)
        
        if compound_id in COMPOUND_EXCEPTIONS:
            return cls(inchi_key, inchi, *COMPOUND_EXCEPTIONS[compound_id],
                       compound_id=compound_id)

        if compute_pkas:
            p_kas, major_ms_smiles = chemaxon.get_dissociation_constants(inchi)
            molecule = pybel.readstring("smi", major_ms_smiles)
            p_kas = sorted([pka for pka in p_kas if MIN_PH < pka < MAX_PH],
                           reverse=True)
        else:
            p_kas = []
            molecule.addh()
            major_ms_smiles = molecule.write('smi')

        atom_bag, major_ms_charge = atom_bag_and_charge(molecule)
        _number_of_protons = atom_bag.get('H', 0)

        n_species = len(p_kas) + 1

        if not p_kas:
            major_microspecies = 0
        else:
            major_microspecies = len([1 for pka in p_kas if pka > 7])

        number_of_protons = []
        charges = []

        for i in range(n_species):
            charges.append((i - major_microspecies) + major_ms_charge)
            number_of_protons.append(
                (i - major_microspecies) + _number_of_protons)

        return cls(inchi_key, inchi, atom_bag, p_kas, major_ms_smiles,
                   major_microspecies, number_of_protons, charges,
                   compound_id=compound_id)

    def to_json_dict(self):
        return {'inchi_key': self.inchi_key,
                'inchi': self.inchi,
                'atom_bag': self.atom_bag,
                'p_kas': self.p_kas,
                'smiles': self.smiles,
                'major_microspecies': self.major_microspecies,
                'number_of_protons': self.number_of_protons,
                'charges': self.charges,
                'compound_id': self.compound_id}

    @classmethod
    def from_json_dict(cls, data):
        return cls(
            data['inchi_key'], data['inchi'], data['atom_bag'], data['p_kas'],
            data['smiles'], data['major_microspecies'],
            data['number_of_protons'], data['charges'], data['compound_id'])

    def __str__(self):
        return "%s\nInChI: %s\npKas: %s\nmajor MS: nH = %d, charge = %d" % \
            (self.compound_id, self.inchi, ', '.join(
                ['%.2f' % p for p in self.p_kas]),
             self.number_of_protons[self.major_microspecies],
             self.charges[self.major_microspecies])

    def _dG0_prime_vector(self, p_h, ionic_strength, temperature):
        """
            Calculates the difference in kJ/mol between dG'0 and
            the dG0 of the MS with the least hydrogens (dG0[0])

            Returns:
                dG'0 - dG0[0]
        """
        if not self.inchi:
            return 0
        elif not self.p_kas:
            dG0s = np.zeros((1, 1))
        else:
            dG0s = -np.cumsum([0] + self.p_kas) * R * temperature * np.log(10)
            dG0s = dG0s
        DH = debye_huckel(ionic_strength, temperature)

        # dG0' = dG0 + nH * (R T ln(10) pH + DH) - charge^2 * DH
        pseudoisomers = np.vstack([dG0s, np.array(self.number_of_protons),
                                   np.array(self.charges)]).T
        dG0_prime_vector = pseudoisomers[:, 0] + \
            pseudoisomers[:, 1] * (R * temperature * np.log(10) * p_h + DH) - \
            pseudoisomers[:, 2]**2 * DH
        return dG0_prime_vector

    def _transform(self, p_h, ionic_strength, temperature):
        return -R * temperature * logsumexp(
            self._dG0_prime_vector(p_h, ionic_strength, temperature) /
            (-R * temperature))

    def _ddG(self, i_from, i_to, temperature):
        """
            Calculates the difference in kJ/mol between two MSs.

            Returns:
                dG0[i_to] - dG0[i_from]
        """
        if not (0 <= i_from <= len(self.p_kas)):
            raise ValueError('MS index is out of bounds: 0 <= %d <= %d' % (
                i_from, len(self.p_kas)))

        if not (0 <= i_to <= len(self.p_kas)):
            raise ValueError('MS index is out of bounds: 0 <= %d <= %d' % (
                i_to, len(self.p_kas)))

        if i_from == i_to:
            return 0
        elif i_from < i_to:
            return sum(self.p_kas[i_from:i_to]) * R * temperature * np.log(10)
        else:
            return -sum(self.p_kas[i_to:i_from]) * R * temperature * np.log(10)

    def transform(self, index, p_h, ionic_strength, temperature):
        """
        Return the difference in kJ/mol between dG'0 and the dG0.

        For a given microspecies at a given index the difference is given by:
                (dG'0 - dG0[0]) + (dG0[0] - dG0[i])  = dG'0 - dG0[i]

        Parameters
        ----------
        index : int
            The index of of the microspecies.
        p_h : float
            The pH value.
        ionic_strength : float
            The ionic strength of the solution.
        temperature : float
            Temperature in Kelvins.

        Returns
        -------
        float
            The difference in kJ/mol.

        """
        return self._transform(
            p_h, ionic_strength, temperature) + self._ddG(0, index, temperature)

    def transform_p_h_7(self, p_h, ionic_strength, temperature):
        """
        Parameters
        ----------
        p_h : float
            The pH value.
        ionic_strength : float
            The ionic strength of the solution.
        temperature : float
            Temperature in Kelvins.

        Returns
        -------
        transform : float
            The transform for the major microspecies at pH 7
        """
        return self.transform(
            self.major_microspecies, p_h, ionic_strength, temperature)

    def transform_neutral(self, p_h, ionic_strength, temperature):
        """
            Returns the transform for the MS with no charge
        """
        try:
            return self.transform(
                self.charges.index(0), p_h, ionic_strength, temperature)
        except ValueError:
            raise ValueError(
                "The compound (%s) does not have a microspecies with 0 "
                "charge." % (self.compound_id,))

    def get_species(self, major_ms_dG0_f, temperature):
        """
        Given the chemical formation energy of the major microspecies,
        uses the pKa values to calculate the chemical formation energies
        of all other species, and returns a list of dictionaries with
        all the relevant data: dG0_f, nH, nMg, z (charge)
        """
        for i, (num_protons, charge) in enumerate(zip(
                self.number_of_protons, self.charges)):
            dG0_f = major_ms_dG0_f + self._ddG(
                i, self.major_microspecies, temperature)
            d = {'phase': 'aqueous', 'dG0_f': np.round(dG0_f, 2),
                 'number_of_protons': num_protons, 'charge': charge,
                 'number_of_magnesium': 0}
            yield d

