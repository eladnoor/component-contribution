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


from collections import defaultdict

import numpy as np
import pandas as pd
import pybel
from pkg_resources import resource_stream
from scipy.special import logsumexp

from component_contribution import chemaxon
from component_contribution.databases import databases
from component_contribution.mol_utils import atom_bag_and_charge
from component_contribution.thermodynamic_constants import R, debye_huckel


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

with resource_stream('component_contribution',
                     '/data/compound_additions.csv') as fp:
    COMPOUND_ADDITIONS = pd.read_csv(fp).set_index('cid')


class MicroSpecie(object):
    """
        A class for storing microspecie data (i.e. data about a single
        protonation state of a single compound)

        ddG_over_T is equal to the sum of pKa values (between this MS and
        the major MS, times Rlog10). If we multiply this by T, we will get
        the relative ddG between this MS and the major MS, in kJ/mol.
    """
    Rlog10 = R * np.log(10)

    def __init__(self, inchi_key, charge, number_of_protons,
                 ddG_over_T, is_major):
        assert type(inchi_key) == str
        assert type(charge) in [int, np.int64]
        assert type(number_of_protons) in [int, np.int64]
        #assert type(ddG_over_T) == float
        assert type(is_major) == bool

        self.inchi_key = inchi_key
        self.charge = charge
        self.number_of_protons = number_of_protons
        self.ddG_over_T = ddG_over_T
        self.is_major = is_major

    def transform(self, p_h, ionic_strength, temperature):
        """
            Use the Legendre transform to convert the ddG_over_T to the
            difference in the transformed energies of this MS and the
            major MS
        """
        DH = debye_huckel(ionic_strength, temperature)
        return (self.ddG_over_T * temperature +
                self.number_of_protons * temperature * self.Rlog10 * p_h +
                self.number_of_protons * DH -
                self.charge**2 * DH)

    @classmethod
    def from_p_kas(cls, inchi_key, major_ms, charges, number_of_protons,
                   p_kas):
        species = []
        for i, (z, nH) in enumerate(zip(charges, number_of_protons)):
            if i == major_ms:
                ms = cls(inchi_key, z, nH, 0.0, True)
            elif i < major_ms:
                ms = cls(inchi_key, z, nH,
                         sum(p_kas[i:major_ms]) * MicroSpecie.Rlog10,
                         False)
            elif i > major_ms:
                ms = cls(inchi_key, z, nH,
                         -sum(p_kas[major_ms:i]) * MicroSpecie.Rlog10,
                         False)
            else:
                raise Exception('Internal error')

            species.append(ms)
        return species

    def to_dict(self):
        return {'inchi_key': self.inchi_key,
                'charge': self.charge,
                'number_of_protons': self.number_of_protons,
                'ddG_over_T': self.ddG_over_T,
                'is_major': self.is_major}

    @classmethod
    def from_dict(cls, d):
        return cls(d['inchi_key'], d['charge'], d['number_of_protons'],
                   d['ddG_over_T'], d['is_major'])


class Compound(object):

    def __init__(self, inchi_key, inchi, smiles, atom_bag, species,
                 compound_id=None):
        assert type(atom_bag) in [dict, defaultdict]
        assert type(species) == list

        self.inchi_key = inchi_key
        self.inchi = inchi
        self.smiles = smiles
        self.atom_bag = atom_bag
        self.species = species
        self.compound_id = compound_id

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

            return cls(inchi_key=compound_id, inchi='', atom_bag={},
                       species=[], smiles=None, compound_id=compound_id)

        if compound_id in COMPOUND_EXCEPTIONS:
            atom_bag, p_kas, major_ms_smiles, major_microspecies, \
                charges, number_of_protons = COMPOUND_EXCEPTIONS[compound_id]
        else:
            if compute_pkas:
                p_kas, major_ms_smiles = \
                    chemaxon.get_dissociation_constants(inchi)
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

        species = MicroSpecie.from_p_kas(inchi_key, major_microspecies,
                                         charges, number_of_protons,
                                         p_kas)

        return cls(inchi_key, inchi, atom_bag, species, major_ms_smiles,
                   compound_id=compound_id)

    def __str__(self):
        return "%s\nInChI: %s\nmajor MS: nH = %d, charge = %d" % \
            (self.compound_id, self.inchi,
             self.number_of_protons[self.major_microspecies],
             self.charges[self.major_microspecies])

    def ddG_primes(self, p_h, ionic_strength, temperature):
        RT = R * temperature
        return [-ms.transform(p_h, ionic_strength, temperature) / RT
                for ms in self.species]

    def transform(self, p_h, ionic_strength, temperature):
        if self.species == []:
            return 0.0
        return -R * temperature * logsumexp(
                self.ddG_primes(p_h, ionic_strength, temperature))

    def get_species(self, major_ms_dG0_f, temperature):
        """
        Given the chemical formation energy of the major microspecies,
        uses the pKa values to calculate the chemical formation energies
        of all other species, and returns a list of dictionaries with
        all the relevant data: dG0_f, nH, nMg, z (charge)
        """
        for s in self.species:
            dG0_f = major_ms_dG0_f + s.ddG_over_T * temperature
            d = {'phase': 'aqueous',
                 'dG0_f': np.round(dG0_f, 2),
                 'number_of_protons': s.number_of_protons,
                 'charge': s.charge,
                 'number_of_magnesium': 0}
            yield d
