import openbabel, urllib, logging
import chemaxon
import numpy as np
from thermodynamic_constants import R, debye_huckel
from scipy.misc import logsumexp

MIN_PH = 0.0
MAX_PH = 14.0

class Compound(object):
    
    def __init__(self, database, compound_id, inchi,
                 atom_bag, pKas, smiles_pH7, majorMSpH7, nHs, zs):
        self.database = database
        self.compound_id = compound_id
        self.inchi = inchi
        self.atom_bag = atom_bag
        self.pKas = pKas
        self.smiles_pH7 = smiles_pH7
        self.majorMSpH7 = majorMSpH7
        self.nHs = nHs
        self.zs = zs
    
    @staticmethod
    def from_kegg(compound_id):
        return Compound.from_inchi('KEGG', compound_id,
                                   Compound.get_inchi(compound_id))

    @staticmethod
    def from_inchi(database, compound_id, inchi):
        if compound_id == 'C00080':
            # We add an exception for H+ (and put nH = 0) in order to eliminate
            # its effect of the Legendre transform
            return Compound(database, compound_id, inchi,
                            {'H' : 1}, [], None, 0, [0], [0])
        elif compound_id == 'C00087':
            # ChemAxon gets confused with the structure of sulfur
            # (returns a protonated form, [SH-], at pH 7).
            # So we implement it manually here.
            return Compound(database, compound_id, inchi,
                            {'S' : 1, 'e-': 16}, [], 'S', 0, [0], [0])
        elif compound_id == 'C00237':
            # ChemAxon gets confused with the structure of carbon monoxide
            # (returns a protonated form, [CH]#[O+], at pH 7).
            # So we implement it manually here.
            return Compound(database, compound_id, inchi,
                            {'C' : 1, 'O': 1, 'e-': 14}, [], '[C-]#[O+]', 0, [0], [0])
        elif compound_id == 'C00282':
            # ChemAxon gets confused with the structure of hydrogen
            # So we implement it manually here.
            return Compound(database, compound_id, inchi,
                            {'H' : 2, 'e-': 2}, [], None, 0, [2], [0])
        elif compound_id == 'C01353':
            # When given the structure of carbonic acid, ChemAxon returns the
            # pKas for CO2(tot), i.e. it assumes the non-hydrated CO2 species is
            # one of the pseudoisomers, and the lower pKa value is 6.05 instead of 
            # 3.78. Here, we introduce a new "KEGG" compound that will represent
            # pure bicarbonate (without CO2(sp)) and therefore plug in the pKa
            # values from Alberty's book.
            return Compound(database, compound_id, inchi,
                            {'C': 1, 'H': 1, 'O': 3, 'e-': 32}, [10.33, 3.43],
                            'OC(=O)[O-]', 1, [0, 1, 2], [-2, -1, 0])            
        # Metal Cations get multiple pKa values from ChemAxon, which is
        # obviously a bug. We override the important ones here:
        elif compound_id == 'C00076': # Ca2+
            return Compound(database, compound_id, inchi,
                            {'Ca' : 1, 'e-': 18}, [], '[Ca++]', 0, [0], [2])
        elif compound_id == 'C00238': # K+
            return Compound(database, compound_id, inchi,
                            {'K' : 1, 'e-': 18}, [], '[K+]', 0, [0], [1])
        elif compound_id == 'C00305': # Mg2+
            return Compound(database, compound_id, inchi,
                            {'Mg' : 1, 'e-': 10}, [], '[Mg++]', 0, [0], [2])
        elif compound_id == 'C14818': # Fe2+
            return Compound(database, compound_id, inchi,
                            {'Fe' : 1, 'e-': 24}, [], '[Fe++]', 0, [0], [2])
        elif compound_id == 'C14819': # Fe3+
            return Compound(database, compound_id, inchi,
                            {'Fe' : 1, 'e-': 23}, [], '[Fe+++]', 0, [0], [3])
        elif compound_id == 'C00138': # ferredoxin(red)
            return Compound(database, compound_id, inchi,
                            {'Fe' : 1, 'e-': 26}, [], None, 0, [0], [0])
        elif compound_id == 'C00139': # ferredoxin(ox)
            return Compound(database, compound_id, inchi,
                            {'Fe' : 1, 'e-': 25}, [], None, 0, [0], [1])
        elif inchi is None:
            # If the compound has no explicit structure, we assume that it has 
            # no proton dissociations in the relevant pH range
            return Compound(database, compound_id, inchi,
                            {}, [], None, 0, [0], [0])
        
        # Otherwise, we use ChemAxon's software to get the pKas and the 
        # properties of all microspecies

        try:
            pKas, major_ms_smiles = chemaxon.GetDissociationConstants(inchi)
            major_ms_smiles = Compound.smiles2smiles(major_ms_smiles)
            pKas = sorted([pka for pka in pKas if pka > MIN_PH and pka < MAX_PH], reverse=True)
        except chemaxon.ChemAxonError:
            logging.warning('chemaxon failed to find pKas for this molecule: ' + inchi)
            # use the original InChI to get the parameters (i.e. assume it 
            # represents the major microspecies at pH 7)
            major_ms_smiles = Compound.inchi2smiles(inchi)
            pKas = []
        
        if major_ms_smiles:
            atom_bag, major_ms_charge = chemaxon.GetAtomBagAndCharge(major_ms_smiles)
            major_ms_nH = atom_bag.get('H', 0)
        else:
            atom_bag = {}
            major_ms_charge = 0
            major_ms_nH = 0

        n_species = len(pKas) + 1
        if pKas == []:
            majorMSpH7 = 0
        else:
            majorMSpH7 = len([1 for pka in pKas if pka > 7])
            
        nHs = []
        zs = []

        for i in xrange(n_species):
            zs.append((i - majorMSpH7) + major_ms_charge)
            nHs.append((i - majorMSpH7) + major_ms_nH)
        
        return Compound(database, compound_id, inchi,
                        atom_bag, pKas, major_ms_smiles, majorMSpH7, nHs, zs)

    def to_json_dict(self):
        return {'database' : self.database,
                'compound_id' : self.compound_id,
                'inchi' : self.inchi,
                'atom_bag' : self.atom_bag,
                'pKas' : self.pKas,
                'smiles_pH7' : self.smiles_pH7,
                'majorMSpH7' : self.majorMSpH7,
                'nHs' : self.nHs,
                'zs' : self.zs}
    
    @staticmethod
    def from_json_dict(d):
        return Compound(d['database'], d['compound_id'], d['inchi'], d['atom_bag'],
                        d['pKas'], d['smiles_pH7'], d['majorMSpH7'],
                        d['nHs'], d['zs'])

    @staticmethod
    def get_inchi(compound_id):
        s_mol = urllib.urlopen('http://rest.kegg.jp/get/cpd:%s/mol' % compound_id).read()
        return Compound.mol2inchi(s_mol)

    @staticmethod
    def mol2inchi(s):
        openbabel.obErrorLog.SetOutputLevel(-1)

        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('mol', 'inchi')
        conv.AddOption("F", conv.OUTOPTIONS)
        conv.AddOption("T", conv.OUTOPTIONS)
        conv.AddOption("x", conv.OUTOPTIONS, "noiso")
        conv.AddOption("w", conv.OUTOPTIONS)
        obmol = openbabel.OBMol()
        if not conv.ReadString(obmol, str(s)):
            return None
        inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if inchi == '':
            return None
        else:
            return inchi

    @staticmethod
    def inchi2smiles(inchi):
        openbabel.obErrorLog.SetOutputLevel(-1)
        
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('inchi', 'smiles')
        #conv.AddOption("F", conv.OUTOPTIONS)
        #conv.AddOption("T", conv.OUTOPTIONS)
        #conv.AddOption("x", conv.OUTOPTIONS, "noiso")
        #conv.AddOption("w", conv.OUTOPTIONS)
        obmol = openbabel.OBMol()
        conv.ReadString(obmol, str(inchi))
        smiles = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if smiles == '':
            return None
        else:
            return smiles
            
    @staticmethod
    def smiles2smiles(smiles_in):
        openbabel.obErrorLog.SetOutputLevel(-1)
        
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('smiles', 'smiles')
        #conv.AddOption("F", conv.OUTOPTIONS)
        #conv.AddOption("T", conv.OUTOPTIONS)
        #conv.AddOption("x", conv.OUTOPTIONS, "noiso")
        #conv.AddOption("w", conv.OUTOPTIONS)
        obmol = openbabel.OBMol()
        conv.ReadString(obmol, str(smiles_in))
        smiles_out = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if smiles_out == '':
            return None
        else:
            return smiles_out
    @staticmethod
    def smiles2inchi(smiles):
        openbabel.obErrorLog.SetOutputLevel(-1)
        
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('smiles', 'inchi')
        conv.AddOption("F", conv.OUTOPTIONS)
        conv.AddOption("T", conv.OUTOPTIONS)
        conv.AddOption("x", conv.OUTOPTIONS, "noiso")
        conv.AddOption("w", conv.OUTOPTIONS)
        obmol = openbabel.OBMol()
        conv.ReadString(obmol, str(smiles))
        inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if inchi == '':
            return None
        else:
            return inchi

    def __str__(self):
        return "%s\nInChI: %s\npKas: %s\nmajor MS: nH = %d, charge = %d" % \
            (self.compound_id, self.inchi, ', '.join(['%.2f' % p for p in self.pKas]),
             self.nHs[self.majorMSpH7], self.zs[self.majorMSpH7])

    def _dG0_prime_vector(self, pH, I, T):
        """
            Calculates the difference in kJ/mol between dG'0 and 
            the dG0 of the MS with the least hydrogens (dG0[0])
            
            Returns:
                dG'0 - dG0[0]
        """
        if self.inchi is None:
            return 0
        elif self.pKas == []:
            dG0s = np.zeros((1, 1))
        else:
            dG0s = -np.cumsum([0] + self.pKas) * R * T * np.log(10)
            dG0s = dG0s
        DH = debye_huckel((I, T))
        
        # dG0' = dG0 + nH * (R T ln(10) pH + DH) - charge^2 * DH
        pseudoisomers = np.vstack([dG0s, np.array(self.nHs), np.array(self.zs)]).T
        dG0_prime_vector = pseudoisomers[:, 0] + \
                           pseudoisomers[:, 1] * (R * T * np.log(10) * pH + DH) - \
                           pseudoisomers[:, 2]**2 * DH
        return dG0_prime_vector
        
    def _transform(self, pH, I, T):

        return -R * T * logsumexp(self._dG0_prime_vector(pH, I, T) / (-R * T))

    def _ddG(self, i_from, i_to, T):
        """
            Calculates the difference in kJ/mol between two MSs.
            
            Returns:
                dG0[i_to] - dG0[i_from]
        """
        if not (0 <= i_from <= len(self.pKas)):
            raise ValueError('MS index is out of bounds: 0 <= %d <= %d' % (i_from, len(self.pKas)))

        if not (0 <= i_to <= len(self.pKas)):
            raise ValueError('MS index is out of bounds: 0 <= %d <= %d' % (i_to, len(self.pKas)))

        if i_from == i_to:
            return 0
        elif i_from < i_to:
            return sum(self.pKas[i_from:i_to]) * R * T * np.log(10)
        else:
            return -sum(self.pKas[i_to:i_from]) * R * T * np.log(10)

    def transform(self, i, pH, I, T):
        """
            Returns the difference in kJ/mol between dG'0 and the dG0 of the 
            MS with index 'i'.
            
            Returns:
                (dG'0 - dG0[0]) + (dG0[0] - dG0[i])  = dG'0 - dG0[i]
        """
        return self._transform(pH, I, T) + self._ddG(0, i, T)

    def transform_pH7(self, pH, I, T):
        """
            Returns the transform for the major MS in pH 7
        """
        return self.transform(self.majorMSpH7, pH, I, T)

    def transform_neutral(self, pH, I, T):
        """
            Returns the transform for the MS with no charge
        """
        try:
            return self.transform(self.zs.index(0), pH, I, T)
        except ValueError:
            raise ValueError("The compound (%s) does not have a microspecies with 0 charge"
                             % self.compound_id)

    def get_species(self, major_ms_dG0_f, T):
        """
            Given the chemical formation energy of the major microspecies,
            uses the pKa values to calculate the chemical formation energies
            of all other species, and returns a list of dictionaries with
            all the relevant data: dG0_f, nH, nMg, z (charge)
        """
        for i, (nH, z) in enumerate(zip(self.nHs, self.zs)):
            dG0_f = major_ms_dG0_f + self._ddG(i, self.majorMSpH7, T)
            d = {'phase': 'aqueous', 'dG0_f': np.round(dG0_f, 2),
                 'nH': nH, 'z': z, 'nMg': 0}
            yield d
        
if __name__ == '__main__':
    import sys, json
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    from python.compound_cacher import CompoundCacher, CompoundEncoder
    from python.molecule import Molecule, OpenBabelError
    ccache = CompoundCacher(cache_fname=None)

    for compound_id in ['C00087', 'C00282', 'C00237']:
        comp = Compound.from_kegg(compound_id)
        try:
            mol = Molecule.FromInChI(str(comp.inchi))
            sys.stderr.write('%s : formula = %s, nE = %s' % 
                             (str(comp.inchi), mol.GetFormula(), mol.GetNumElerctons()))
        except OpenBabelError:
            pass
        ccache.add(comp)
        sys.stderr.write('\ncompound id = %s, nH = %s, z = %s, pKa = %s, bag = %s\n\n\n' % 
                         (compound_id, str(comp.nHs), str(comp.zs), str(comp.pKas), str(comp.atom_bag)))
    
    ccache.dump()
