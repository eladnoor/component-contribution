import openbabel, urllib, re, logging
import chemaxon
import numpy as np
from thermodynamic_constants import R, debye_huckel
from scipy.misc import logsumexp

MIN_PH = 0.0
MAX_PH = 14.0

class Compound(object):
    
    _obElements = openbabel.OBElementTable()

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
    def from_kegg(cid):
        inchi = Compound.get_inchi(cid)
        database = 'KEGG'
        compound_id = 'C%05d' % cid
        if cid == 80:
            # We add an exception for H+ (and put nH = 0) in order to eliminate
            # its effect of the Legendre transform
            return Compound(database, compound_id, inchi,
                            {'H' : 1}, [], None, 0, [0], [0])
        elif cid == 237:
            # ChemAxon gets confused with the structure of carbon monoxide
            # (returns a protonated form, [CH]#[O+] at pH 7).
            # So we implement it manually here.
            return Compound(database, compound_id, inchi,
                            {'C' : 1, 'O': 1}, [], '[C-]#[O+]', 0, [0], [0])
        elif inchi is None:
            # If the compound has no explicit structure, we assume that it has 
            # no proton dissociations in the relevant pH range
            return Compound(database, compound_id, inchi,
                            {}, [], None, 0, [0], [0])
        else:
            # Otherwise, we use ChemAxon's software to get the pKas and the 
            # properties of all microspecies
            return Compound.from_inchi(database, compound_id, inchi)

    @staticmethod
    def from_inchi(database, compound_id, inchi):
        try:
            pKas, major_ms_smiles = chemaxon.GetDissociationConstants(inchi)
            pKas = sorted([pka for pka in pKas if pka > MIN_PH and pka < MAX_PH], reverse=True)
            atom_bag, major_ms_charge = chemaxon.GetAtomBagAndCharge(major_ms_smiles)
        except chemaxon.ChemAxonError:
            logging.warning('chemaxon failed to find pKas for this molecule: ' + inchi)
            # use the original InChI to get the parameters (i.e. assume it 
            # represents the major microspecies at pH 7)
            major_ms_smiles = Compound.inchi2smiles(inchi)
            pKas = []
            atom_bag, major_ms_charge = chemaxon.GetAtomBagAndCharge(inchi)

        major_ms_nH = atom_bag.get('H', 0)

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
                'id' : self.compound_id,
                'inchi' : self.inchi,
                'atom_bag' : self.atom_bag,
                'pKas' : self.pKas,
                'smiles_pH7' : self.smiles_pH7,
                'majorMSpH7' : self.majorMSpH7,
                'nHs' : self.nHs,
                'zs' : self.zs}
    
    @staticmethod
    def from_json_dict(d):
        return Compound(d['database'], d['id'], d['inchi'], d['atom_bag'],
                        d['pKas'], d['smiles_pH7'], d['majorMSpH7'],
                        d['nHs'], d['zs'])

    @staticmethod
    def get_inchi(cid):
        s_mol = urllib.urlopen('http://rest.kegg.jp/get/cpd:C%05d/mol' % cid).read()
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
        if not conv.ReadString(obmol, s):
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
        conv.ReadString(obmol, inchi)
        smiles = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if smiles == '':
            return None
        else:
            return smiles

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
        conv.ReadString(obmol, smiles)
        inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if inchi == '':
            return None
        else:
            return inchi

    def __str__(self):
        return "%s\nInChI: %s\npKas: %s\nmajor MS: nH = %d, charge = %d" % \
            (self.compound_id, self.inchi, ', '.join(['%.2f' % p for p in self.pKas]),
             self.nHs[self.majorMSpH7], self.zs[self.majorMSpH7])
    
    def transform(self, pH, I, T):
        """
            Returns the difference between the dG0 of the major microspecies at pH 7
            and the transformed dG0' (in kJ/mol)
        """
        ddG0 = sum(self.pKas[:self.majorMSpH7]) * R * T * np.log(10)
        return self._transform(pH, I, T) + ddG0

    def transform_neutral(self, pH, I, T):
        """
            Returns the difference between the dG0 of microspecies with the 0 charge
            and the transformed dG0' (in kJ/mol)
        """
        try:
            MS_ind = self.zs.index(0)
        except ValueError:
            raise ValueError("The compound (%s) does not have a microspecies with 0 charge"
                             % self.compound_id)
        
        # calculate the difference between the microspecies with the least hydrogens
        # and the microspecies with 0 charge 
        ddG0 = sum(self.pKas[:MS_ind]) * R * T * np.log(10)
        return self._transform(pH, I, T) + ddG0

    def _transform(self, pH, I, T):
        """
            Returns the difference between the dG0 of microspecies with the least hydrogens
            and the transformed dG0' (in kJ/mol)
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

        return -R * T * logsumexp(dG0_prime_vector / (-R * T))

if __name__ == '__main__':
    import sys
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)

    for cid in [138, 282, 237, 87, 32, 1137, 2191, 15670, 15672]:
        comp = Compound.from_kegg(cid)
        sys.stderr.write('\ncompound id = C%05d, nH = %s, z = %s\n\n\n' % 
                         (cid, str(comp.nHs), str(comp.zs)))
