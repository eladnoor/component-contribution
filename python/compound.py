import openbabel, urllib, re, string, json, logging
from chemaxon import GetDissociationConstants, ChemAxonError
import numpy as np
MIN_PH = 0.0
MAX_PH = 14.0

class Compound(object):
    
    _obElements = openbabel.OBElementTable()

    def __init__(self, database, compound_id, inchi, pKas, majorMSpH7, nHs, zs):
        self.database = database
        self.compound_id = compound_id
        self.inchi = inchi
        self.pKas = pKas
        self.majorMSpH7 = majorMSpH7
        self.nHs = nHs
        self.zs = zs
    
    @staticmethod
    def from_kegg(cid):
        inchi = Compound.get_inchi(cid)
        if inchi == '':
            pKas = np.zeros((1, 1), dtype=float)
            majorMSpH7 = np.ones((1, 1), dtype=int)
            nHs = np.zeros((1, 1), dtype=int)
            zs = np.zeros((1, 1), dtype=int)
        else:
            pKas, majorMSpH7, nHs, zs = Compound.get_species_pka(inchi)
        return Compound('KEGG', 'C%05d' % cid, inchi,
                            pKas, majorMSpH7, nHs, zs)

    def to_json_dict(self):
        return {'database' : self.database,
                'id' : self.compound_id,
                'inchi' : self.inchi,
                'pKas' : self.pKas.tolist(),
                'majorMSpH7' : list(self.majorMSpH7.flat),
                'nHs' : list(self.nHs.flat),
                'zs' : list(self.zs.flat)}
    
    @staticmethod
    def from_json_dict(d):
        pKas       = np.array(d['pKas'],       dtype=float)
        majorMSpH7 = np.array(d['majorMSpH7'], dtype=int)
        nHs        = np.array(d['nHs'],        dtype=int)
        zs         = np.array(d['zs'],         dtype=int)
        return Compound(d['database'], d['id'], d['inchi'],
                            pKas, majorMSpH7, nHs, zs)

    @staticmethod
    def get_inchi(cid):
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('mol', 'inchi')
        #conv.AddOption("X", conv.OUTOPTIONS, "noiso")
        #conv.AddOption("X", conv.OUTOPTIONS, "nochg")
        #conv.AddOption("X", conv.OUTOPTIONS, "nostereo")
        conv.AddOption("w", conv.OUTOPTIONS)
        obmol = openbabel.OBMol()
        s = urllib.urlopen('http://rest.kegg.jp/get/cpd:C%05d/mol' % cid).read()
        conv.ReadString(obmol, s)
        inchi = conv.WriteString(obmol).strip()
        return inchi
        
    @staticmethod
    def smiles2inchi(smiles):
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('smiles', 'inchi')
        conv.AddOption("X", conv.OUTOPTIONS, "FixedH")
        conv.AddOption("X", conv.OUTOPTIONS, "RecMet")
        conv.AddOption("w", conv.OUTOPTIONS)
        obmol = openbabel.OBMol()
        conv.ReadString(obmol, smiles)
        inchi = conv.WriteString(obmol).strip()
        return inchi

    @staticmethod
    def get_formula_from_inchi(inchi):
        tokens = re.findall('/f([0-9A-Za-z\.]+/)', inchi)
        if len(tokens) == 0:
            tokens = re.findall('InChI=1S?/([0-9A-Za-z\.]+)', inchi)

        if len(tokens) == 1:
            return tokens[0]
        elif len(tokens) > 1:
            raise ValueError('Bad InChI: ' + inchi)
        else:
            return ''
                
    @staticmethod
    def get_atom_bag_and_charge_from_inchi(inchi):
        fixed_charge = 0
        for q in re.findall('/q([0-9\+\-\;]+)', inchi):
            for s in q.split(';'): 
                if s:
                    fixed_charge += int(s)

        fixed_protons = 0
        for p in re.findall('/p([0-9\+\-\;]+)', inchi):
            for s in p.split(';'):
                if s:
                    fixed_protons += int(s)
        
        formula = Compound.get_formula_from_inchi(inchi)

        atom_bag = {}
        for mol_formula_times in formula.split('.'):
            for times, mol_formula in re.findall('^(\d+)?(\w+)', mol_formula_times):
                if not times:
                    times = 1
                else:
                    times = int(times)
                for atom, count in re.findall("([A-Z][a-z]*)([0-9]*)", mol_formula):
                    if count == '':
                        count = 1
                    else:
                        count = int(count)
                    atom_bag[atom] = atom_bag.get(atom, 0) + count * times
        
        if fixed_protons:
            atom_bag['H'] = atom_bag.get('H', 0) + fixed_protons
            fixed_charge += fixed_protons
        return atom_bag, fixed_charge
    
    @staticmethod
    def get_species_pka(inchi):
        try:
            pka_list, major_ms = GetDissociationConstants(inchi)
            pka_list = sorted([pka for pka in pka_list if pka > MIN_PH and pka < MAX_PH], reverse=True)
            major_ms_inchi = Compound.smiles2inchi(major_ms)
        except ChemAxonError:
            logging.warning('chemaxon failed to find pKas for this inchi: ' + inchi)
            pka_list = []
            major_ms_inchi = inchi

        atom_bag, major_ms_charge = Compound.get_atom_bag_and_charge_from_inchi(major_ms_inchi)
        major_ms_nH = atom_bag.get('H', 0)

        n_species = len(pka_list) + 1
        if pka_list == []:
            major_ms_index = 0
        else:
            major_ms_index = len([1 for pka in pka_list if pka > 7])
            
        pKas = np.zeros((n_species, n_species), dtype=float)
        for i, pka in enumerate(pka_list):
            pKas[i+1, i] = pka
            pKas[i, i+1] = pka

        majorMSpH7 = np.zeros((n_species, 1), dtype=int)
        nHs = np.zeros((n_species, 1), dtype=int)
        zs = np.zeros((n_species, 1), dtype=int)

        majorMSpH7[major_ms_index, 0] = 1
        for i in xrange(n_species):
            zs[i, 0] = (i - major_ms_index) + major_ms_charge
            nHs[i, 0] = (i - major_ms_index) + major_ms_nH
        
        return pKas, majorMSpH7, nHs, zs
    
    def __str__(self):
        return "CID: C%05d\nInChI: %s\npKas:%s\nmajor MS:%s\nnHs: %s\nzs: %s" % \
            (self.cid, self.inchi, str(self.pKas), str(self.majorMSpH7.T), str(self.nHs.T), str(self.zs.T))
    
    def get_atom_bag_with_electrons(self):
        """
            Calculates the number of electrons in a given molecule
            Returns:
                a dictionary of all element counts and also electron count ('e-')
        """
        if self.inchi == '':
            return None
        atom_bag, charge = Compound.get_atom_bag_and_charge_from_inchi(self.inchi)
        n_protons = sum([count * Compound._obElements.GetAtomicNum(str(elem))
                         for (elem, count) in atom_bag.iteritems()])
        atom_bag['e-'] = n_protons - charge
        return atom_bag

if __name__ == '__main__':
    print get_inchi(10)
