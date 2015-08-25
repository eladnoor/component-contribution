import openbabel
import types
import re
import chemaxon
from thermodynamic_constants import default_T, default_pH

class OpenBabelError(Exception):
    pass

class Molecule(object):

    # for more rendering options visit:
    # http://www.ggasoftware.com/opensource/indigo/api/options#rendering
    _obElements = openbabel.OBElementTable()
    _obSmarts = openbabel.OBSmartsPattern()
    
    @staticmethod
    def GetNumberOfElements():
        return Molecule._obElements.GetNumberOfElements()
    
    @staticmethod
    def GetAllElements():
        return [Molecule._obElements.GetSymbol(i) for i in 
                xrange(Molecule.GetNumberOfElements())]

    @staticmethod
    def GetSymbol(atomic_num):
        return Molecule._obElements.GetSymbol(atomic_num)
            
    @staticmethod
    def GetAtomicNum(elem):
        if type(elem) == types.UnicodeType:
            elem = str(elem)
        return Molecule._obElements.GetAtomicNum(elem)
    
    @staticmethod
    def VerifySmarts(smarts):
        return Molecule._obSmarts.Init(smarts)
    
    def __init__(self):
        self.title = None
        self.obmol = openbabel.OBMol()
        self.smiles = None
        self.inchi = None

    def __str__(self):
        return self.title or self.smiles or self.inchi or ""
        
    def __len__(self):
        return self.GetNumAtoms()
    
    def Clone(self):
        tmp = Molecule()
        tmp.title = self.title
        tmp.obmol = openbabel.OBMol(self.obmol)
        tmp.smiles = self.smiles
        tmp.inchi = self.inchi
        return tmp
    
    def SetTitle(self, title):
        self.title = title 
    
    @staticmethod
    def FromSmiles(smiles):
        m = Molecule()
        m.smiles = smiles
        obConversion = openbabel.OBConversion()
        obConversion.AddOption("w", obConversion.OUTOPTIONS)
        obConversion.SetInFormat("smiles")
        if not obConversion.ReadString(m.obmol, m.smiles):
            raise OpenBabelError("Cannot read the SMILES string: " + smiles)
        try:
            m.UpdateInChI()
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from SMILES: " + smiles)
        m.SetTitle(smiles)
        return m
        
    @staticmethod
    def FromInChI(inchi):
        m = Molecule()
        m.inchi = inchi
        obConversion = openbabel.OBConversion()
        obConversion.AddOption("w", obConversion.OUTOPTIONS)
        obConversion.SetInFormat("inchi")
        obConversion.ReadString(m.obmol, m.inchi)
        try:
            m.UpdateSmiles()
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from InChI: " + inchi)
        m.SetTitle(inchi)
        return m
    
    @staticmethod
    def FromMol(mol):
        m = Molecule()
        obConversion = openbabel.OBConversion()
        obConversion.AddOption("w", obConversion.OUTOPTIONS)
        obConversion.SetInFormat("mol")
        obConversion.ReadString(m.obmol, mol)
        try:
            m.UpdateInChI()
            m.UpdateSmiles()
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from MOL file:\n" + mol)
        m.SetTitle("")
        return m

    @staticmethod
    def FromOBMol(obmol):
        m = Molecule()
        m.obmol = obmol
        try:
            m.UpdateInChI()
            m.UpdateSmiles()
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from OBMol")
        m.SetTitle("")
        return m
    
    @staticmethod
    def _FromFormat(s, fmt='inchi'):
        if fmt == 'smiles' or fmt == 'smi':
            return Molecule.FromSmiles(s)
        if fmt == 'inchi':
            return Molecule.FromInChI(s)
        if fmt == 'mol':
            return Molecule.FromMol(s)
        if fmt == 'obmol':
            return Molecule.FromOBMol(s)
    
    @staticmethod
    def _ToFormat(obmol, fmt='inchi'):
        obConversion = openbabel.OBConversion()
        obConversion.AddOption("w", obConversion.OUTOPTIONS)
        obConversion.SetOutFormat(fmt)
        res = obConversion.WriteString(obmol)
        if not res:
            raise OpenBabelError("Cannot convert OBMol to %s" % fmt)
        if fmt == 'smiles' or fmt == 'smi':
            res = res.split()
            if res == []:
                raise OpenBabelError("Cannot convert OBMol to %s" % fmt)
            else:
                return res[0]
        elif fmt == 'inchi':
            return res.strip()
        else:
            return res
        
    @staticmethod
    def Smiles2InChI(smiles):
        obConversion = openbabel.OBConversion()
        obConversion.AddOption("w", obConversion.OUTOPTIONS)
        obConversion.SetInAndOutFormats("smiles", "inchi")
        obmol = openbabel.OBMol()
        if not obConversion.ReadString(obmol, smiles):
            raise OpenBabelError("Cannot read the SMILES string: " + smiles)
        return obConversion.WriteString(obmol).strip()

    @staticmethod
    def InChI2Smiles(inchi):
        obConversion = openbabel.OBConversion()
        obConversion.AddOption("w", obConversion.OUTOPTIONS)
        obConversion.SetInAndOutFormats("inchi", "smiles")
        obmol = openbabel.OBMol()
        if not obConversion.ReadString(obmol, inchi):
            raise OpenBabelError("Cannot read the InChI string: " + inchi)
        return obConversion.WriteString(obmol).split()[0]
        
    def RemoveHydrogens(self):
        self.obmol.DeleteHydrogens()
    
    def RemoveAtoms(self, indices):
        self.obmol.BeginModify()
        for i in sorted(indices, reverse=True):
            self.obmol.DeleteAtom(self.obmol.GetAtom(i+1))
        self.obmol.EndModify()
        self.smiles = None
        self.inchi = None
        
    def SetAtomicNum(self, index, new_atomic_num):
        self.obmol.GetAtom(index+1).SetAtomicNum(new_atomic_num)
        self.smiles = None
        self.inchi = None
        
    def ToOBMol(self):
        return self.obmol
    
    def ToFormat(self, fmt='inchi'):
        return Molecule._ToFormat(self.obmol, fmt=fmt)
    
    def ToMolfile(self):
        return self.ToFormat('mol')

    def UpdateInChI(self):
        self.inchi = Molecule._ToFormat(self.obmol, 'inchi')

    def ToInChI(self):
        """ 
            Lazy storage of the InChI identifier (calculate once only when 
            asked for and store for later use).
        """
        if not self.inchi:
            self.UpdateInChI()
        return self.inchi
    
    def UpdateSmiles(self):
        self.smiles = Molecule._ToFormat(self.obmol, 'smiles')
    
    def ToSmiles(self):
        """ 
            Lazy storage of the SMILES identifier (calculate once only when 
            asked for and store for later use).
        """
        if not self.smiles:
            self.UpdateSmiles()
        return self.smiles
    
    def GetFormula(self):
        tokens = re.findall('InChI=1S?/([0-9A-Za-z\.]+)', self.ToInChI())
        if len(tokens) == 1:
            return tokens[0]
        elif len(tokens) > 1:
            raise ValueError('Bad InChI: ' + self.ToInChI())
        else:
            return ''
    
    def GetExactMass(self):
        return self.obmol.GetExactMass()
    
    def GetAtomBagAndCharge(self):
        inchi = self.ToInChI()
        atom_bag, major_ms_charge = chemaxon.GetAtomBagAndCharge(inchi)
        return atom_bag, major_ms_charge

    def GetHydrogensAndCharge(self):
        atom_bag, charge = self.GetAtomBagAndCharge()
        return atom_bag.get('H', 0), charge
        
    def GetNumElectrons(self):
        """Calculates the number of electrons in a given molecule."""
        atom_bag, fixed_charge = self.GetAtomBagAndCharge()
        return atom_bag.get('e-', 0)
    
    def GetNumAtoms(self):
        return self.obmol.NumAtoms()

    def GetAtoms(self):
        return [self.obmol.GetAtom(i+1) for i in xrange(self.obmol.NumAtoms())]
    
    def FindSmarts(self, smarts):
        """
        Corrects the pyBel version of Smarts.findall() which returns results as tuples,
        with 1-based indices even though Molecule.atoms is 0-based.
    
        Args:
            mol: the molecule to search in.
            smarts_str: the SMARTS query to search for.
        
        Returns:
            The re-mapped list of SMARTS matches.
        """
        Molecule._obSmarts.Init(smarts)
        if Molecule._obSmarts.Match(self.obmol):
            match_list = Molecule._obSmarts.GetMapList()
            shift_left = lambda m: [(n - 1) for n in m] 
            return map(shift_left, match_list)
        else:
            return []

    def GetAtomCharges(self):
        """
            Returns:
                A list of charges, according to the number of atoms
                in the molecule
        """
        return [atom.GetFormalCharge() for atom in self.GetAtoms()]

if __name__ == '__main__':
    mol = Molecule.FromInChI('InChI=1S/H2/h1H')
    print mol.ExactMass()
