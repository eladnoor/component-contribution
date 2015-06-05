import logging, csv, re, platform
import StringIO
from subprocess import Popen, PIPE
import openbabel

if platform.system() == 'Windows':
    CXCALC_BIN = 'C:\\Program Files (x86)\\ChemAxon\\MarvinBeans\\bin\\cxcalc.bat'
    use_shell_for_echo = True
else:
    CXCALC_BIN = 'cxcalc'
    use_shell_for_echo = False

MID_PH = 7.0
N_PKAS = 20
_obElements = openbabel.OBElementTable()

class ChemAxonError(Exception):
    pass

def RunCxcalc(molstring, args):
    with open(platform.DEV_NULL, 'w') as dev_null:
        try:
            logging.debug("INPUT: echo %s | %s" % (molstring, ' '.join([CXCALC_BIN] + args)))
            p1 = Popen(["echo", molstring], stdout=PIPE, shell=use_shell_for_echo)
            p2 = Popen([CXCALC_BIN] + args, stdin=p1.stdout,
                       executable=CXCALC_BIN, stdout=PIPE, stderr=dev_null, shell=False)
            #p.wait()
            #os.remove(temp_fname)
            res = p2.communicate()[0]
            if p2.returncode != 0:
                raise ChemAxonError(str(args))
            logging.debug("OUTPUT: %s" % res)
            return res
        except OSError:
            raise Exception("Marvin (by ChemAxon) must be installed to calculate pKa data.")

def ParsePkaOutput(s, n_acidic, n_basic):
    """
        Returns:
            A dictionary that maps the atom index to a list of pKas
            that are assigned to that atom.
    """
    atom2pKa = {}

    pkaline = s.split('\n')[1]
    splitline = pkaline.split('\t')
    splitline.pop(0)
    
    if n_acidic + n_basic > 0:
        if len(splitline) != (n_acidic + n_basic + 2):
            raise ChemAxonError('ChemAxon failed to find any pKas')
        
        pKa_list = []
        acid_or_base_list = []
        for i in range(n_acidic + n_basic):
            x = splitline.pop(0) 
            if x == '':
                continue
            
            pKa_list.append(float(x))
            if i < n_acidic:
                acid_or_base_list.append('acid')
            else:
                acid_or_base_list.append('base')
        
        atom_list = splitline.pop(0)

        if atom_list: # a comma separated list of the deprotonated atoms
            atom_numbers = [int(y)-1 for y in atom_list.split(',')]
            for i, j in enumerate(atom_numbers):
                atom2pKa.setdefault(j, [])
                atom2pKa[j].append((pKa_list[i], acid_or_base_list[i]))
    
    smiles_list = splitline
    return atom2pKa, smiles_list

def _GetDissociationConstants(molstring, n_acidic=N_PKAS, n_basic=N_PKAS,
                              pH=MID_PH):
    """
        Returns:
            A pair of (pKa list, major pseudoisomer)
            
            - the pKa list is of the pKa values in ascending order.
            - the major pseudoisomer is a SMILES string of the major species
              at the given pH.
    """
    args = []
    if n_acidic + n_basic > 0:
        args += ['pka', '-a', str(n_acidic), '-b', str(n_basic),
                 'majorms', '-M', 'true', '--pH', str(pH)]
    
    output = RunCxcalc(molstring, args)
    atom2pKa, smiles_list = ParsePkaOutput(output, n_acidic, n_basic)
    
    all_pKas = []
    for pKa_list in atom2pKa.values():
        all_pKas += [pKa for pKa, _ in pKa_list]
    
    return sorted(all_pKas), smiles_list

def GetDissociationConstants(molstring, n_acidic=N_PKAS, n_basic=N_PKAS,
                             pH=MID_PH):
    """
        Arguments:
            molstring - a text description of the molecule (SMILES or InChI)
            n_acidic  - the max no. of acidic pKas to calculate
            n_basic   - the max no. of basic pKas to calculate
            pH        - the pH for which the major pseudoisomer is calculated

        Returns a pair:
            (all_pKas, major_ms)
            
        - all_pKas is a list of floats (pKa values)
        - major_ms is a SMILES string of the major pseudoisomer at pH_mid 
    """
    all_pKas, smiles_list = _GetDissociationConstants(molstring, n_acidic, 
                                                      n_basic, pH)
    major_ms = smiles_list[0]
    return all_pKas, major_ms

def GetFormulaAndCharge(molstring):
    """
        Arguments:
            molstring - a text description of the molecule (SMILES or InChI)

        Returns:
            chemical formula of the molecule
    """
    args = ['formula', 'formalcharge']
    output = RunCxcalc(molstring, args)
    # the output is a tab separated table whose columns are:
    # id, Formula, Formal charge
    f = StringIO.StringIO(output)
    tsv_output = csv.reader(f, delimiter='\t')
    headers = tsv_output.next()
    if headers != ['id', 'Formula', 'Formal charge']:
        raise ChemAxonError('cannot get the formula and charge for: ' + molstring)
    _, formula, formal_charge = tsv_output.next()

    try:
        formal_charge = int(formal_charge)
    except ValueError:
        formal_charge = 0
    
    return formula, formal_charge

def GetAtomBagAndCharge(molstring):
    formula, formal_charge = GetFormulaAndCharge(molstring)

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
    
    n_protons = sum([c * _obElements.GetAtomicNum(str(elem))
                     for (elem, c) in atom_bag.iteritems()])
    atom_bag['e-'] = n_protons - formal_charge

    return atom_bag, formal_charge

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.WARNING)
    from molecule import Molecule
    compound_list = [('D-Erythrulose', 'InChI=1S/C4H8O4/c5-1-3(7)4(8)2-6/h3,5-7H,1-2H2/t3-/m1/s1')]
    
    for name, inchi in compound_list:
        print "Formula: %s\nCharge: %d" % GetFormulaAndCharge(inchi)
        diss_table, major_ms = GetDissociationConstants(inchi)
        m = Molecule.FromSmiles(major_ms)
        print "Name: %s\nInChI: %s\npKas: %s" % (name, m.ToInChI(), str(diss_table))
