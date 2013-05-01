import logging, itertools, os
import numpy as np
from subprocess import Popen, PIPE

CXCALC_BIN = "cxcalc"
MID_PH = 7.0
N_PKAS = 20

class ChemAxonError(Exception):
    pass

def RunCxcalc(molstring, args):
    devnull = open('/dev/null', 'w')
    try:
        p1 = Popen(["echo", molstring], stdout=PIPE)
        p2 = Popen([CXCALC_BIN] + args, stdin=p1.stdout,
                   executable=CXCALC_BIN, stdout=PIPE, stderr=devnull)
        logging.debug("INPUT: echo %s | %s" % (molstring, ' '.join([CXCALC_BIN] + args)))
        #p.wait()
        #os.remove(temp_fname)
        res = p2.communicate()[0]
        if p2.returncode != 0:
            raise ChemAxonError(debug_args)
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
            atom_numbers = [int(x)-1 for x in atom_list.split(',')]
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

if __name__ == "__main__":
    
    from molecule import Molecule
    compound_list = [('glycine', 'C(=O)(O)CN'),
                     ('CO2', 'O=C=O'),
                     ('ATP', 'Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])(=O)OP(O)([O-])=O)C(O)C1O'),
                     ('3-Ketoarabinitol', 'OCC(O)C(C(O)CO)=O')]
    
    for name, smiles in compound_list:
        diss_table, major_ms = GetDissociationConstants(smiles)
        m = Molecule.FromSmiles(major_ms)
        print name, m.ToInChI(), str(diss_table)
