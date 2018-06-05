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

import logging
import platform

import pandas

from six import StringIO
from subprocess import Popen, PIPE

from component_contribution import exceptions

if platform.system() == 'Windows':
    CXCALC_BIN = 'C:\\Program Files (x86)\\ChemAxon\\MarvinBeans\\bin\\cxcalc.bat'
    use_shell_for_echo = True
else:
    CXCALC_BIN = 'cxcalc'
    use_shell_for_echo = False


LOGGER = logging.getLogger(__name__)

MID_PH = 7.0
N_PKAS = 20


def run_cxcalc(molstring, args):
    """
    Runs cxcalc.

    Parameters
    ----------
    molstring : str
        A text description of the molecule (SMILES or InChI).
    args : list
        A list of arguments for cxcalc.

    Returns
    -------
    output : str
        The cxcalc output.

    Raises
    ------
    ChemAxonRuntimeError
        If the command fails.

    """

    with open(platform.DEV_NULL, 'w') as dev_null:
        try:
            LOGGER.debug("INPUT: echo %s | %s" % (molstring, ' '.join([CXCALC_BIN] + args)))
            p1 = Popen(["echo", molstring], stdout=PIPE, shell=use_shell_for_echo)
            p2 = Popen([CXCALC_BIN] + args, stdin=p1.stdout,
                       executable=CXCALC_BIN, stdout=PIPE, stderr=dev_null, shell=False)

            res = p2.communicate()[0].decode('utf-8')
            if p2.returncode != 0:
                raise exceptions.ChemAxonRuntimeError(str(args))

            LOGGER.debug("OUTPUT: %s" % res)
            return res
        except OSError:
            raise Exception("Marvin (by ChemAxon) must be installed to calculate pKa data.")


def parse_pka_output(output, n_acidic, n_basic):
    """
    Parses the pKa command output.

    Parameters
    ----------
    output : string
        cxcalc output string.
    n_acidic : int
        The max number of acidic pKas to calculate.
    n_basic : int
        The max number of basic pKas to calculate.

    Returns
    -------
    atom2pka : dict
        A dictionary that maps the atom index to a list of pKas that are assigned to that atom.
    """
    atom2pKa = {}

    try:
        pkaline = output.split('\n')[1]
    except IndexError:
        raise ValueError('Cannot parse pKa output: ' + output)
    splitline = pkaline.split('\t')
    splitline.pop(0)

    if n_acidic + n_basic > 0:
        if len(splitline) != (n_acidic + n_basic + 2):
            raise exceptions.ChemAxonRuntimeError('ChemAxon failed to find any pKas')

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

        if atom_list:  # a comma separated list of the deprotonated atoms
            atom_numbers = [int(y)-1 for y in atom_list.split(',')]
            for i, j in enumerate(atom_numbers):
                atom2pKa.setdefault(j, [])
                atom2pKa[j].append((pKa_list[i], acid_or_base_list[i]))

    smiles_list = splitline
    return atom2pKa, smiles_list


def get_dissociation_constants(molstring, n_acidic=N_PKAS, n_basic=N_PKAS, p_h=MID_PH):
    """

    Computes the dissociation constants and major microspecies at a defined pH.

    Parameters
    ----------
    molstring : str
        A text description of the molecule (SMILES or InChI).
    n_acidic : int
        The max number of acidic pKas to calculate.
    n_basic : int
        The max number of basic pKas to calculate.
    p_h : float
        The pH for which the major pseudoisomer is calculated.

    Returns
    -------
    all_p_kas : list
        A list of floats (pKa values).
    major_ms : string
        SMILES string of the major pseudoisomer at pH.
    """
    args = []
    if n_acidic + n_basic > 0:
        args += ['pka', '-a', str(n_acidic), '-b', str(n_basic), 'majorms', '-M', 'true', '--pH', str(p_h)]

    output = str(run_cxcalc(molstring, args))
    atom2pka, smiles_list = parse_pka_output(output, n_acidic, n_basic)

    all_p_kas = []
    for p_ka_list in atom2pka.values():
        all_p_kas += [pka for pka, _ in p_ka_list]

    return sorted(all_p_kas), smiles_list[0]


def get_formula_and_charge(molstring):
    """
    Calculates the formula and charge of a molecule.

    Parameters
    ----------
    molstring : str
        A text description of the molecule (SMILES or InChI).

    Returns
    -------
    formula : str
        The chemical formula of the molecule.
    formal_charge : int
        The formal charge of the molecule.
    """
    args = ['formula', 'formalcharge']
    output = run_cxcalc(molstring, args)
    # the output is a tab separated table whose columns are:
    # id, Formula, Formal charge
    output_df = pandas.read_csv(StringIO(output), delimiter='\t')
    if output_df.columns.tolist() != ['id', 'Formula', 'Formal charge']:
        raise exceptions.ChemAxonRuntimeError('cannot get the formula and charge for: ' + molstring)

    formula = output_df['Formula'].iat[0]
    try:
        formal_charge = int(output_df['Formal charge'].iat[0])
    except ValueError:
        formal_charge = 0

    return formula, formal_charge


try:
    run_cxcalc("CCCC", ["formalcharge"])
except Exception as exception:
    new_exception = exceptions.ChemAxonNotAvailable()

    raise new_exception


if __name__ == "__main__":
    LOGGER.setLevel(logging.WARNING)
    from component_contribution import Molecule
    compound_list = [('orthophosphate', 'InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-3'),
                     ('D-Erythrulose', 'InChI=1S/C4H8O4/c5-1-3(7)4(8)2-6/h3,5-7H,1-2H2/t3-/m1/s1')]
    for name, inchi in compound_list:
        print("Formula: %s\nCharge: %d" % get_formula_and_charge(inchi))
        diss_table, major_ms = get_dissociation_constants(inchi)
        m = Molecule.FromSmiles(major_ms)
        print("Name: %s\nInChI: %s\npKas: %s" % (name, m.ToInChI(), str(diss_table)))
