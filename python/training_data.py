import os, logging, csv
import numpy as np
from thermodynamic_constants import R, F
from compound_cacher import CompoundCacher
import kegg_reaction

class TrainingData(object):
    
    # a dictionary of the filenames of the training data and the relative 
    # weight of each one    
    FNAME_DICT = {'TECRDB' : ('../data/TECRDB.tsv', 1.0),
                  'FORMATION' : ('../data/formation_energies_transformed.tsv', 1.0),
                  'REDOX' : ('../data/redox.tsv', 1.0)}

    def __del__(self):
        self.ccache.dump()

    def __init__(self):
        self.ccache = CompoundCacher.getInstance()
    
        # verify that the files exist
        for fname, _ in TrainingData.FNAME_DICT.values():
            if not os.path.exists(fname):
                raise Exception('file not found: ' + fname)
        
        tecrdb_params = TrainingData.read_tecrdb()
        
        formation_params, cids_that_dont_decompose = TrainingData.read_formations()
        
        redox_params = TrainingData.read_redox()
        
        thermo_params = tecrdb_params + formation_params + redox_params
        
        cids = set()
        for d in thermo_params:
            cids = cids.union(d['reaction'].keys())
        cids = sorted(cids)
        
        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids'.
        self.S = np.zeros((len(cids), len(thermo_params)))
        for k, d in enumerate(thermo_params):
            for cid, coeff in d['reaction'].iteritems():
                self.S[cids.index(cid), k] = coeff
            
        self.cids = cids;
        self.cids_that_dont_decompose = cids_that_dont_decompose

        # get the InChIs and pKas for all the compounds in the training data
        # (note that all of them have KEGG IDs)
        self.cid2compound = {}
        for cid in self.cids:
            self.cid2compound[cid] = self.ccache.get_kegg_compound(cid)
        
        self.dG0_prime = np.array([d['dG\'0'] for d in thermo_params])
        self.T = np.array([d['T'] for d in thermo_params])
        self.I = np.array([d['I'] for d in thermo_params])
        self.pH = np.array([d['pH'] for d in thermo_params])
        self.pMg = np.array([d['pMg'] for d in thermo_params])
        self.weight = np.array([d['weight'] for d in thermo_params])
        rxn_inds_to_balance = [i for i in xrange(len(thermo_params))
                               if thermo_params[i]['balance']]

        self.balance_reactions(rxn_inds_to_balance)

    @staticmethod
    def read_tecrdb():       
        """Read the raw data of TECRDB (NIST)"""
        fname, weight = TrainingData.FNAME_DICT['TECRDB']

        thermo_params = [] # columns are: reaction, dG'0, T, I, pH, pMg, weight, balance?

        headers = ["URL", "REF_ID", "METHOD", "EVAL", "EC", "ENZYME NAME",
                   "REACTION IN KEGG IDS", "REACTION IN COMPOUND NAMES",
                   "K", "K'", "T", "I", "pH", "pMg"]

        for row_list in csv.reader(open(fname, 'r'), delimiter='\t'):
            if row_list == []:
                continue
            row = dict(zip(headers, row_list))
            if (row['K\''] == '') or (row['T'] == '') or (row['pH'] == ''):
                continue
            
            # parse the reaction
            sparse = kegg_reaction.parse_kegg_reaction(row['REACTION IN KEGG IDS'], arrow='=')

            # calculate dG'0
            dG0_prime = -R * float(row['T']) * np.log(float(row['K\''])) 
            try:
                thermo_params.append({'reaction': sparse,
                                      'dG\'0' : dG0_prime,
                                      'T': float(row['T']), 
                                      'I': float(row['I'] or 0),
                                      'pH': float(row['pH']),
                                      'pMg': float(row['pMg'] or 0),
                                      'weight': weight,
                                      'balance': True})
            except ValueError:
                raise Exception('Cannot parse row: ' + str(row))

        logging.info('Successfully added %d reactions from TECRDB' % len(thermo_params))
        return thermo_params
        
    @staticmethod
    def read_formations():
        """Read the Formation Energy data"""
        fname, weight = TrainingData.FNAME_DICT['FORMATION']
        # columns are: reaction, dG'0, T, I, pH, pMg, weight, balance?
        thermo_params = []
        cids_that_dont_decompose = set()
        
        # fields are: cid, name, dG'0, pH, I, pMg, T, decompose?,
        #             compound_ref, remark
        for row in csv.DictReader(open(fname, 'r'), delimiter='\t'):
            if row['dG\'0'] == '':
                continue
            cid = int(row['cid'])
            if int(row['decompose']) == 0:
                cids_that_dont_decompose.add(cid)

            thermo_params.append({'reaction': {cid : 1},
                                  'dG\'0' : float(row['dG\'0']),
                                  'T': float(row['T']), 
                                  'I': float(row['I'] or 0),
                                  'pH': float(row['pH']),
                                  'pMg': float(row['pMg'] or 0),
                                  'weight': weight,
                                  'balance': False})

        logging.info('Successfully added %d formation energies' % len(thermo_params))
        return thermo_params, cids_that_dont_decompose
        
    @staticmethod
    def read_redox():
        """Read the Reduction potential data"""
        
        fname, weight = TrainingData.FNAME_DICT['REDOX']
        # columns are: reaction, dG'0, T, I, pH, pMg, weight, balance?
        thermo_params = []
        
        # fields are: name, CID_ox, nH_ox, charge_ox, CID_red,
        #             nH_red, charge_red, E'0, pH, I, pMg, T, ref
        for row in csv.DictReader(open(fname, 'r'), delimiter='\t'):
            cid_ox = int(row['CID_ox'])
            cid_red = int(row['CID_red'])
            delta_nH = float(row['nH_red']) - float(row['nH_ox'])
            delta_charge = float(row['charge_red']) - float(row['charge_ox'])
            delta_e = delta_nH - delta_charge
            dG0_prime = -F * float(row['E\'0']) * delta_e
            
            thermo_params.append({'reaction': {cid_ox : -1, cid_red : 1},
                                  'dG\'0' : dG0_prime,
                                  'T': float(row['T']), 
                                  'I': float(row['I'] or 0),
                                  'pH': float(row['pH']),
                                  'pMg': float(row['pMg'] or 0),
                                  'weight': weight,
                                  'balance': False})        

        logging.info('Successfully added %d redox potentials' % len(thermo_params))
        return thermo_params
    
    def balance_reactions(self, rxn_inds_to_balance):
        """
            use the chemical formulas from the InChIs to verify that each and every
            reaction is balanced
        """
        elements = set()
        atom_bag_list = []
        for cid in self.cids:
            atom_bag = self.cid2compound[cid].get_atom_bag_with_electrons()
            if atom_bag is not None:
                elements = elements.union(atom_bag.keys())
            atom_bag_list.append(atom_bag)
        elements.discard('H') # no need to balance H atoms (balancing electrons is sufficient)
        elements = sorted(elements)
        
        Ematrix = np.zeros((len(atom_bag_list), len(elements)))
        cpd_inds_without_formula = [] # we will have to skip reactions that contain them
        for i, atom_bag in enumerate(atom_bag_list):
            if atom_bag is None:
                cpd_inds_without_formula.append(i)
                continue
            for j, elem in enumerate(elements):
                Ematrix[i, j] = atom_bag.get(elem, 0)

        rxn_inds_without_formula = np.nonzero(np.sum(np.abs(self.S[cpd_inds_without_formula, :]), 0))[0]
        rxn_inds_to_balance = set(rxn_inds_to_balance).difference(rxn_inds_without_formula)

        # need to check that all elements are balanced (except H, but including e-)
        # if only O is not balanced, add water molecules
        if 'O' in elements:
            i_H2O = self.cids.index(1)
            j_O = elements.index('O')
            conserved = np.dot(Ematrix.T, self.S)
            for k in rxn_inds_to_balance:
                self.S[i_H2O, k] = self.S[i_H2O, k] - conserved[j_O, k]

        # recalculate conservation matrix
        conserved = np.dot(Ematrix.T, self.S)
        
        rxn_inds_to_remove = [k for k in rxn_inds_to_balance 
                              if not np.all(conserved[:, k] == 0)]
        rxn_inds_to_keep = \
            set(range(self.S.shape[1])).difference(rxn_inds_to_remove)
        
        rxn_inds_to_keep = sorted(rxn_inds_to_keep)
        
        self.S = self.S[:, rxn_inds_to_keep]
        self.dG0_prime = self.dG0_prime[:, rxn_inds_to_keep]
        self.T = self.T[:, rxn_inds_to_keep]
        self.I = self.I[:, rxn_inds_to_keep]
        self.pH = self.pH[:, rxn_inds_to_keep]
        self.pMg = self.pMg[:, rxn_inds_to_keep]
        self.weight = self.weight[:, rxn_inds_to_keep]

        logging.info('After removing %d unbalanced reactions, the stoichiometric '
                     'matrix contains: '
                     '%d compounds and %d reactions' %
                     (len(rxn_inds_to_remove), self.S.shape[0], self.S.shape[1]))


if __name__ == '__main__':
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)

    td = TrainingData()
