import os, logging, csv
import numpy as np
import pandas as pd
from scipy.io import savemat
from .thermodynamic_constants import R, F
from .compound_cache import CompoundCache
from .kegg_reaction import KeggReaction
from pkg_resources import resource_stream

class TrainingData(object):
    
    # a dictionary of the filenames of the training data and the relative 
    # weight of each one
    FNAME_DICT = {'TECRDB' : ('data/TECRDB.csv', 1.0),
                  'FORMATION' : ('data/formation_energies_transformed.csv', 1.0),
                  'REDOX' : ('data/redox.csv', 1.0)}

    def __del__(self):
        self.ccache.dump()

    def __init__(self):
        self.ccache = CompoundCache()
        
        thermo_df, self.cids_that_dont_decompose = \
            TrainingData.get_all_thermo_params()
        
        cids = set()
        for rxn in thermo_df['reaction']:
            cids = cids.union(rxn.keys())
        cids = sorted(cids)
        
        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids'.
        self.S = np.zeros((len(cids), thermo_df.shape[0]))
        for k, rxn in enumerate(thermo_df['reaction'].values):
            for cid, coeff in rxn.items():
                self.S[cids.index(cid), k] = coeff
            
        self.cids = cids

        self.dG0_prime = np.array(thermo_df["dG'0"])
        self.T = np.array(thermo_df["T"])
        self.I = np.array(thermo_df["I"])
        self.pH = np.array(thermo_df["pH"])
        self.pMg = np.array(thermo_df["pMg"])
        self.weight = np.array(thermo_df["weight"])
        self.reference = thermo_df["reference"].tolist()
        self.description = thermo_df["description"].tolist()
        
        balance = thermo_df["balance"].tolist()
        rxn_inds_to_balance = \
            [i for i in range(thermo_df.shape[0]) if balance[i]]

        self.balance_reactions(rxn_inds_to_balance)
        
        self.reverse_transform()

    def savemat(self, file_name):
        """
            Write all training data to a Matlab file.
            
            Arguments:
                file_name - str or file-like object to which data will be written
        """
        d = {'dG0_prime': self.dG0_prime,
             'dG0': self.dG0,
             'T': self.T,
             'I': self.I,
             'pH': self.pH,
             'pMg': self.pMg,
             'weight': self.weight,
             'cids': self.cids,
             'S': self.S}
        savemat(file_name, d, oned_as='row')

    def savecsv(self, fname):
        csv_output = csv.writer(open(fname, 'w'))
        csv_output.writerow(['reaction', 'T', 'I', 'pH', 'reference', 'dG0', 'dG0_prime'])
        for j in range(self.S.shape[1]):
            sparse = {self.cids[i]: self.S[i, j] for i in range(self.S.shape[0])}
            r_string = KeggReaction(sparse).write_formula()
            csv_output.writerow([r_string, self.T[j], self.I[j], self.pH[j],
                                 self.reference[j], self.dG0[j], self.dG0_prime[j]])

    @staticmethod
    def str2double(s):
        """
            casts a string to float, but if the string is empty return NaN
        """
        if s == '':
            return np.nan
        else:
            return float(s)

    @staticmethod
    def read_tecrdb(resource_name, weight):
        """
            Read the raw data of TECRDB (NIST)
        """
        
        tecr_df = pd.read_csv(resource_stream('component_contribution',
                                              '/data/TECRDB.csv'))

        for col in ["T", "I", "pH", "pMg", "K'"]:
            tecr_df[col] = tecr_df[col].apply(float)
        
        # remove rows with missing crucial data
        tecr_df = tecr_df[~pd.isnull(tecr_df[["K'", "T", "pH"]]).any(axis=1)]
        
        # calculate the dG'0 from the Keq and T
        tecr_df["dG'0"] = -R * tecr_df["T"] * np.log(tecr_df["K'"])
        
        # parse the reaction
        tecr_df['reaction'] = tecr_df['REACTION IN KEGG IDS'].apply(lambda r:
            KeggReaction.parse_formula(r, arrow='='))

        tecr_df.rename(columns={'REF_ID': 'reference',
                                'REACTION IN COMPOUND NAMES': 'description'},
                       inplace=True)
        tecr_df['weight'] = weight
        tecr_df['balance'] = True
        tecr_df.drop(["URL", "METHOD", "K", "K'", "EVAL", "EC", "ENZYME NAME",
                      "REACTION IN KEGG IDS"],
                     axis=1, inplace=True)

        logging.debug('Successfully added %d reactions from TECRDB' % 
                      tecr_df.shape[0])
        return tecr_df
        
    @staticmethod
    def read_formations(resource_name, weight):
        """
            Read the Formation Energy data
        """
        
        formation_df = pd.read_csv(resource_stream('component_contribution',
                                                   '/data/formation_energies_transformed.csv'))

        cids_that_dont_decompose = set(
            formation_df.loc[formation_df['decompose'] == 0, 'cid'])

        for col in ["dG'0", "T", "I", "pH", "pMg"]:
            formation_df[col] = formation_df[col].apply(float)

        formation_df = formation_df[~pd.isnull(formation_df["dG'0"])]
        formation_df['reaction'] = formation_df['cid'].apply(
                lambda c: KeggReaction({c: 1}))

        formation_df['weight'] = weight
        formation_df['balance'] = False
        formation_df['description'] = formation_df['name'] + ' formation'
        formation_df.rename(columns={'compound_ref': 'reference'}, inplace=True)
        formation_df.drop(['name', 'cid', 'remark', 'decompose'],
                          axis=1, inplace=True)

        logging.debug('Successfully added %d formation energies' % 
                      formation_df.shape[0])
        return formation_df, cids_that_dont_decompose
        
    @staticmethod
    def read_redox(resource_name, weight):
        """
            Read the Reduction potential data
        """
        redox_df = pd.read_csv(resource_stream('component_contribution',
                                               '/data/redox.csv'))
        
        delta_nH = redox_df['nH_red'] - redox_df['nH_ox']
        delta_charge = redox_df['charge_red'] - redox_df['charge_ox']
        delta_e = delta_nH - delta_charge
        redox_df["dG'0"] = -F * redox_df["E'0"] * delta_e
        redox_df['reaction'] = [KeggReaction({row['CID_ox'] : -1, row['CID_red'] : 1})
                                for _, row in redox_df.iterrows()]
        redox_df['weight'] = weight
        redox_df['balance'] = False
        redox_df['description'] = redox_df['name'] + ' redox'
        redox_df.rename(columns={'ref': 'reference'}, inplace=True)
        redox_df.drop(['name', 'CID_ox', 'CID_red', 'charge_ox', 'charge_red',
                       'nH_ox', 'nH_red', "E'0"], axis=1, inplace=True)

        logging.debug('Successfully added %d redox potentials' % 
                      redox_df.shape[0])
        return redox_df
    
    @staticmethod
    def get_all_thermo_params():
        base_path = os.path.split(os.path.realpath(__file__))[0]
    
        fname, weight = TrainingData.FNAME_DICT['TECRDB']
        fname = os.path.join(base_path, fname)
        tecr_df = TrainingData.read_tecrdb(fname, weight)
        
        fname, weight = TrainingData.FNAME_DICT['FORMATION']
        fname = os.path.join(base_path, fname)
        formation_df, cids_that_dont_decompose = TrainingData.read_formations(fname, weight)
        
        fname, weight = TrainingData.FNAME_DICT['REDOX']
        fname = os.path.join(base_path, fname)
        redox_df = TrainingData.read_redox(fname, weight)
        
        thermo_df = pd.concat([tecr_df, formation_df, redox_df], sort=False)
        return thermo_df, cids_that_dont_decompose
    
    def balance_reactions(self, rxn_inds_to_balance):
        """
            use the chemical formulas from the InChIs to verify that each and every
            reaction is balanced
        """
        elements, Ematrix = self.ccache.get_element_matrix(
            map(lambda s: 'KEGG:' + s, self.cids))
        cpd_inds_without_formula = list(np.nonzero(np.any(np.isnan(Ematrix), 1))[0].flat)
        Ematrix[np.isnan(Ematrix)] = 0

        S_without_formula = self.S[cpd_inds_without_formula, :]
        rxn_inds_without_formula = np.nonzero(np.any(S_without_formula != 0, 0))[0]
        rxn_inds_to_balance = set(rxn_inds_to_balance).difference(rxn_inds_without_formula)

        # need to check that all elements are balanced (except H, but including e-)
        # if only O is not balanced, add water molecules
        if 'O' in elements:
            i_H2O = self.cids.index('C00001')
            j_O = elements.index('O')
            conserved = np.dot(Ematrix.T, self.S)
            for k in rxn_inds_to_balance:
                self.S[i_H2O, k] = self.S[i_H2O, k] - conserved[j_O, k]

        # recalculate conservation matrix
        conserved = Ematrix.T * self.S
        
        rxn_inds_to_remove = [k for k in rxn_inds_to_balance 
                              if np.any(conserved[:, k] != 0, 0)]
        
        for k in rxn_inds_to_remove:
            sprs = {}
            for i in np.nonzero(self.S[:, k])[0]:
                sprs[self.cids[i]] = self.S[i, k]
            reaction = KeggReaction(sprs)
            logging.debug('unbalanced reaction #%d: %s' %
                          (k, reaction.write_formula()))
            for j in np.where(conserved[:, k])[0].flat:
                logging.debug('there are %d more %s atoms on the right-hand side' %
                              (conserved[j, k], elements[j]))
        
        rxn_inds_to_keep = \
            set(range(self.S.shape[1])).difference(rxn_inds_to_remove)
        
        rxn_inds_to_keep = sorted(rxn_inds_to_keep)
        
        self.S = self.S[:, rxn_inds_to_keep]
        self.dG0_prime = self.dG0_prime[rxn_inds_to_keep]
        self.T = self.T[rxn_inds_to_keep]
        self.I = self.I[rxn_inds_to_keep]
        self.pH = self.pH[rxn_inds_to_keep]
        self.pMg = self.pMg[rxn_inds_to_keep]
        self.weight = self.weight[rxn_inds_to_keep]
        self.reference = [self.reference[i] for i in rxn_inds_to_keep]
        self.description = [self.description[i] for i in rxn_inds_to_keep]

        logging.debug('After removing %d unbalanced reactions, the stoichiometric '
                      'matrix contains: '
                      '%d compounds and %d reactions' %
                      (len(rxn_inds_to_remove), self.S.shape[0], self.S.shape[1]))

    def reverse_transform(self):
        """
            Calculate the reverse transform for all reactions in training_data.
        """
        n_rxns = self.S.shape[1]
        reverse_ddG0 = np.zeros(n_rxns)
        self.I[np.isnan(self.I)] = 0.25 # default ionic strength is 0.25M
        self.pMg[np.isnan(self.pMg)] = 14 # default pMg is 14
        for i in range(n_rxns):
            for j in np.nonzero(self.S[:, i])[0]:
                cid = self.cids[j]
                if cid == 'C00080': # H+ should be ignored in the Legendre transform
                    continue
                comp = self.ccache.get_compound(cid)
                ddG0 = comp.transform_p_h_7(self.pH[i], self.I[i], self.T[i])
                reverse_ddG0[i] = reverse_ddG0[i] + ddG0 * self.S[j, i]

        self.dG0 = self.dG0_prime - reverse_ddG0
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=
        'Prepare all thermodynamic training data in a .mat file for running '
        'component contribution.')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                       help='the path to the .mat file that should be written '
                       'containing the training data')
    
    args = parser.parse_args()
    td = TrainingData()
    td.savemat(args.outfile)