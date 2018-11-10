import logging
import numpy as np
import pandas as pd
from scipy.io import savemat
from .thermodynamic_constants import R, F
from .compound_cache import ccache
from . import Reaction
from pkg_resources import resource_stream

class TrainingData(object):
    
    def __init__(self):
        self.S = None # a DataFrame containing the stoichiometric matrix
        self.reaction_df = None # a DataFrame containing all the reaction data 
    
    @property
    def stoichiometric_matrix(self):
        return self.S.values
    
    @property
    def cids(self):
        return list(self.S.index)
    
    @property
    def dG0(self):
        return self.reaction_df['dG0'].values

    @property
    def weight(self):
        return self.reaction_df['weight'].values

    def create_stoichiometric_matrix_from_reactions(self, reactions):
        cids = set(['C00001', 'C00080'])
        for rxn in reactions:
            cids = cids.union(rxn.keys())
        cids = sorted(map(lambda s: 'KEGG:' + s, cids))
        
        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids' and the columns match the reaction indices
        self.S = pd.DataFrame(index=cids,
                              columns=reactions.index,
                              dtype=float).fillna(0)
        
        for k, rxn in reactions.items():
            for cid, coeff in rxn.items():
                self.S.at['KEGG:' + cid, k] = coeff

    def balance_reactions(self):
        """
            use the chemical formulas from the InChIs to verify that each and
            every reaction is balanced
        """

        element_df = ccache.get_element_data_frame(self.cids)
        
        # find all reactions that contain only compounds that have formulae
        cpd_with_formulae = (element_df != 0).any(axis=1)
        logging.info('# compounds without a formula: %d'
                     % sum(~cpd_with_formulae))
        
        rxn_with_formulae = (self.S.loc[~cpd_with_formulae, :] == 0).all(axis=0)
        logging.info('# reactions with full formulae: %d'
                     % sum(rxn_with_formulae))

        # recalculate final conservation matrix
        to_balance = self.reaction_df['balance'].copy()
        logging.info('# reactions we need to check for balacne: %d'
                     % to_balance.sum())
        
        to_balance = to_balance & rxn_with_formulae
        logging.info('# -> of which also have a formulae: %d'
                     % to_balance.sum())

        # balance O atoms using water
        self.S.loc['KEGG:C00001', to_balance] -= element_df['O'].T @ self.S.loc[:, to_balance]

        # balance H atoms using protons
        self.S.loc['KEGG:C00080', to_balance] -= element_df['H'].T @ self.S.loc[:, to_balance]
        
        imbalance_matrix = element_df.T @ self.S
        to_remove = to_balance & imbalance_matrix.any(axis=0)
        logging.info('# --> of which are not balanced and should be removed: %d'
                     % to_remove.sum())
        
        if to_remove.sum() > 0:
            for i, row in self.S.loc[:, to_remove].T.iterrows():
                sprs = {cid: coeff for cid,coeff in row.items() if coeff != 0}
                reaction = Reaction(sprs)
                logging.warning('unbalanced reaction #%s: %s' %
                                (i, reaction.write_formula()))
                for j, v in imbalance_matrix[i].items():
                    logging.warning('there are %d more %s atoms on the right-hand side' %
                                  (v, j))
            self.S = self.S.loc[:, ~to_remove]
            self.S.columns = range(self.S.shape[1])
        
        self.reaction_df = self.reaction_df.loc[self.S.columns, :]
        
        # now get rid of the protons, since we are applying Alberty's 
        # framework where their potential is set to 0, and the pH is held
        # as a controlled parameter
        self.S.drop('KEGG:C00080', axis=0, inplace=True)

        logging.info('After removing %d unbalanced reactions, the stoichiometric '
                     'matrix contains: '
                     '%d compounds and %d reactions' %
                     (sum(to_remove), self.S.shape[0], self.S.shape[1]))

    def iterreactions(self):
        for i, row in self.S.T.iterrows():
            sprs = {cid: coeff for cid,coeff in row.items() if coeff != 0}
            rxn = Reaction(sprs)
            yield i, rxn

    def reverse_transform(self):
        """
            Calculate the reverse transform for all reactions in training_data.
        """
        self.reaction_df['dG0'] = self.reaction_df['dG0_prime']
        
        for i, rxn in self.iterreactions():
            aq_cond = self.reaction_df.loc[i, ['pH', 'I', 'T']]
            self.reaction_df.at[i, 'dG0'] -= rxn.get_transform_ddG0(*aq_cond)

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
             'S': self.S.values}
        savemat(file_name, d, oned_as='row')

    def to_data_frame(self):
        df = self.reaction_df.copy()
        
        for i, rxn in self.iterreactions():
            df.at[i, 'formula'] = str(rxn)

        df.index.name = 'index'
        return df

class ToyTrainingData(TrainingData):
    
    def __init__(self):
        super().__init__()
        
        self.reaction_df = pd.read_csv(resource_stream('component_contribution',
                                                       '/data/toy_training_data.csv'))
        reactions = self.reaction_df['KEGG reaction'].apply(
                lambda r: Reaction.parse_formula(r, arrow='='))
        self.create_stoichiometric_matrix_from_reactions(reactions)
        self.reverse_transform()

class FullTrainingData(TrainingData):

    def __init__(self):
        super().__init__()
        logging.info('Reading the training data files')
        self.reaction_df, cids_that_dont_decompose = self.get_all_thermo_params()
        
        self.create_stoichiometric_matrix_from_reactions(self.reaction_df['reaction'])
        
        cids_that_dont_decompose.update(['C00001', 'C00080'])
        self.cids_that_dont_decompose = list(map(lambda s: 'KEGG:' + s,
                                                 cids_that_dont_decompose))
        
        logging.info('Balancing reactions in the training dataset with H2O and H+')
        self.balance_reactions()
        
        logging.info('Applying the reverse Legendre transform on all dG0_prime values')
        self.reverse_transform()

    @staticmethod
    def read_tecrdb():
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
            Reaction.parse_formula(r, arrow='='))

        tecr_df.rename(columns={'REF_ID': 'reference',
                                'REACTION IN COMPOUND NAMES': 'description'},
                       inplace=True)
        tecr_df['balance'] = True
        tecr_df.drop(["URL", "METHOD", "K", "K'", "EVAL", "EC", "ENZYME NAME",
                      "REACTION IN KEGG IDS"],
                     axis=1, inplace=True)

        logging.debug('Successfully added %d reactions from TECRDB' % 
                      tecr_df.shape[0])
        return tecr_df
        
    @staticmethod
    def read_formations():
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
                lambda c: Reaction({c: 1}))

        formation_df['balance'] = False
        formation_df['description'] = formation_df['name'] + ' formation'
        formation_df.rename(columns={'compound_ref': 'reference'}, inplace=True)
        formation_df.drop(['name', 'cid', 'remark', 'decompose'],
                          axis=1, inplace=True)

        logging.debug('Successfully added %d formation energies' % 
                      formation_df.shape[0])
        return formation_df, cids_that_dont_decompose
        
    @staticmethod
    def read_redox():
        """
            Read the Reduction potential data
        """
        redox_df = pd.read_csv(resource_stream('component_contribution',
                                               '/data/redox.csv'))
        
        delta_nH = redox_df['nH_red'] - redox_df['nH_ox']
        delta_charge = redox_df['charge_red'] - redox_df['charge_ox']
        delta_e = delta_nH - delta_charge
        redox_df["dG'0"] = -F * redox_df["E'0"] * delta_e
        redox_df['reaction'] = [Reaction({row['CID_ox'] : -1, row['CID_red'] : 1})
                                for _, row in redox_df.iterrows()]
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
        tecr_df = FullTrainingData.read_tecrdb()
        tecr_df['weight'] = 1.0
        
        formation_df, cids_that_dont_decompose = FullTrainingData.read_formations()
        formation_df['weight'] = 1.0
        
        redox_df = FullTrainingData.read_redox()
        redox_df['weight'] = 1.0
        
        reaction_df = pd.concat([tecr_df, formation_df, redox_df], sort=False)
        reaction_df.reset_index(drop=True, inplace=True)

        reaction_df['I'].fillna(0.25, inplace=True) # default ionic strength is 0.25M
        reaction_df['pMg'].fillna(14, inplace=True) # default pMg is 14
        
        reaction_df.rename(columns={"dG'0": 'dG0_prime'}, inplace=True)

        return reaction_df, cids_that_dont_decompose
    
        
if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    import argparse
    parser = argparse.ArgumentParser(description=
        'Prepare all thermodynamic training data in a .mat file for running '
        'component contribution.')
    parser.add_argument('-t', action='store_true', dest='toy',
                       help='use the toy example rather than the full database')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                       help='the path to the .csv file that should be written '
                       'containing the training data')
    
    args = parser.parse_args()
    if args.toy:
        td = ToyTrainingData()
    else:
        td = FullTrainingData()
    td.to_data_frame().to_csv(args.outfile)
