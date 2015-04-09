import csv
import numpy as np

class redox_carrier(object):

    # dG0 =  -E'*F * deltaE - R*T*ln(10)*pH * deltaH
    # Where: 
    #    F = 96.48 # kC/mol
    #    R*T*ln(10) = 5.7 kJ/mol
    #    deltaE - change in e-
    #    deltaH - change in H+
    #    pH - the conditions in which the E' was measured
    
    def __init__(self, cid_ox, cid_red, nH_ox, nH_red, z_ox, z_red, E0_prime, pH, ref):
        self.cid_ox = cid_ox
        self.cid_red = cid_red
        self.nH_ox = nH_ox
        self.nH_red = nH_red
        self.z_ox = z_ox
        self.z_red = z_red
        self.E0_prime = E0_prime
        self.pH = pH
        self.ref = ref
        self.delta_H = nH_red - nH_ox
        self.delta_e = (nH_red - nH_ox) - (z_red - z_ox) # difference in no. of electrons
        self.ddG0_prime = -E0_prime * F * self.delta_e
        
        # this calculation is not correct, one must use the reverse Lagendre transform
        # in order to convert G' to G.
        self.ddG0 = self.ddG0_prime - R * default_T * np.log(10) * pH * self.delta_H

class redox_carrier_dict(dict):
    
    def __init__(self):
        for row in csv.DictReader(open('../data/redox.tsv', 'r'), delimiter='\t'):
            name = row['name']
            cid_ox = int(row['CID_ox'])
            cid_red = int(row['CID_red'])
            nH_ox = int(row['nH_ox'])
            z_ox = int(row['charge_ox'])
            cids_ox = [val.cid_ox for val in self.values()]
            nH_red = int(row['nH_red'])
            z_red = int(row['charge_red'])
            E0_prime = float(row["E'0"])
            pH = float(row['pH'])
            ref = row['ref']
            self[name] = redox_carrier(cid_ox, cid_red, nH_ox, nH_red, 
                                       z_ox, z_red, E0_prime, pH, ref)

    def get_all_cids(self):
        return set(cids_ox + cids_red)
        cids_red = [val.cid_red for val in self.values()]