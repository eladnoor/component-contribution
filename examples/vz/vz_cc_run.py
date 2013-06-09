import os	
os.chdir('/home/vitoz/Dropbox/Rabino/git/component-contribution/python')

import sys
import numpy as np
from inchi2gv import GroupsData, InChI2GroupVector, GROUP_CSV, GroupDecompositionError
from training_data import TrainingData
from kegg_model import KeggModel
from compound_cacher import CompoundCacher
from component_contribution import ComponentContribution
import csv 

##

with open('../examples/vz/StoiMat.tsv','r') as tsv:
    tabinput = [line.strip().split('\t') for line in tsv]

reactions = tabinput[0]
reactions = reactions[1:]
tabinput = tabinput[1:]

S = np.zeros((len(tabinput), len(tabinput[0])-1))
cids = list()
for i in range(len(tabinput)):
    cids.append(float(tabinput[i][0].replace('"','')))
    for j in range(len(tabinput[0])-1):
        S[i,j] = tabinput[i][j+1]

np.array()
model = KeggModel(S, cids)
model = model.check_S_balance()
td = TrainingData()
cc = ComponentContribution(td)
model_dG0, model_cov_dG0 = cc.estimate_kegg_model(model)
model_dG0_prime = model_dG0 + model.get_transform_ddG0(pH=7.5, I=0.2, T=298.15)

dG0_std = np.sqrt(model_cov_dG0.diagonal())

out_file = '../examples/vz/output/cc_dG.tsv'
outf = open(out_file, "w")
outf.write("reaction\tdGr\tdGrSD\n")
for i in range(len(reactions)):
    outf.write(reactions[i] + "\t" + str(float(model_dG0_prime[i])) + "\t" +str(float(dG0_std[0,i]))+ "\n")

outf.close()

out_file = '../examples/vz/output/cc_dG_cov.tsv'
outf = open(out_file, "w")
outf.write("reaction\tdGr\tdGrSD\n")
for i in range(len(reactions)):
    for j in range(len(reactions)):
        if j+1 < len(reactions):
            outf.write(str(model_cov_dG0[i,j]) + "\t" )
        else:
            outf.write(str(model_cov_dG0[i,j]) + "\n" )

outf.close()