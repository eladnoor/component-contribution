import sys
sys.path.append('../python')
from kegg_model import KeggModel

model = KeggModel.load_kegg_model('../examples/wolf_reactions.txt')
print model.S
print model.cids
