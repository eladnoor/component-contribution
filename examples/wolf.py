import sys, logging
sys.path.append('../python')
from kegg_model import KeggModel
from training_data import TrainingData
from component_contribution import ComponentContribution
from scipy.io import savemat

logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
td = TrainingData()
model = KeggModel.load_kegg_model('../examples/wolf_reactions.txt')
cc = ComponentContribution(model, td)

#savemat('../examples/groups.mat', {'G' : cc.G, 'dG0_cc' : cc.dG0}, oned_as='row')
