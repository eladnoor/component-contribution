import sys, logging
sys.path.append('../python')
from kegg_model import KeggModel
from training_data import TrainingData
from component_contribution import ComponentContribution
from compound import Compound

logger = logging.getLogger('')
logger.setLevel(logging.INFO)

td = TrainingData()
cc = ComponentContribution(td)

model = KeggModel.load_kegg_model('../examples/wolf_reactions.txt')
model_dG0, model_cov_dG0 = cc.estimate_kegg_model(model)

