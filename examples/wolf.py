from subprocess import Popen, PIPE
import numpy as np

REACTION_FNAME = 'wolf_reactions.txt'
PYTHON_BIN = 'python'
PYTHON_SCRIPT_FNAME = '../python/component_contribution.py'

p1 = Popen(['cat', REACTION_FNAME], stdout=PIPE)
p2 = Popen([PYTHON_BIN, PYTHON_SCRIPT_FNAME], stdin=p1.stdout,
           executable=PYTHON_BIN, stdout=PIPE)
res = p2.communicate()[0]
model_dG0 = np.array(res)
print str(model_dG0)