import sys, logging
sys.path.append('../python')
from compound import Compound

logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)

for cid in [282]:
    comp = Compound.from_kegg(cid)
    print comp
    print "C%05d: ddG0_prime = %.2f" % (cid, comp.transform(7.0, 0.2, 300))
