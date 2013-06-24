import sys, logging
sys.path.append('../python')
from compound import Compound

logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)

for cid in [2]:
    comp = Compound.from_kegg(cid)
    print comp
    print "C%05d: ddG0_prime = %.2f" % (cid, comp.transform(7.0, 0.2, 300))
    print "C%05d: ddG0_prime(z=0) = %.2f" % (cid, comp.transform_neutral(7.0, 0.2, 300))
    
