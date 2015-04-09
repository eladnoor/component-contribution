from compound_cacher import CompoundCacher
import csv, sys, re

def main(fname, pH, I, T):
    ccache = CompoundCacher()
    for row in csv.reader(open(fname, 'r'), delimiter='\t'):
        compound_id = re.findall('(C[0-9]+)_10', row[0])[0]
        dG0 = float(row[1])
        comp = ccache.get_compound(compound_id)
        dG0_prime = dG0 + comp.transform_neutral(pH, I, T)
        print '%s\t%f\t%f' % (compound_id, dG0, dG0_prime)
    ccache.dump()
        
if __name__ == '__main__':
    pH = 7
    I = 0.2
    T = 298.15
    main(sys.argv[1], pH, I, T)
