import json, os, logging, csv, gzip, sys
from compound import Compound
import numpy as np

DEFAULT_CACHE_FNAME = '../cache/compounds.json.gz'
COMPOUND_TSV_FNAME = '../data/compounds.tsv'

class CompoundEncoder(json.JSONEncoder):
    def default(self, obj):
        if (isinstance(obj, Compound)):
            return obj.to_json_dict()
        return json.JSONEncoder.default(self, obj)

class Singleton(type):
    def __init__(cls,name,bases,dic):
        super(Singleton,cls).__init__(name,bases,dic)
        cls.instance=None
    def __call__(cls,*args,**kw):
        if cls.instance is None:
            cls.instance=super(Singleton,cls).__call__(*args,**kw)
        return cls.instance
    
class CompoundCacher(object):
    """
        CompoundCacher is a singleton that handles caching of Compound objects
        for the component-contribution package. The Compounds are retrieved by
        their ID (which is the KEGG ID in most cases).
        The first time a Compound is requested, it is obtained from the relevant
        database and a Compound object is created (this takes a while because
        it usually involves internet communication and then invoking the ChemAxon
        plugin for calculating the pKa values for that structure).
        Any further request for the same Compound ID will draw the object from
        the cache. When the method dump() is called, all cached data is written
        to a file that will be loaded in future python sessions.
    """
    __metaclass__ = Singleton
    
    def __init__(self, cache_fname=None):
        base_path = os.path.split(os.path.realpath(__file__))[0]
        self.cache_fname = cache_fname
        if self.cache_fname is None:
            self.cache_fname = os.path.join(base_path, DEFAULT_CACHE_FNAME)
        
        compound_tsv_fname = os.path.join(base_path, COMPOUND_TSV_FNAME)
        
        compound_csv = csv.DictReader(open(compound_tsv_fname, 'r'),
                                      delimiter='\t')        
        self.compound_id2inchi = { d['compound_id']: d['inchi'] 
                                   for d in compound_csv }
        self.need_to_update_cache_file = False
        self.load()
    
    def get_all_compound_ids(self):
        return sorted(self.compound_id2inchi.keys())
    
    def load(self):
        # parse the JSON cache file and store in a dictionary 'compound_dict'
        self.compound_dict = {}
        self.compound_ids = []
        if os.path.exists(self.cache_fname):
            for d in json.load(gzip.open(self.cache_fname, 'r')):
                self.compound_ids.append(d['compound_id'])
                self.compound_dict[d['compound_id']] = Compound.from_json_dict(d)

    def dump(self):
        if self.need_to_update_cache_file:
            fp = gzip.open(self.cache_fname, 'w')
            data = sorted(self.compound_dict.values(),
                          key=lambda d:d.compound_id)
            json.dump(data, fp, cls=CompoundEncoder, 
                      sort_keys=True, indent=4,  separators=(',', ': '))
            fp.close()
            self.need_to_update_cache_file = False
        
    def get_compound(self, compound_id):
        if compound_id not in self.compound_dict:
            logging.debug('Cache miss: %s' % compound_id)
            inchi = self.compound_id2inchi[compound_id]
            comp = Compound.from_inchi('KEGG', compound_id, inchi)
            self.add(comp)

        logging.debug('Cache hit: %s' % compound_id)
        return self.compound_dict[compound_id]
            
    def add(self, comp):
        self.compound_dict[comp.compound_id] = comp
        self.need_to_update_cache_file = True
            
    def get_element_matrix(self, compound_ids):
        # gather the "atom bags" of all compounds in a list 'atom_bag_list'
        elements = set()
        atom_bag_list = []
        for compound_id in compound_ids:
            comp = self.get_compound(compound_id)
            atom_bag = comp.atom_bag
            if atom_bag is not None:
                elements = elements.union(atom_bag.keys())
            atom_bag_list.append(atom_bag)
        elements.discard('H') # don't balance H (it's enough to balance e-)
        elements = sorted(elements)
        
        # create the elemental matrix, where each row is a compound and each
        # column is an element (or e-)
        Ematrix = np.matrix(np.zeros((len(atom_bag_list), len(elements))))
        for i, atom_bag in enumerate(atom_bag_list):
            if atom_bag is None:
                Ematrix[i, :] = np.nan
            else:
                for j, elem in enumerate(elements):
                    Ematrix[i, j] = atom_bag.get(elem, 0)
        return elements, Ematrix

###############################################################################
        
if __name__ == '__main__':
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    
    if False:
        # this is legacy code which was used for creating the compound.tsv
        # file.
    
        COMPOUND_JSON_FNAME = 'data/kegg_compounds.json.gz'
        KEGG_ADDITIONS_TSV_FNAME = 'data/kegg_additions.tsv'
    
        kegg_dict = {}
        for d in json.load(gzip.open(COMPOUND_JSON_FNAME, 'r')):
            kegg_dict[d['CID']] = (d['name'], d['InChI'])
        
        for d in csv.DictReader(open(KEGG_ADDITIONS_TSV_FNAME, 'r'),
                                delimiter='\t'):
            kegg_dict[d['cid']] = (d['name'], d['inchi'])
        
        kegg_csv = csv.writer(open(COMPOUND_TSV_FNAME, 'w'), delimiter='\t')
        kegg_csv.writerow(['compound_id', 'name', 'inchi'])
        for compound_id in sorted(kegg_dict.keys()):
            name, inchi = kegg_dict[compound_id]
            kegg_csv.writerow([compound_id, name, inchi])
        sys.exit(0)
    
    ccache = CompoundCacher(cache_fname=None)
    
    i = 0
    for d in csv.DictReader(open(COMPOUND_TSV_FNAME, 'r'), delimiter='\t'):
        logging.info('Caching %s' % d['compound_id'])    
        comp = ccache.get_compound(d['compound_id'])
        logging.info(str(comp))
        i += 1
        if i % 100 == 0:
            logging.info('Dumping Cache ...')
            ccache.dump()
    
    ccache.dump()
