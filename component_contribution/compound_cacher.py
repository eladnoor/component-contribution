import json, os, logging, csv, gzip, numpy
from compound import Compound
base_path = os.path.split(os.path.realpath(__file__))[0]

### Input Files:
# original version of the KEGG compound file 
OLD_COMPOUND_JSON_FNAME = os.path.join(base_path, '../data/equilibrator_compounds.json.gz')

# a CSV file with additional names and InChIs (mostly compounds missing from KEGG
# and added manually)
KEGG_ADDITIONS_TSV_FNAME = os.path.join(base_path, '../data/kegg_additions.tsv')

### Files created by this module:
# names and InChIs only
KEGG_COMPOUND_JSON_FNAME = os.path.join(base_path, '../data/kegg_compounds.json.gz') 

# names, InChIs and pKa data
DEFAULT_CACHE_FNAME = os.path.join(base_path, '../cache/compounds.json.gz') 


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
        self.cache_fname = cache_fname
        if self.cache_fname is None:
            self.cache_fname = DEFAULT_CACHE_FNAME
        
        compounds = json.load(gzip.open(KEGG_COMPOUND_JSON_FNAME, 'r'))
        self.compound_id2inchi = { d['compound_id']: d['inchi'] 
                                   for d in compounds }
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
            dict_data = map(lambda x:x.to_json_dict(), data)
            json.dump(dict_data, fp, cls=CompoundEncoder, 
                      sort_keys=True, indent=4,  separators=(',', ': '))
            fp.close()
            self.need_to_update_cache_file = False
        
    def get_compound(self, compound_id):
        if compound_id not in self.compound_dict:
            logging.debug('Cache miss: %s' % str(compound_id))
            inchi = self.compound_id2inchi[compound_id]
            comp = Compound.from_inchi('KEGG', compound_id, inchi)
            self.add(comp)

        logging.debug('Cache hit: %s' % str(compound_id))
        return self.compound_dict[compound_id]

    def remove(self, compound_id):
        if compound_id in self.compound_dict:
            del self.compound_dict[compound_id]
        else:
            logging.debug('%s is not cached, cannot remove it' % str(compound_id))
    
    def add(self, comp):
        self.compound_dict[comp.compound_id] = comp
        self.need_to_update_cache_file = True
            
    def get_element_matrix(self, compound_ids):
        if type(compound_ids) == str:
            compound_ids = [compound_ids]
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
        Ematrix = numpy.matrix(numpy.zeros((len(atom_bag_list), len(elements))))
        for i, atom_bag in enumerate(atom_bag_list):
            if atom_bag is None:
                Ematrix[i, :] = numpy.nan
            else:
                for j, elem in enumerate(elements):
                    Ematrix[i, j] = atom_bag.get(elem, 0)
        return elements, Ematrix

###############################################################################

    @staticmethod
    def RebuildCompoundJSON():
        
        kegg_dict = {}
        for d in json.load(gzip.open(OLD_COMPOUND_JSON_FNAME, 'r')):
            cid = d['CID']
            kegg_dict[cid] = {'compound_id': cid,
                              'name': d['name'],
                              'names': d['names'],
                              'inchi': d['InChI']}
        
        # override some of the compounds or add new ones with 'fake' IDs,
        # i.e. C80000 or higher.
        for d in csv.DictReader(open(KEGG_ADDITIONS_TSV_FNAME, 'r'),
                                delimiter='\t'):
            cid = 'C%05d' % int(d['cid'])
            kegg_dict[cid] = {'compound_id': cid,
                              'name': d['name'],
                              'names': [d['name']],
                              'inchi': d['inchi']}
        
        compound_json = [kegg_dict[compound_id] for compound_id in sorted(kegg_dict.keys())]

        new_json = gzip.open(KEGG_COMPOUND_JSON_FNAME, 'w')
        json.dump(compound_json, new_json, sort_keys=True, indent=4)
        new_json.close()
    
###############################################################################

    @staticmethod
    def BuildCache(start_from_scratch=False):
        if start_from_scratch and os.path.exists(DEFAULT_CACHE_FNAME):
            os.remove(DEFAULT_CACHE_FNAME)
            
        ccache = CompoundCacher(cache_fname=DEFAULT_CACHE_FNAME)
        
        i = 0
        for compound_id in ccache.get_all_compound_ids():
            logging.debug('Caching %s' % compound_id)
            comp = ccache.get_compound(compound_id)
            logging.debug(str(comp))
            i += 1
            if i % 100 == 0:
                logging.debug('Dumping Cache ...')
                ccache.dump()
        
        ccache.dump()

###############################################################################
        
if __name__ == '__main__':
    logger = logging.getLogger('')
    logger.setLevel(logging.WARNING)
    
    CompoundCacher.RebuildCompoundJSON()
    CompoundCacher.BuildCache()
