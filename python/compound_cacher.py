import json, os, logging, csv
from singletonmixin import Singleton
from compound import Compound

CACHE_FNAME = '../cache/compounds.json'
KEGG_ADDITIONS_TSV_FNAME = '../data/kegg_additions.tsv'

class CompoundEncoder(json.JSONEncoder):
    def default(self, obj):
        if (isinstance(obj, Compound)):
            return obj.to_json_dict()
        return json.JSONEncoder.default(self, obj)

class CompoundCacher(Singleton):
    def __init__(self):
        self.compound_dict = {}
        self.need_to_update_cache_file = False

        if os.path.exists(CACHE_FNAME):
            for d in json.load(open(CACHE_FNAME, 'r')):
                self.compound_dict[d['id']] = Compound.from_json_dict(d)
        else:
            path = os.path.dirname(CACHE_FNAME)
            try:
                os.makedirs(path)
            except OSError:
                if not os.path.isdir(path):
                    raise Exception('Cannot create cache directory: ' + path)
        self.get_kegg_additions()
    
    def get_kegg_additions(self):
        # fields are: name, cid, inchi
        for row in csv.DictReader(open(KEGG_ADDITIONS_TSV_FNAME, 'r'), delimiter='\t'):
            cid = int(row['cid'])
            compound_id = 'C%05d' % cid
            if compound_id not in self.compound_dict:
                logging.info('Calculating pKa for additional compound: ' + compound_id)
                inchi = row['inchi']
                pKas, majorMSpH7, nHs, zs = Compound.get_species_pka(inchi)
                comp = Compound('KEGG', compound_id, inchi,
                                pKas, majorMSpH7, nHs, zs)
                self.compound_dict[comp.compound_id] = comp
                self.need_to_update_cache_file = True            

    def dump(self):
        if self.need_to_update_cache_file:
            fp = open(CACHE_FNAME, 'w')
            json.dump(self.compound_dict.values(), fp, cls=CompoundEncoder)
            fp.close()
            self.need_to_update_cache_file = False
        
    def get_kegg_compound(self, cid):
        compound_id = 'C%05d' % cid
        if compound_id in self.compound_dict:
            return self.compound_dict[compound_id]
        
        logging.info('Downloading structure and calculating pKa for: ' + compound_id)
        comp = Compound.from_kegg(cid)
        self.compound_dict[comp.compound_id] = comp
        self.need_to_update_cache_file = True
        return comp
        
if __name__ == '__main__':
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    ccache = CompoundCacher.getInstance()
    comp = ccache.get_kegg_compound(25)
    print comp.pKas
    ccache.dump()
