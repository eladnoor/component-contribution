import six

import numpy

from weakref import ref as weakref

from pandas import DataFrame, read_csv, concat

from component_contribution.compound import Compound
from component_contribution.singleton import Singleton


class CompoundCache(Singleton):

    COLUMNS = ["inchi", "name", "atom_bag", "p_kas", "smiles", "major_ms", "number_of_protons", "charges"]

    """
    CompoundCache is a singleton that handles caching of Compound objects for the component-contribution package. 
    The Compounds are retrieved by their ID (e.g., KEGG COMPOUND ID, ChEBI Id or HMDB in most cases) or InChI Key.
    The first time a Compound is requested, it is obtained from the relevant database and a Compound object is 
    created (this takes a while because it usually involves internet communication and then invoking the ChemAxon
    plugin for calculating the pKa values for that structure). Any further request for the same Compound ID will 
    draw the object from the cache. When the method dump() is called, all cached data is written to a file that 
    will be loaded in future python sessions.
    """

    @staticmethod
    def _read_atom_bag(serialized):
        return {atom: int(count) for atom_count in serialized.split(";") for atom, count in atom_count.split(":")}

    @staticmethod
    def _serialize_atom_bag(atom_bag):
        return ";".join("%s:%i" % (atom, count) for atom, count in six.iteritems(atom_bag))

    @staticmethod
    def _read_int_list(serialized):
        return [int(i) for i in serialized.split(";")]

    @staticmethod
    def _serialize_int_list(int_list):
        return ";".join(str(i) for i in int_list)

    @staticmethod
    def _serialize_cross_refs(data, inchi_key_to_compound_ids):
        return [";".join(inchi_key_to_compound_ids[inchi_key]) for inchi_key in data.index]

    @staticmethod
    def _read_cross_refs(data):
        compound_id_to_inchi_key = {}
        for inchi_key, row in data.iterrows():
            compound_ids = row.cross_references
            for compound_id in compound_ids.split(";"):
                compound_id_to_inchi_key[compound_id] = inchi_key

        return compound_id_to_inchi_key

    def __init__(self):
        self._data = DataFrame(columns=self.COLUMNS)
        self._compound_id_to_inchi_key = {}
        self._inchi_key_to_compound_ids = {}
        self.compound_dict = {}

    def get_all_compound_ids(self):
        return sorted(self._data.index)

    def load(self, file_name):
        """
        Loads a cache dump from a file.

        Parameters
        ----------
        file_name : str
            The file name.

        """
        data = read_csv(file_name, index_col=0)
        data['atom_bag'] = data.atom_bag.apply(self._read_atom_bag)
        data["number_of_protons"] = data.number_of_protons.apply(self._read_int_list)
        data["charges"] = data.charges.apply(self._serialize_int_list)

        del data['cross_references']
        self._data = concat([self._data, data])
        compound_id_to_inchi_key = self._read_cross_refs(data)

        self._compound_id_to_inchi_key.update(compound_id_to_inchi_key)

    def dump(self, file_name):
        """
        Writes the cache to a file.

        Parameters
        ----------
        file_name : str
            The file name.

        """
        to_write = DataFrame(self._data)
        to_write['atom_bag'] = to_write.atom_bag.apply(self._serialize_atom_bag)
        to_write["number_of_protons"] = to_write.number_of_protons.apply(self._serialize_int_list)
        to_write["charges"] = to_write.charges.apply(self._serialize_int_list)
        to_write["cross_references"] = self._serialize_cross_refs(to_write, self._inchi_key_to_compound_ids)
        to_write.to_csv(file_name)

    def get_compound(self, compound_id):
        if compound_id in self._compound_id_to_inchi_key:  # compound exists
            return self.get(self._compound_id_to_inchi_key[compound_id])
        elif compound_id in self._data.index:  # compound_id is an InChI Key and exists
            return self.get(compound_id)
        else:  # compound does not exist.
            cpd = Compound.get(compound_id)
            if cpd.inchi_key in self._data.index:
                cpd = self.get(cpd.inchi_key)
                self._compound_id_to_inchi_key[compound_id] = cpd.inchi_key
            else:
                self.add(cpd)

            return cpd
    
    def remove(self, inchi_key):
        if inchi_key not in self._data.index:
            raise KeyError(inchi_key)

        del self.compound_dict[inchi_key]
        del self._inchi_key_to_compound_ids
        self._data.drop(inchi_key)

        for compound_id in self._compound_id_to_inchi_key:
            if self._compound_id_to_inchi_key[compound_id] == inchi_key:
                del self._compound_id_to_inchi_key[compound_id]

    def get(self, inchi_key):
        if inchi_key not in self._data.index:
            raise KeyError(inchi_key)

        if inchi_key in self.compound_dict:
            cpd = self.compound_dict[inchi_key]
        else:
            data = self._data.loc[inchi_key]
            cpd = Compound(inchi_key, data.inchi, data.atom_bag, data.p_kas, data.smiles,
                           data.major, data.number_of_protons, data.charges)
            self.compound_dict[inchi_key] = cpd

        return cpd
    
    def add(self, cpd):
        if cpd.inchi_key not in self._inchi_key_to_compound_ids:
            self._inchi_key_to_compound_ids[cpd.inchi_key] = set()

        self._inchi_key_to_compound_ids[cpd.cpd.inchi_key].add(cpd.compound_id)
        self._compound_id_to_inchi_key[cpd.compound_id] = cpd.inchi_key
        self.compound_dict[cpd.inchi_key] = weakref(cpd)
        self._data.loc[cpd.inchi_key] = [cpd.inchi, cpd.name, cpd.atom_bag, cpd.p_kas, cpd.smiles,
                                         cpd.major_ms, cpd.number_of_protons, cpd.charges]

    def get_element_matrix(self, compound_ids):
        if isinstance(compound_ids, six.string_types):
            compound_ids = [compound_ids]
        # gather the "atom bags" of all compounds in a list 'atom_bag_list'
        elements = set()
        atom_bag_list = []
        for compound_id in compound_ids:
            cpd = self.get_compound(compound_id)
            atom_bag = cpd.atom_bag
            if atom_bag is not None:
                elements = elements.union(atom_bag.keys())
            atom_bag_list.append(atom_bag)
        elements.discard('H')  # don't balance H (it's enough to balance e-)
        elements = sorted(elements)
        
        # create the elemental matrix, where each row is a compound and each
        # column is an element (or e-)
        element_matrix = numpy.matrix(numpy.zeros((len(atom_bag_list), len(elements))))
        for i, atom_bag in enumerate(atom_bag_list):
            if atom_bag is None:
                element_matrix[i, :] = numpy.nan
            else:
                for j, elem in enumerate(elements):
                    element_matrix[i, j] = atom_bag.get(elem, 0)
        return elements, element_matrix
