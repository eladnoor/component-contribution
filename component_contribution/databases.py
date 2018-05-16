import requests
import pybel
from component_contribution.singleton import Singleton


class Databases(Singleton):
    _registry = {}

    @classmethod
    def register(cls, database_name, mol_retrieving_function):
        cls._registry[database_name] = mol_retrieving_function

    @classmethod
    def get_molecule(cls, database_name, accession):
        if database_name in cls._registry:
            return cls._registry[database_name](accession)
        else:
            raise KeyError(database_name)


def get_kegg_molecule(accession):
    response = requests.get('http://rest.kegg.jp/get/cpd:%s/mol' % accession)
    try:
        response.raise_for_status()
        return pybel.readstring("mol", response.text)
    except requests.HTTPError as e:
        exception = KeyError(accession)
        exception.__cause__ = e
        raise exception


def get_hmdb_molecule(accession):
    accession = "HMDB%07d" % int(accession.replace("HMDB", ""))
    response = requests.get("http://www.hmdb.ca/structures/metabolites/%s.mol" % accession)
    try:
        response.raise_for_status()
        return pybel.readstring("mol", response.text)
    except requests.HTTPError as e:
        exception = KeyError(accession)
        exception.__cause__ = e
        raise exception


def get_inchi_molecule(accession):
    return pybel.readstring("inchi", accession)


def get_chebi_molecule(accession):
    raise NotImplementedError


Databases.register("KEGG", get_kegg_molecule)
Databases.register("ChEBI", get_chebi_molecule)
Databases.register("HMDB", get_hmdb_molecule)
Databases.register("InChI", get_inchi_molecule)
