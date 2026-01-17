"""
Simplified identifier normalization for SPOKE - focuses on case normalization and simple heuristics
"""

import logging
import sys
from typing import Any, List, Dict, Tuple, Union

import requests
from biomapper2.core.normalizer import Normalizer


class SpokeIDNormalizer:
    """Simple, robust identifier normalization focused on case handling"""
    
    def __init__(self, biolink_version: str):
        self.biolink_version = biolink_version
        self.normalizer = Normalizer(biolink_version=biolink_version)
        self.curie_construction_map = self._load_curie_construction_map()
        self.spoke_underscore_prefixes = {'SNOMED_', 'CLO_', 'ENVO_', 'CHR_', 'HPS_', 'CVCL_', 'BFO_'}

    @staticmethod
    def _load_curie_construction_map() -> Dict[Tuple[str, str], Union[str, List[str]]]:
        curie_construction_map = {
            ('Anatomy', 'mesh_id'): 'mesh',
            ('Anatomy', 'uberon'): 'uberon',
            ('BiologicalProcess', 'go'): 'go',
            ('Blend', 'nhanes'): 'nhanes',
            ('CellLine', 'celllineontology'): 'clo',
            ('CellLine', 'clo'): 'clo',
            ('CellLine', 'cvcl'): 'cvcl',
            ('CellType', 'cl'): 'cl',
            ('CellularComponent', 'go'): 'go',
            ('ClinicalLab', 'unknown'): ('loinc', 'umls'),
            ('Complex', 'complexportal'): 'complexportal',
            ('Compound', 'chebi'): 'chebi',
            ('Compound', 'chembl_ids'): 'chembl.compound',
            ('Compound', 'chembl.compound'): 'chembl.compound',
            ('Compound', 'drugbank_ids'): 'drugbank',
            ('Compound', 'inchikey'): 'inchikey',
            ('Compound', 'kegg_compound_ids'): 'kegg.compound',
            ('Compound', 'kegg_drug_ids'): 'kegg.drug',
            ('Compound', 'pubchem_compound_ids'): 'pubchem.compound',
            # ('Compound', 'standardized_smiles'): 'smiles',  #  Some of these are incorrect; skipping for now
            ('Cytoband', 'unknown'): 'cytoband',
            ('DietarySupplement', 'nhanes'): 'nhanes',
            ('Disease', 'doid'): 'doid',
            ('Disease', 'icd9'): 'icd9',
            ('Disease', 'icd10'): 'icd10',
            ('Disease', 'mesh_list'): 'mesh',
            ('Disease', 'omim_list'): 'omim',
            ('Disease', 'snomedct'): 'snomedct',
            ('EC', 'explorenz'): 'ec',
            ('EC', 'metacyc'): 'metacyc.ec',
            ('Environment', 'bfo'): 'bfo',
            ('Environment', 'envo'): 'envo',
            ('ExtracellularParticle', 'vesiclepedia'): 'vesiclepedia',
            ('Gene', 'accession'): 'mirbase',
            ('Gene', 'chembl_id'): 'chembl.target',
            ('Gene', 'entrezgene'): 'ncbigene',
            ('Gene', 'mirbase'): 'ncbigene',  # Main 'identifier' given for these is ncbigene ID
            ('Haplotype', 'unknown'): ('pharmvar', 'dbsnp'),
            ('Location', 'geonames'): 'geonames',
            ('Location', 'unitedstateszipcode_database'): ('uszipcode', 'fips.place', 'fips.state'),
            ('MiRNA', 'accession'): 'mirbase',
            ('MiRNA', 'mirdb'): 'mirdb',
            ('MolecularFunction', 'go'): 'go',
            ('Organism', 'bv-brc'): 'bvbrc',
            ('Organism', 'ncbi-taxonomy'): 'ncbitaxon',
            ('Pathway', 'reactome'): 'react',
            ('Pathway', 'unknown'): 'metacyc.pathway',
            ('Pathway', 'wikipathways'): 'wikipathways',
            ('PharmacologicClass', 'fdaviadrugcentral'): ('ndfrt', 'mesh'),
            ('Protein', 'chembl_id'): 'chembl.target',
            ('Protein', 'uniprot'): 'uniprotkb',
            ('ProteinDomain', 'pfam'): 'pfam',
            ('ProteinFamily', 'pfam'): 'pfam',
            ('PwGroup', 'reactome'): 'react',
            ('Reaction', 'kegg'): 'kegg.reaction',
            ('Reaction', 'metacyc'): 'metacyc.reaction',
            ('Reaction', 'reactome'): 'react',
            ('SDoH', 'ahrqsdohdatabase'): 'ahrq',
            ('SDoH', 'cdc/atsdrsocialvulnerabilityindex'): 'cdcsvi',
            ('SDoH', 'chr'): 'chr',
            ('SDoH', 'hps'): 'hps',
            ('SDoH', 'mesh_ids'): 'mesh',
            ('SDoH', 'snomed'): 'snomedct',
            ('SideEffect', 'sider4.1'): 'umls',  # They give CUIs for nodes with SIDER source
            ('Symptom', 'hpo'): 'mesh',  # They give MeSH IDs for nodes with HPO source
            ('Symptom', 'icd9'): 'icd9',
            ('Symptom', 'icd10'): 'icd10',
            ('Symptom', 'mesh'): 'mesh',
            ('Symptom', 'snomedct'): 'snomedct',
            ('Variant', 'dbsnp'): 'dbsnp',
            ('Variant', 'unknown'): 'dbsnp',
        }
        return curie_construction_map


    def normalize_spoke_identifier(self, node_type: str, source: str, identifier: Any, properties: dict) -> Tuple[str, str]:
        """
        Converts SPOKE 'identifiers' and other miscellaneous ID properties to Biolink standardized curies (and IRIs).
        """
        identifier = str(identifier)
        spoke_prefix, local_id = self.get_curie_parts(identifier)

        # If the identifier was a curie, its prefix should override the source
        source = spoke_prefix if spoke_prefix else source

        source_cleaned = self.clean_name(source)
        lookup_key = (node_type, source_cleaned)

        if local_id and lookup_key in self.curie_construction_map:
            # Grab the curie construction info
            prefix_entry = self.curie_construction_map[lookup_key]
            if self.is_known_invalid_id(local_id, prefix_entry):
                logging.warning(f"Skipping known invalid ID for {prefix_entry}: {local_id}.")
                return KNOWN_INVALID, ''
            else:
                # Actually construct the curie
                curie_dict, _ = self.normalizer.get_curies({prefix_entry: local_id}, stop_on_invalid_id=False)
                if curie_dict:
                    curie, iri = next(iter(curie_dict.items()))  # There can only be one entry in here
                    return curie, iri

        return '', ''


    def get_curie_parts(self, identifier: str) -> Tuple[str, str]:
        num_colons = identifier.count(':')
        if num_colons == 0:
            if any(identifier.startswith(underscore_prefix) for underscore_prefix in self.spoke_underscore_prefixes):
                parts = identifier.split('_', 1)
                return parts[0].strip(), parts[1].strip()
            else:
                return '', identifier.strip()  # Some metacyc.ec IDs have trailing space
        elif num_colons == 1:
            parts = identifier.split(':')
            if parts[0].startswith('http'):  # This isn't a curie colon..
                return '', identifier.strip()
            else:
                return parts[0].strip(), parts[1].strip()
        else:
            logging.error(f"An identifier has more than one colon in it: {identifier}. Not sure what to do.")
            sys.exit(1)
    

    def extract_equivalent_identifiers(self, node_type: str, properties: Dict) -> List[str]:
        """Extract equivalent IDs from properties - simplified approach"""
        equivalent_ids = set()
        none_strings = {'null', 'none', 'nan'}
        
        # Figure out which properties probably contain an identifier
        relevant_properties = set()
        for property_name in properties.keys():
            if self.is_trusted_id_property(property_name):
                relevant_properties.add(property_name)

        # Construct proper curie(s) for each of those properties
        for id_prop_name in relevant_properties:
            id_prop_value = properties[id_prop_name]
            if id_prop_value:                
                # Handle list values
                if isinstance(id_prop_value, list):
                    for equiv_id in id_prop_value:
                        if equiv_id and str(equiv_id).strip() and equiv_id.lower() not in none_strings:
                            equiv_curie, iri = self.normalize_spoke_identifier(node_type, id_prop_name, str(equiv_id), properties)
                            if equiv_curie and equiv_curie != KNOWN_INVALID:
                                equivalent_ids.add(equiv_curie)
                # Handle string values
                elif isinstance(id_prop_value, str) and id_prop_value.strip() and id_prop_value.lower() not in none_strings:
                    equiv_curie, iri = self.normalize_spoke_identifier(node_type, id_prop_name, id_prop_value, properties)
                    if equiv_curie and equiv_curie != KNOWN_INVALID:
                        equivalent_ids.add(equiv_curie)
        
        # NOTE: Skipping xrefs field for now; quite complicated to determine correct prefix.

        return list(equivalent_ids)

    @staticmethod
    def clean_name(prop_or_source_name: str) -> str:
        return prop_or_source_name.lower().replace(' ', '')

    @staticmethod
    def is_trusted_id_property(property_name: str) -> bool:
        equiv_id_sources = {'chembl', 'drugbank', 'chebi', 'pubchem', 'kegg', 'mesh', 'ensembl', 'omim'}
        exact_fields = {'snomedct', 'icd10', 'icd9', 'accession'}
        prop_name_lower = property_name.lower()
        first_word = prop_name_lower.split('_')[0]
        if prop_name_lower in exact_fields or (first_word in equiv_id_sources and ('id' in prop_name_lower or '_list' in prop_name_lower)):
            return True
        else:
            return False

    def is_known_invalid_id(self, local_id: str, standard_prefix: str | List[str]):
        standard_prefixes = [standard_prefix] if isinstance(standard_prefix, str) else standard_prefix
        if local_id.strip() in self.normalizer.dashes:
            # Some items have an id of just '-', we can ignore these
            return True
        elif 'metacyc.pathway' in standard_prefixes and not local_id.isupper():
            # Meant to catch english names given as identifiers, like Glycan biosynthesis - 2
            return True
        elif 'snomedct' in standard_prefixes and 'e+' in local_id:
            return True
        elif 'dbsnp' in standard_prefixes and ',' in local_id:
            # Detect when multiple dbsnp IDs are concatenated into one ID (skip these for now)
            return True
        elif ('react' in standard_prefixes or 'chebi' in standard_prefixes) and local_id == 'root':
            # We don't want weird abstract 'root' nodes that are present in SPOKE
            return True
        elif 'clo' in standard_prefixes and local_id.startswith('http://'):
            # A couple CLO nodes have old URLs as IDs
            return True
        else:
            return False
