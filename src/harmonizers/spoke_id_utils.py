"""
Simplified identifier normalization for SPOKE - focuses on case normalization and simple heuristics
"""

import json
import logging
import re
from pathlib import Path
import sys
from typing import Any, List, Dict, Set, Tuple, Union, Callable, Optional

import requests
from ..utils.constants import *
from ..utils.general import load_biolink_file
from ..utils.identifiers import IdentifierNorm


class SpokeIDNormalizer:
    """Simple, robust identifier normalization focused on case handling"""
    
    def __init__(self, biolink_version: str):
        self.biolink_version = biolink_version
        self.id_norm = IdentifierNorm(biolink_version)
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
            ('ClinicalLab', 'unknown'): ['loinc', 'umls'],
            ('Complex', 'complexportal'): 'complexportal',
            ('Compound', 'chebi'): 'chebi',
            ('Compound', 'chembl_ids'): 'chembl.compound',
            ('Compound', 'chembl.compound'): 'chembl.compound',
            ('Compound', 'drugbank_ids'): 'drugbank',
            ('Compound', 'inchikey'): 'inchikey',
            ('Compound', 'kegg_compound_ids'): 'kegg.compound',
            ('Compound', 'kegg_drug_ids'): 'kegg.drug',
            ('Compound', 'pubchem_compound_ids'): 'pubchem.compound',
            # ('Compound', 'standardized_smiles'): 'smiles',
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
            ('Haplotype', 'unknown'): ['pharmvar', 'dbsnp'],
            ('Location', 'geonames'): 'geonames',
            ('Location', 'unitedstateszipcode_database'): ['uszipcode', 'fips.place', 'fips.state'],
            ('MiRNA', 'accession'): 'mirbase',
            ('MiRNA', 'mirdb'): 'mirdb',
            ('MolecularFunction', 'go'): 'go',
            ('Organism', 'bv-brc'): 'bvbrc',
            ('Organism', 'ncbi-taxonomy'): 'ncbitaxon',
            ('Pathway', 'reactome'): 'react',
            ('Pathway', 'unknown'): 'metacyc.pathway',
            ('Pathway', 'wikipathways'): 'wikipathways',
            ('PharmacologicClass', 'fdaviadrugcentral'): ['ndfrt', 'mesh'],
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

        source_cleaned = source.lower().replace(' ', '')
        lookup_key = (node_type, source_cleaned)

        if local_id and lookup_key in self.curie_construction_map:
            # Grab the curie construction info
            prefix_entry = self.curie_construction_map[lookup_key]
            chosen_prefix = None

            # Handle when multiple prefixes are listed for this node type/source pair
            if isinstance(prefix_entry, list):
                # We need to figure out which prefix applies based on the local ID's format/structure
                for prefix in prefix_entry:
                    is_valid, _ = self.id_norm.is_valid_id(local_id, prefix)
                    if is_valid in {True, None}:  # It's ok if it's a 'known invalid' format
                        chosen_prefix = prefix
                        break
            else:
                chosen_prefix = prefix_entry

            # Actually construct the curie
            if chosen_prefix:
                curie, iri = self.id_norm.construct_curie(chosen_prefix, local_id)
                if curie:
                    return curie, iri


        logging.error(f"Could not determine proper curie for lookup key: {lookup_key}:\n   type: {node_type}, "
                      f"source: {source_cleaned} ({source}), identifier: {identifier}, local_id: {local_id}"
                      f"\n   Properties: {properties}")
        sys.exit(1)


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
        equiv_id_sources = {'chembl', 'drugbank', 'chebi', 'pubchem', 'kegg', 'mesh', 'ensembl', 'omim'}
        exact_fields = {'snomedct', 'icd10', 'icd9', 'accession'}
        
        # Figure out which properties probably contain an identifier
        relevant_properties = set()
        for property_name in properties.keys():
            prop_name_lower = property_name.lower()
            first_word = prop_name_lower.split('_')[0]
            if prop_name_lower in exact_fields or (first_word in equiv_id_sources and ('id' in prop_name_lower or '_list' in prop_name_lower)):
                relevant_properties.add(property_name)

        # Construct proper curie(s) for each of those properties
        for id_prop_name in relevant_properties:
            id_prop_value = properties[id_prop_name]
            if id_prop_value:                
                # Handle list values
                if isinstance(id_prop_value, list):
                    for equiv_id in id_prop_value:
                        if equiv_id and str(equiv_id).strip() and equiv_id.lower() not in none_strings:
                            equivalent_ids.add(self.normalize_spoke_identifier(node_type, id_prop_name, str(equiv_id), properties)[0])
                # Handle string values
                elif isinstance(id_prop_value, str) and id_prop_value.strip() and id_prop_value.lower() not in none_strings:
                    equivalent_ids.add(self.normalize_spoke_identifier(node_type, id_prop_name, id_prop_value, properties)[0])
        
        # NOTE: Skipping xrefs for now; quite complicated to determine correct prefix. 
        
        return [equiv_id for equiv_id in equivalent_ids if equiv_id]  # Filter out ones that were invalid (returned as empty string)
