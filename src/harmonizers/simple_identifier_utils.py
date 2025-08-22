"""
Simplified identifier normalization for SPOKE - focuses on case normalization and simple heuristics
"""

import json
import logging
import re
from pathlib import Path
import sys
from typing import Any, List, Dict, Set, Tuple, Union

import requests


class SimpleIdentifierNormalizer:
    """Simple, robust identifier normalization focused on case handling"""
    
    def __init__(self, biolink_version: str):
        self.biolink_version = biolink_version
        self.biolink_prefixes = self._load_biolink_prefixes()
        self.prefix_lowercase_map = self._load_prefix_lowercase_map()
    
    def _load_biolink_prefixes(self) -> Dict[str, str]:
        """Load Biolink model prefix map"""
        logging.info(f"Grabbing biolink prefix map for version: {self.biolink_version}")
        url = f"https://raw.githubusercontent.com/biolink/biolink-model/refs/tags/v{self.biolink_version}/project/prefixmap/biolink-model-prefix-map.json"
        response = requests.get(url)
        if response.ok:
            prefix_map = response.json()

            # Remove prefixes as needed
            if 'KEGG' in prefix_map:
                del prefix_map['KEGG']  # We want to use only KEGG.COMPOUND, KEGG.REACTION, etc.
            # Add prefixes as needed (ones we're making up, that don't exist in biolink)
            prefix_map['USZIPCODE'] = "https://www.unitedstateszipcodes.org/"
            prefix_map['SMILES'] = "https://pubchem.ncbi.nlm.nih.gov/compound/"
            prefix_map['CELLOSAURUS'] = "https://web.expasy.org/cellosaurus/"
            prefix_map['VESICLEPEDIA'] = "http://microvesicles.org/exp_summary?exp_id="
            prefix_map['NDFRT'] = "http://purl.bioontology.org/ontology/NDFRT/"

            return prefix_map
        else:
            logging.error(f"Failed to download Biolink prefix map ({response.status_code} error). {response.text}")
            sys.exit(1)
    
    def _load_prefix_lowercase_map(self) -> Dict[str, str]:
        prefix_lowercase_map = {prefix.lower(): prefix for prefix in self.biolink_prefixes.keys()}
        direct_shortnames = {  # Need to preload those relevant to the equiv_id extraction (otherwise not guaranteed to be there)
            'entrez': 'ncbigene',
            'uniprot': 'uniprotkb',
            'pubchem': 'pubchem.compound',
            'chembl': 'chembl.compound',
            'reactome': 'react'
        }
        for shortname, prefix_lower in direct_shortnames.items():
            prefix_lowercase_map[shortname] = prefix_lowercase_map[prefix_lower]
        return prefix_lowercase_map
    
    def normalize_spoke_identifier(self, node_type: str, source: str, identifier: Any) -> str:
        """
        Simple identifier normalization - handles most cases with minimal complexity
        """
        identifier = str(identifier)  # Sometimes SPOKE gives ints here

        # 1. If it's already a CURIE, use the prefix as our source (prefix overrides source info for determining what vocabulary this is from)
        if ":" in identifier:
            source, identifier = identifier.split(":", 1)
        
        # 2. Construct the normalized curie
        source_lower = source.lower()
        if source_lower in self.prefix_lowercase_map:
            normalized_prefix = self.prefix_lowercase_map[source_lower]
            curie = f"{normalized_prefix}:{identifier}"
        elif (source_lower, node_type) in self.prefix_lowercase_map:
            normalized_prefix = self.prefix_lowercase_map[(source_lower, node_type)]
            curie = f"{normalized_prefix}:{identifier}"
        else:
            curie = self._derive_curie(node_type, source, identifier)
        
        return curie
    
    def _derive_curie(self, node_type: str, source: str, identifier: str) -> str:
        """Assign prefix based on node source, type and simple heuristics"""
        
        # Source/type-based assignments
        if source:  # TODO: if these prefixes/their shortcuts are already in prefix map, can remove the if block? or want the flexible 'if in'?
            source_lower = source.lower()
            if "entrez" in source_lower or "ncbi" in source_lower:
                return self._construct_normalized_curie(source_lower, 'ncbigene', identifier)
            elif "uniprot" in source_lower:
                return self._construct_normalized_curie(source_lower, 'uniprotkb', identifier)
            elif "dbsnp" in source_lower:
                return self._construct_normalized_curie(source_lower, 'dbsnp', identifier)
            elif "mesh" in source_lower:
                return self._construct_normalized_curie(source_lower, 'mesh', identifier)
            elif "omim" in source_lower:
                return self._construct_normalized_curie(source_lower, 'omim', identifier)
            elif "chebi" in source_lower:
                return self._construct_normalized_curie(source_lower, 'chebi', identifier)
            elif "drugbank" in source_lower:
                return self._construct_normalized_curie(source_lower, 'drugbank', identifier)
            elif "pubchem" in source_lower:
                return self._construct_normalized_curie(source_lower, 'pubchem.compound', identifier)
            elif "chembl" in source_lower:
                return self._construct_normalized_curie(source_lower, 'chembl.compound', identifier)
            elif "uberon" in source_lower:
                return self._construct_normalized_curie(source_lower, 'uberon', identifier)
            elif "metacyc" in source_lower and node_type == 'Reaction':
                return self._construct_normalized_curie((source_lower, node_type), 'metacyc.reaction', identifier)
            elif node_type == 'EC':
                return self._construct_normalized_curie((source_lower, node_type), 'ec', identifier)
            elif node_type == 'ClinicalLab':
                return self._construct_normalized_curie((source_lower, node_type), 'loinc', identifier)
            elif "disease" in source_lower and "ontology" in source_lower:
                return self._construct_normalized_curie(source_lower, 'doid', identifier)
            elif "pfam" in source_lower:
                return self._construct_normalized_curie(source_lower, 'pfam', identifier)
            elif "taxonomy" in source_lower:
                return self._construct_normalized_curie(source_lower, 'ncbitaxon', identifier)
            elif 'kegg' in source_lower:
                if "drug" in source_lower:
                    return self._construct_normalized_curie(source_lower, 'kegg.drug', identifier)
                elif node_type == "Compound":
                    return self._construct_normalized_curie((source_lower, node_type), 'kegg.compound', identifier)
                elif node_type == "Reaction":
                    return self._construct_normalized_curie((source_lower, node_type), 'kegg.reaction', identifier)
            elif "ensembl" in source_lower:
                return self._construct_normalized_curie(source_lower, 'ensembl', identifier)
            elif "icd10" in source_lower:
                return self._construct_normalized_curie(source_lower, 'icd10', identifier)
            elif 'snomedct' in source_lower:
                return self._construct_normalized_curie(source_lower, 'snomedct', identifier)
            elif 'react' in source_lower or 'reactome' in source_lower:
                return self._construct_normalized_curie(source_lower, 'react', identifier)
            elif source_lower == 'unitedstateszipcode_database':
                return self._construct_normalized_curie(source_lower, 'uszipcode', identifier)
            elif source_lower == 'fda via drugcentral':
                return self._construct_normalized_curie(source_lower, 'ndfrt', identifier)
            elif 'smiles' in source_lower:
                return self._construct_normalized_curie(source_lower, 'smiles', identifier)
            elif 'vesiclepedia' in source_lower:
                return self._construct_normalized_curie(source_lower, 'vesiclepedia', identifier)
            elif "ontology" in source_lower and node_type == "CellType":
                return self._construct_normalized_curie(source_lower, 'cl', identifier)
            elif "ontology" in source_lower and node_type in ["BiologicalProcess", "MolecularFunction", "CellularComponent"]:
                return self._construct_normalized_curie(source_lower, 'go', identifier)
        
        # Back up to identifier pattern-based heuristics (simple ones only)
        types_mesh_used_for = {'Symptom', 'SideEffect'}
        if identifier.startswith("rs") and identifier[2:].isdigit():
            return f"{self.prefix_lowercase_map['dbsnp']}:{identifier}"
        elif node_type in types_mesh_used_for and self._is_mesh_id(identifier):
            return f"{self.prefix_lowercase_map['mesh']}:{identifier}"
        elif node_type == 'Pathway' and source_lower == 'unknown' and self._is_metacyc_pathway_id(identifier):
            return f"{self.prefix_lowercase_map['metacyc.reaction']}:{identifier}"
        elif node_type == 'EC' and re.match(r'^\d+\.\d+\.\d+\.\d+$', identifier):  # EC number
            return f"{self.prefix_lowercase_map['ec']}:{identifier}"
        elif node_type == 'CellLine' and self.is_cellosaurus_id(identifier):
            return f"{self.prefix_lowercase_map['cellosaurus']}:{identifier}"
        elif re.match(r'^\d+$', identifier):  # Pure number
            if node_type == "Gene":
                return f"{self.prefix_lowercase_map['ncbigene']}:{identifier}"
            elif node_type == "Organism":
                return f"{self.prefix_lowercase_map['ncbitaxon']}:{identifier}"
        elif node_type == 'Pathway' and source_lower == 'unknown':
            return f"SPOKE:{re.sub(r'\s+', '_', identifier.strip())}"  # For identifiers like 'Glycan biosynthesis - 2' 
        
        logging.error(f"Could not determine prefix for identifier: type: {node_type}, source: {source}, identifier: {identifier}")
        sys.exit(1)
    
    def extract_equivalent_identifiers(self, node_type: str, properties: Dict) -> List[str]:
        """Extract equivalent IDs from properties - simplified approach"""
        equivalent_ids = set()
        none_strings = {'null', 'none', 'nan'}
        equiv_id_sources = {'chembl', 'drugbank', 'chebi', 'pubchem', 'kegg', 'mesh', 'ensembl', 'omim', 'icd10', 'snomedct'}
        exact_fields = {'standardized_smiles'}
        
        # Figure out which properties probably contain an identifier
        relevant_properties = set()
        for property_name in properties.keys():
            property_name_lower = property_name.lower()
            first_word = property_name.split('_')[0].lower()
            if property_name in exact_fields or (first_word in equiv_id_sources and ('id' in property_name_lower or '_list' in property_name_lower)):
                relevant_properties.add(property_name)
        
        # Construct proper curie(s) for each of those properties
        for id_prop_name in relevant_properties:
            id_prop_value = properties[id_prop_name]
            if id_prop_value:                
                # Handle list values
                if isinstance(id_prop_value, list):
                    for item in id_prop_value:
                        if item and str(item).strip() and item.lower() not in none_strings:
                            equivalent_ids.add(self.normalize_spoke_identifier(node_type, id_prop_name, str(item)))
                # Handle string values
                elif isinstance(id_prop_value, str) and id_prop_value.strip() and id_prop_value.lower() not in none_strings:
                    equivalent_ids.add(self.normalize_spoke_identifier(node_type, id_prop_name, id_prop_value))
        
        # NOTE: Skipping xrefs for now; quite complicated to determine correct prefix. 
        
        return list(equivalent_ids)
    
    def _construct_normalized_curie(self, prefix_key: Union[str, Tuple[str, str]], prefix_lowercase: str, identifier: str) -> str:
        normalized_prefix = self.prefix_lowercase_map[prefix_lowercase]
        # Cache this mapping (source --> normalized prefix), for faster processing later
        self.prefix_lowercase_map[prefix_key] = normalized_prefix
        return f"{normalized_prefix}:{identifier}"
    
    @staticmethod
    def _is_mesh_id(local_id: str) -> bool:
        return bool(re.match(r'^D\d+$', local_id)) or bool(re.match(r'^C\d+$', local_id))

    @staticmethod
    def _is_metacyc_pathway_id(local_id: str) -> bool:
        return bool(re.match(r'^[A-Z0-9][A-Z0-9-]*$', local_id)) and 'PWY' in local_id

    @staticmethod
    def is_cellosaurus_id(local_id: str) -> bool:
        return bool(re.match(r'^CVCL_\d+$', local_id))
