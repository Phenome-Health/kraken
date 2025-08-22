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
            if 'KEGG' in prefix_map:
                del prefix_map['KEGG']  # We want to use only KEGG.COMPOUND, KEGG.REACTION, etc.
            return prefix_map
        else:
            logging.error(f"Failed to download Biolink prefix map ({response.status_code} error). {response.text}")
            sys.exit(1)
    
    def _load_prefix_lowercase_map(self) -> Dict[str, str]:
        prefix_lowercase_map = {prefix.lower(): prefix for prefix in self.biolink_prefixes.keys()}
        direct_shortnames = {
            'entrez': 'ncbigene',
            'uniprot': 'uniprotkb',
            'pubchem': 'pubchem.compound',
            'chembl': 'chembl.compound',
            'reactome': 'react'
        }
        for shortname, prefix_lower in direct_shortnames.items():
            prefix_lowercase_map[shortname] = prefix_lowercase_map[prefix_lower]
        return prefix_lowercase_map

    def normalize_prefix_case(self, prefix: str) -> str:
        """Normalize prefix to correct Biolink capitalization"""
        return self.prefix_lowercase_map[prefix.lower()]
    
    def normalize_spoke_identifier(self, node_type: str, source: str, identifier: str) -> str:
        """
        Simple identifier normalization - handles most cases with minimal complexity
        """
        # 1. If it's already a CURIE, use the prefix as our source (prefix overrides source info for determining what vocabulary this is from)
        if ":" in identifier:
            source, identifier = identifier.split(":", 1)
        
        # 2. Construct the normalized curie
        source_lower = source.lower()
        if source_lower in self.prefix_lowercase_map or (source_lower, node_type) in self.prefix_lowercase_map:
            normalized_prefix = self.normalize_prefix_case(source_lower)
            curie = f"{normalized_prefix}:{identifier}"
        else:
            curie = self._derive_curie(node_type, source, identifier)
        return curie
    
    def _derive_curie(self, node_type: str, source: str, identifier: str) -> str:
        """Assign prefix based on node source, type and simple heuristics"""
        
        # Source/type -based assignments
        if source:
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
            elif "ontology" in source_lower and node_type == "CellType":
                return self._construct_normalized_curie(source_lower, 'cl', identifier)
            elif "ontology" in source_lower and node_type in ["BiologicalProcess", "MolecularFunction", "CellularComponent"]:
                return self._construct_normalized_curie(source_lower, 'go', identifier)
        
        # Back up to identifier pattern-based heuristics (simple ones only)
        if identifier.startswith("rs") and identifier[2:].isdigit():
            return f"{self.prefix_lowercase_map['dbsnp']}:{identifier}"
        elif re.match(r'^\d+\.\d+\.\d+\.\d+$', identifier):  # EC number
            return f"{self.prefix_lowercase_map['ec']}:{identifier}"
        elif re.match(r'^\d+$', identifier):  # Pure number
            if node_type == "Gene":
                return f"{self.prefix_lowercase_map['ncbigene']}:{identifier}"
            elif node_type == "Organism":
                return f"{self.prefix_lowercase_map['ncbitaxon']}:{identifier}"
        
        logging.error(f"Could not determine prefix for identifier: type: {node_type}, source: {source}, identifier: {identifier}")
        sys.exit(1)
    
    def extract_equivalent_identifiers(self, node_type: str, properties: Dict) -> List[str]:
        """Extract equivalent IDs from properties - simplified approach"""
        equivalent_ids = set()
        none_strings = {'null', 'none', 'nan'}
        
        # Figure out which properties probably contain an identifier
        relevant_properties = set()
        for property_name in properties.keys():
            property_name_lower = property_name.lower()
            if property_name_lower in self.prefix_lowercase_map:
                relevant_properties.add(property_name)
            else:
                first_word = property_name.split('_')[0].lower()
                if first_word in self.prefix_lowercase_map and ('id' in property_name_lower or '_list' in property_name_lower):
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
