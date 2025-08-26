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


class SpokeIDNormalizer:
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
            prefix_map['BVBRC'] = "https://www.bv-brc.org/view/Genome/"
            prefix_map['GeoNames'] = "http://www.geonames.org/search.html?q="  # Note: this doesn't go exactly to page for item, but closest I could find
            prefix_map['NHANES'] = "https://dsld.od.nih.gov/label/"  # These IRIs work, but weirdly SPOKE's identifiers for these nodes don't match what they have..
            prefix_map['MIRDB'] = "https://mirdb.org/cgi-bin/mature_mir.cgi?name="
            prefix_map['CYTOBAND'] = ""  # Haven't found good iri for these yet..
            prefix_map['CHR'] = ""  # Haven't found good iri for these yet..
            prefix_map['AHRQ'] = ""
            prefix_map['HPS'] = ""  # Household Pulse Survey

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
    

    def normalize_spoke_identifier(self, node_type: str, source: str, identifier: Any, properties: dict) -> str:
        """
        Simple identifier normalization - handles most cases with minimal complexity
        """
        identifier = str(identifier)

        # 1. Extract the properly-formatted local ID
        if ':' in identifier and not identifier.startswith('http'):  # Sometimes SPOKE gives full URL as node 'identifier'
            # If it's already a CURIE, use the prefix as our source (prefix overrides source info for determining what vocabulary this is from)
            source, local_id = identifier.split(':', 1)
            source_cleaned = source.lower().replace(' ', '')
        else:
            # Extract the local ID for the curie, reformatting it as necessary
            source_cleaned = source.lower().replace(' ', '')
            local_id = self._get_local_id(identifier, source_cleaned, properties)
        local_id = local_id.strip()
        
        # 3. Construct the normalized curie (using cached prefix mapping if available)
        if source_cleaned in self.prefix_lowercase_map:
            normalized_prefix = self.prefix_lowercase_map[source_cleaned]
            curie = f"{normalized_prefix}:{local_id}"
        elif (source_cleaned, node_type) in self.prefix_lowercase_map:
            normalized_prefix = self.prefix_lowercase_map[(source_cleaned, node_type)]
            curie = f"{normalized_prefix}:{local_id}"
        else:
            curie = self._derive_curie(node_type, source_cleaned, local_id)
        
        return curie
    

    def _get_local_id(self, identifier: str, source_cleaned: str, properties: dict) -> str:
        if source_cleaned == 'complexportal':
            return str(properties['complex_portal'])  # The 'identifier' SPOKE gives isn't a ComplexPortal ID for some reason; given here instead
        elif source_cleaned == 'celllineontology':
            return identifier.removeprefix('CLO_')
        elif source_cleaned == 'environmentontology':
            return identifier.removeprefix('ENVO_')
        elif 'snomedct' in source_cleaned:
            return identifier.removeprefix('SNOMED_')
        elif source_cleaned == 'countyhealthrankings':
            return identifier.removeprefix('CHR_')
        elif 'householdpulsesurvey' in source_cleaned:
            return identifier.removeprefix('HPS_')
        else:
            return str(identifier)  # Sometimes SPOKE gives ints here


    def _derive_curie(self, node_type: str, source_cleaned: str, local_id: str) -> str:
        """Assign prefix based on node source, type and simple heuristics"""
        
        # Source/type-based assignments
        if source_cleaned:  # TODO: if these prefixes/their shortcuts are already in prefix map, can remove the if block? or want the flexible 'if in'?
            if "entrez" in source_cleaned or "ncbi" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'ncbigene', local_id)
            elif "uniprot" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'uniprotkb', local_id)
            elif "dbsnp" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'dbsnp', local_id)
            elif "mesh" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'mesh', local_id)
            elif "omim" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'omim', local_id)
            elif "chebi" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'chebi', local_id)
            elif "drugbank" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'drugbank', local_id)
            elif "pubchem" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'pubchem.compound', local_id)
            elif "chembl" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'chembl.compound', local_id)
            elif "uberon" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'uberon', local_id)
            elif "metacyc" in source_cleaned and node_type == 'Reaction':
                return self._construct_normalized_curie((source_cleaned, node_type), 'metacyc.reaction', local_id)
            elif node_type == 'EC':
                return self._construct_normalized_curie((source_cleaned, node_type), 'ec', local_id)
            elif node_type == 'ClinicalLab':
                return self._construct_normalized_curie((source_cleaned, node_type), 'loinc', local_id)
            elif "disease" in source_cleaned and "ontology" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'doid', local_id)
            elif "pfam" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'pfam', local_id)
            elif "taxonomy" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'ncbitaxon', local_id)
            elif 'kegg' in source_cleaned:
                if "drug" in source_cleaned:
                    return self._construct_normalized_curie(source_cleaned, 'kegg.drug', local_id)
                elif node_type == "Compound":
                    return self._construct_normalized_curie((source_cleaned, node_type), 'kegg.compound', local_id)
                elif node_type == "Reaction":
                    return self._construct_normalized_curie((source_cleaned, node_type), 'kegg.reaction', local_id)
            elif "ensembl" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'ensembl', local_id)
            elif "icd10" in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'icd10', local_id)
            elif 'snomedct' in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'snomedct', local_id)
            elif 'react' in source_cleaned or 'reactome' in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'react', local_id)
            elif source_cleaned == 'unitedstateszipcode_database':
                return self._construct_normalized_curie(source_cleaned, 'uszipcode', local_id)
            elif source_cleaned == 'fdaviadrugcentral':
                return self._construct_normalized_curie(source_cleaned, 'ndfrt', local_id)
            elif source_cleaned == 'bv-brc':
                return self._construct_normalized_curie(source_cleaned, 'bvbrc', local_id)
            elif source_cleaned == 'nhanes':
                return self._construct_normalized_curie(source_cleaned, 'nhanes', local_id)
            elif source_cleaned == 'mirdb':
                return self._construct_normalized_curie(source_cleaned, 'mirdb', local_id)
            elif source_cleaned == 'celllineontology':
                return self._construct_normalized_curie(source_cleaned, 'clo', local_id)
            elif source_cleaned == 'environmentontology':
                return self._construct_normalized_curie(source_cleaned, 'envo', local_id)
            elif 'geonames' in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'geonames', local_id)
            elif 'smiles' in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'smiles', local_id)
            elif 'vesiclepedia' in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'vesiclepedia', local_id)
            elif 'complexportal' in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'complexportal', local_id)
            elif source_cleaned == 'countyhealthrankings':
                return self._construct_normalized_curie(source_cleaned, 'chr', local_id)
            elif source_cleaned == 'ahrqsdohdatabase':
                return self._construct_normalized_curie(source_cleaned, 'ahrq', local_id)
            elif 'householdpulsesurvey' in source_cleaned:
                return self._construct_normalized_curie(source_cleaned, 'hps', local_id)
            elif node_type == 'Cytoband' and source_cleaned == 'unknown':
                return self._construct_normalized_curie(source_cleaned, 'cytoband', local_id)  # TODO: Seems like we're getting a lot of CYTOBAND:rs1591517484 (fix)
            elif "ontology" in source_cleaned and node_type == "CellType":
                return self._construct_normalized_curie(source_cleaned, 'cl', local_id)
            elif "ontology" in source_cleaned and node_type in ["BiologicalProcess", "MolecularFunction", "CellularComponent"]:
                return self._construct_normalized_curie(source_cleaned, 'go', local_id)
        
        # Back up to identifier pattern-based heuristics (simple ones only)
        types_mesh_used_for = {'Symptom', 'SideEffect'}
        if local_id.startswith("rs") and local_id[2:].isdigit():
            return f"{self.prefix_lowercase_map['dbsnp']}:{local_id}"
        elif node_type in types_mesh_used_for and self._is_mesh_id(local_id):
            return f"{self.prefix_lowercase_map['mesh']}:{local_id}"
        elif node_type == 'Pathway' and source_cleaned == 'unknown' and self._is_metacyc_pathway_id(local_id):
            return f"{self.prefix_lowercase_map['metacyc.reaction']}:{local_id}"
        elif node_type == 'EC' and re.match(r'^\d+\.\d+\.\d+\.\d+$', local_id):  # EC number
            return f"{self.prefix_lowercase_map['ec']}:{local_id}"
        elif node_type == 'CellLine' and self.is_cellosaurus_id(local_id):
            return f"{self.prefix_lowercase_map['cellosaurus']}:{local_id}"
        elif re.match(r'^\d+$', local_id):  # Pure number
            if node_type == "Gene":
                return f"{self.prefix_lowercase_map['ncbigene']}:{local_id}"
            elif node_type == "Organism":
                return f"{self.prefix_lowercase_map['ncbitaxon']}:{local_id}"
        elif node_type == 'Pathway' and source_cleaned == 'unknown':
            return f"SPOKE:{re.sub(r'\s+', '_', local_id.strip())}"  # For identifiers like 'Glycan biosynthesis - 2' 
        
        logging.error(f"Could not determine prefix for identifier: type: {node_type}, source: {source_cleaned}, identifier: {local_id}")
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
                    for equiv_id in id_prop_value:
                        if equiv_id and str(equiv_id).strip() and equiv_id.lower() not in none_strings:
                            equivalent_ids.add(self.normalize_spoke_identifier(node_type, id_prop_name, str(equiv_id), properties))
                # Handle string values
                elif isinstance(id_prop_value, str) and id_prop_value.strip() and id_prop_value.lower() not in none_strings:
                    equivalent_ids.add(self.normalize_spoke_identifier(node_type, id_prop_name, id_prop_value, properties))
        
        # NOTE: Skipping xrefs for now; quite complicated to determine correct prefix. 
        
        return list(equivalent_ids)
    

    def _construct_normalized_curie(self, prefix_key: Union[str, Tuple[str, str]], prefix_lowercase: str, local_id: str) -> str:
        normalized_prefix = self.prefix_lowercase_map[prefix_lowercase]
        # Cache this mapping (source --> normalized prefix), for faster processing later
        self.prefix_lowercase_map[prefix_key] = normalized_prefix
        return f"{normalized_prefix}:{local_id}"
    

    @staticmethod
    def _is_mesh_id(local_id: str) -> bool:
        return bool(re.match(r'^D\d+$', local_id)) or bool(re.match(r'^C\d+$', local_id))


    @staticmethod
    def _is_metacyc_pathway_id(local_id: str) -> bool:
        return bool(re.match(r'^[A-Z0-9][A-Z0-9-]*$', local_id)) and 'PWY' in local_id


    @staticmethod
    def is_cellosaurus_id(local_id: str) -> bool:
        # Allows: CVCL_[digits or letters]
        return bool(re.match(r'^CVCL_[A-Z0-9]+$', local_id))
