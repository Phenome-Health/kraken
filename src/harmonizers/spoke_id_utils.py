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
        self.prefix_lowercase_map = {prefix.lower(): prefix for prefix in self.biolink_prefixes.keys()}
    

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
            prefix_map['mirbase'] = "https://mirbase.org/hairpin/"  # Biolink has mirbase in here, but their iri doesn't work

            return prefix_map
        else:
            logging.error(f"Failed to download Biolink prefix map ({response.status_code} error). {response.text}")
            sys.exit(1)
    

    def normalize_spoke_identifier(self, node_type: str, source: str, identifier: Any, properties: dict) -> str:
        """
        Simple identifier normalization - handles most cases with minimal complexity
        """
        identifier = str(identifier)
        prefix, local_id = self.get_curie_parts(identifier)
        # If the identifier was a curie, its prefix should override the source
        source = prefix if prefix else source

        source_cleaned = source.lower().replace(' ', '')
        lookup_key = (node_type, source_cleaned)

        if lookup_key == ('Anatomy', 'uberon'):
            if self.is_uberon_id(local_id):
                return self.construct_curie('uberon', local_id)
        elif lookup_key == ('Anatomy', 'mesh_id'):
            if self.is_mesh_id(local_id):
                return self.construct_curie('mesh', local_id)
        elif lookup_key == ('Variant', 'unknown'):
            if self.is_dbsnp_id(local_id):
                return self.construct_curie('dbsnp', local_id)
        elif lookup_key == ('EC', 'explorenz'):
            if self.is_ec_id(local_id):
                return self.construct_curie('ec', local_id)
        elif lookup_key == ('Pathway', 'reactome'):
            if self.is_reactome_id(local_id):
                return self.construct_curie('react', local_id)
        elif lookup_key == ('Pathway', 'wikipathways'):
            local_id_cleaned = local_id.split('_')[0]  # Get rid of version info, like in WP5395_r126912
            if self.is_wikipathways_id(local_id_cleaned):
                return self.construct_curie('wikipathways', local_id_cleaned)



        # if lookup_key not in curie_map:
        logging.error(f"Could not determine prefix for identifier: type: {node_type}, source: {source_cleaned} "
                      f"({source}), identifier: {identifier}, local_id: {local_id}.\n   Properties: {properties}")
        sys.exit(1)


    def get_curie_parts(self, identifier: str) -> Tuple[str, str]:
        num_colons = identifier.count(':')
        if num_colons == 0:
            return '', identifier
        elif num_colons == 1:
            parts = identifier.split(':')
            return parts[0], parts[1]
        else:
            logging.error(f"An identifier has more than one colon in it: {identifier}. Not sure what to do.")
            sys.exit(1)


    def construct_curie(self, prefix_lowercase: str, local_id: str) -> str:
        return f"{self.prefix_lowercase_map[prefix_lowercase]}:{local_id}"


    def _get_local_id(self, identifier: str, source_cleaned: str, properties: dict) -> str:
        if source_cleaned == 'complexportal':
            return str(properties['complex_portal'])  # The 'identifier' SPOKE gives isn't a ComplexPortal ID for some reason; given here instead
        elif source_cleaned == 'celllineontology':
            return identifier.removeprefix('CLO_')
        elif source_cleaned == 'environmentontology':
            return identifier.removeprefix('ENVO_')
        elif 'snomedct' in source_cleaned:
            return str(int(identifier.removeprefix('SNOMED_')))  # They sometimes look like this: "SNOMEDCT": "9209005.0" or "identifier": "SNOMED_1186610007"
        elif source_cleaned == 'countyhealthrankings':
            return identifier.removeprefix('CHR_')
        elif 'householdpulsesurvey' in source_cleaned:
            return identifier.removeprefix('HPS_')
        else:
            return str(identifier)  # Sometimes SPOKE gives ints here


    def _derive_curie_old(self, node_type: str, source_cleaned: str, local_id: str) -> str:
        """Assign prefix based on node source, type and simple heuristics"""
        
        # Source/type-based assignments
        if source_cleaned:  # TODO: if these prefixes/their shortcuts are already in prefix map, can remove the if block? or want the flexible 'if in'?
            if node_type in {'Gene', 'Protein', 'MiRNA'} and ("entrez" in source_cleaned or "ncbi" in source_cleaned or 'mirbase' in source_cleaned):
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
            elif node_type == 'Organism' and 'taxon' in source_cleaned:
                return self._construct_normalized_curie((source_cleaned, node_type), 'ncbitaxon', local_id)
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
            elif source_cleaned == 'accession' and (node_type == 'Gene' or node_type == 'MiRNA') and self.is_mirbase_id(local_id):
                return self._construct_normalized_curie((source_cleaned, node_type), 'mirbase', local_id)
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
                return self._construct_normalized_curie((source_cleaned, node_type), 'cytoband', local_id)  # TODO: Seems like we're getting a lot of CYTOBAND:rs1591517484 (fix)
            elif "ontology" in source_cleaned and node_type == "CellType":
                return self._construct_normalized_curie((source_cleaned, node_type), 'cl', local_id)
            elif "ontology" in source_cleaned and node_type in ["BiologicalProcess", "MolecularFunction", "CellularComponent"]:
                return self._construct_normalized_curie((source_cleaned, node_type), 'go', local_id)
        
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
    def is_mesh_id(local_id: str) -> bool:
        return bool(re.match(r'^D\d+$', local_id)) or bool(re.match(r'^C\d+$', local_id))


    @staticmethod
    def is_metacyc_pathway_id(local_id: str) -> bool:
        return bool(re.match(r'^[A-Z0-9][A-Z0-9-]*$', local_id)) and 'PWY' in local_id


    @staticmethod
    def is_cellosaurus_id(local_id: str) -> bool:
        # Allows: CVCL_[digits or letters]
        return bool(re.match(r'^CVCL_[A-Z0-9]+$', local_id))

    @staticmethod
    def is_mirbase_id(local_id: str) -> bool:
        # Allows: MI[digits] or MIMAT[digits]
        return bool(re.match(r'^MI[0-9]+$', local_id)) or bool(re.match(r'^MIMAT[0-9]+$', local_id))

    @staticmethod
    def is_uberon_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0003233 from UBERON:0003233)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_dbsnp_id(local_id: str) -> bool:
        # Allows: rs followed by digits (e.g., rs1060501038)
        return bool(re.match(r'^rs[0-9]+$', local_id))

    @staticmethod
    def is_ec_id(local_id: str) -> bool:
        # Allows: EC number format (e.g., 3.1.7.2, 1.14.13.M81)
        parts = local_id.split('.')
        if len(parts) < 2 or len(parts) > 4:
            return False

        # Check each part
        for part in parts:
            # Each part must be either:
            # - A number (including 0)
            # - A letter followed by numbers (like M81, B1)
            # - Just a dash (for unspecified sub-subclasses)
            if not re.match(r'^([0-9]+|[A-Z]+[0-9]*|-)$', part):
                return False

        return True

    @staticmethod
    def is_reactome_id(local_id: str) -> bool:
        # Allows: R-HSA-digits (e.g., R-HSA-162582)
        return bool(re.match(r'^R-[A-Z]{3}-[0-9]+$', local_id))

    @staticmethod
    def is_wikipathways_id(local_id: str) -> bool:
        # Allows: WP followed by digits
        return bool(re.match(r'^WP[0-9]+$', local_id))

    @staticmethod
    def is_doid_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0070557 from DOID:0070557)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_ncbigene_id(local_id: str) -> bool:
        # Allows: pure digits (Entrez Gene IDs)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_metacyc_reaction_id(local_id: str) -> bool:
        # Allows: uppercase letters, digits, and hyphens (e.g., R13147, RXN-15029)
        return bool(re.match(r'^R[A-Z0-9-]*$', local_id)) or bool(re.match(r'^RXN-[0-9]+$', local_id))

    @staticmethod
    def is_inchikey_id(local_id: str) -> bool:
        # Allows: standard InChI key format (e.g., AMOFQIUOTAJRKS-UHFFFAOYSA-N)
        return bool(re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', local_id))

    @staticmethod
    def is_go_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0004339 from GO:0004339)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_cl_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0000540 from CL:0000540)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_hpo_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0001234 from HP:0001234)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_uszipcode_id(local_id: str) -> bool:
        # Allows: 5-digit US ZIP codes
        return bool(re.match(r'^[0-9]{5}$', local_id))

    @staticmethod
    def is_ndfrt_id(local_id: str) -> bool:
        # Allows: NDFRT identifiers (typically alphanumeric)
        return bool(re.match(r'^[A-Z0-9_]+$', local_id))

    @staticmethod
    def is_sider_id(local_id: str) -> bool:
        # Allows: SIDER identifiers (typically alphanumeric with possible special chars)
        return bool(re.match(r'^[A-Z0-9._-]+$', local_id))
