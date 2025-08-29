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


class SpokeIDNormalizer:
    """Simple, robust identifier normalization focused on case handling"""
    
    def __init__(self, biolink_version: str):
        self.biolink_version = biolink_version
        self.normalized_prefixes_to_iris = self._load_prefix_to_iri_map()
        self.prefix_lowercase_map = {prefix.lower(): prefix for prefix in self.normalized_prefixes_to_iris.keys()}
        self.validator_map = self._load_validator_map()
        self.curie_construction_map, self.prefix_prop, self.cleaner_prop = self._load_curie_construction_map()
        self.known_invalid = "KNOWN_INVALID"
    

    def _load_prefix_to_iri_map(self) -> Dict[str, str]:
        """Load Biolink model prefix map"""
        logging.info(f"Grabbing biolink prefix map for version: {self.biolink_version}")
        url = f"https://raw.githubusercontent.com/biolink/biolink-model/refs/tags/v{self.biolink_version}/project/prefixmap/biolink-model-prefix-map.json"
        response = requests.get(url)
        if response.ok:
            prefix_to_iri_map = response.json()

            # Remove prefixes as needed
            if 'KEGG' in prefix_to_iri_map:
                del prefix_to_iri_map['KEGG']  # We want to use only KEGG.COMPOUND, KEGG.REACTION, etc.

            # Add prefixes as needed (ones we're making up, that don't exist in biolink)
            prefix_to_iri_map['USZIPCODE'] = "https://www.unitedstateszipcodes.org/"
            prefix_to_iri_map['SMILES'] = "https://pubchem.ncbi.nlm.nih.gov/compound/"
            prefix_to_iri_map['CELLOSAURUS'] = "https://web.expasy.org/cellosaurus/"
            prefix_to_iri_map['VESICLEPEDIA'] = "http://microvesicles.org/exp_summary?exp_id="
            prefix_to_iri_map['NDFRT'] = "http://purl.bioontology.org/ontology/NDFRT/"
            prefix_to_iri_map['BVBRC'] = "https://www.bv-brc.org/view/Genome/"
            prefix_to_iri_map['GeoNames'] = "http://www.geonames.org/search.html?q="  # Note: this doesn't go exactly to page for item, but closest I could find
            prefix_to_iri_map['NHANES'] = "https://dsld.od.nih.gov/label/"  # These IRIs work, but weirdly SPOKE's identifiers for these nodes don't match what they have..
            prefix_to_iri_map['MIRDB'] = "https://mirdb.org/cgi-bin/mature_mir.cgi?name="
            prefix_to_iri_map['CYTOBAND'] = ""  # Haven't found good iri for these yet..
            prefix_to_iri_map['CHR'] = ""  # Haven't found good iri for these yet..
            prefix_to_iri_map['AHRQ'] = ""
            prefix_to_iri_map['HPS'] = ""  # Household Pulse Survey
            prefix_to_iri_map['mirbase'] = "https://mirbase.org/hairpin/"  # Biolink has mirbase in here, but their iri doesn't work
            prefix_to_iri_map['metacyc.pathway'] = "https://metacyc.org/pathway?orgid=META&id="  # Biolink has metacyc.reaction, but not pathway

            # Override prefixes as needed (if Biolink's iri is broken)
            prefix_to_iri_map['omim'] = "https://omim.org/entry/"

            return prefix_to_iri_map
        else:
            logging.error(f"Failed to download Biolink prefix map ({response.status_code} error). {response.text}")
            sys.exit(1)


    def _load_validator_map(self) -> Dict[str, Callable]:
        return {
            'chebi': self.is_chebi_id,
            'chembl.compound': self.is_chembl_id,
            'cl': self.is_cl_id,
            'dbsnp': self.is_dbsnp_id,
            'doid': self.is_doid_id,
            'ec': self.is_ec_id,
            'icd9': self.is_icd9_id,
            'icd10': self.is_icd10_id,
            'inchikey': self.is_inchikey_id,
            'kegg.reaction': self.is_kegg_reaction_id,
            'loinc': self.is_loinc_id,
            'mesh': self.is_mesh_id,
            'metacyc.pathway': self.is_metacyc_pathway_id,
            'ncbigene': self.is_ncbigene_id,
            'ncbitaxon': self.is_ncbitaxon_id,
            'omim': self.is_omim_id,
            'pfam': self.is_pfam_id,
            'pubchem.compound': self.is_pubchem_compound_id,
            'react': self.is_reactome_id,
            'smiles': self.is_smiles_string,
            'snomedct': self.is_snomedct_id,
            'uberon': self.is_uberon_id,
            'umls': self.is_umls_id,
            'wikipathways': self.is_wikipathways_id,
        }

    def _load_curie_construction_map(self) -> Tuple[Dict[Tuple[str, str], Dict[str, Union[str, List[str], Callable]]], str, str]:
        prefix = 'prefix'
        cleaner = 'local_id_cleaner'
        curie_construction_map = {
            ('Anatomy', 'mesh_id'): {prefix: 'mesh'},
            ('Anatomy', 'uberon'): {prefix: 'uberon'},
            ('CellType', 'cl'): {prefix: 'cl'},
            ('ClinicalLab', 'unknown'): {prefix: ['loinc', 'umls']},
            ('Compound', 'chebi'): {prefix: 'chebi'},
            ('Compound', 'chembl_ids'): {prefix: 'chembl.compound'},
            ('Compound', 'chembl.compound'): {prefix: 'chembl.compound'},
            ('Compound', 'inchikey'): {prefix: 'inchikey'},
            ('Compound', 'pubchem_compound_ids'): {prefix: 'pubchem.compound'},
            ('Compound', 'standardized_smiles'): {prefix: 'smiles'},
            ('Disease', 'doid'): {prefix: 'doid'},
            ('Disease', 'mesh_list'): {prefix: 'mesh'},
            ('Disease', 'omim_list'): {prefix: 'omim'},
            ('EC', 'explorenz'): {prefix: 'ec'},
            ('Gene', 'entrezgene'): {prefix: 'ncbigene'},
            ('Organism', 'ncbi-taxonomy'): {prefix: 'ncbitaxon'},
            ('Pathway', 'reactome'): {prefix: 'react'},
            ('Pathway', 'unknown'): {prefix: 'metacyc.pathway'},
            ('Pathway', 'wikipathways'): {prefix: 'wikipathways', cleaner: self.clean_wikipathways_id},
            ('ProteinDomain', 'pfam'): {prefix: 'pfam'},
            ('Reaction', 'kegg'): {prefix: 'kegg.reaction'},
            ('SideEffect', 'sider4.1'): {prefix: 'umls'},  # They give CUIs for nodes with SIDER source
            ('Symptom', 'hpo'): {prefix: 'mesh'},  # They give MeSH IDs for nodes with HPO source
            ('Symptom', 'icd9'): {prefix: 'icd9'},
            ('Symptom', 'icd10'): {prefix: 'icd10'},
            ('Symptom', 'mesh'): {prefix: 'mesh'},
            ('Symptom', 'snomedct'): {prefix: 'snomedct', cleaner: self.convert_float_to_int_str},
            ('Variant', 'unknown'): {prefix: 'dbsnp'},
        }
        return curie_construction_map, prefix, cleaner


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
            curie_info = self.curie_construction_map[lookup_key]
            if isinstance(curie_info[self.prefix_prop], list):
                # We need to figure out which prefix applies based on the local ID's format/structure
                prefixes_lowercase = [prefix for prefix in curie_info[self.prefix_prop]
                                      if self.validator_map[prefix](local_id)]
                prefix_lowercase = prefixes_lowercase[0] if prefixes_lowercase else None
            else:
                prefix_lowercase = curie_info[self.prefix_prop]

            # Clean up the local ID, if needed
            if curie_info.get(self.cleaner_prop):
                local_id = curie_info[self.cleaner_prop](local_id)

            if prefix_lowercase:
                curie, iri = self.construct_curie(prefix_lowercase, local_id)
                if curie == self.known_invalid:
                    return "", ""
                elif curie:
                    print(f"curie: {curie}, name: {properties.get('name')}, iri: {iri}")
                    return curie, iri


        logging.error(f"Could not determine proper curie for lookup key: {lookup_key}:\n   type: {node_type}, "
                      f"source: {source_cleaned} ({source}), identifier: {identifier}, local_id: {local_id}"
                      f"\n   Properties: {properties}")
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


    def construct_curie(self, prefix_lowercase: str, local_id: str) -> Tuple[str, str]:
        is_valid_id_validator = self.validator_map[prefix_lowercase]
        is_valid_id = is_valid_id_validator(local_id)
        if is_valid_id:
            prefix_normalized = self.prefix_lowercase_map[prefix_lowercase]
            iri_root = self.normalized_prefixes_to_iris[prefix_normalized]
            iri = f"{iri_root}{local_id}" if iri_root else ""
            return f"{prefix_normalized}:{local_id}", iri
        elif is_valid_id is None:
            # Indicates this is a known invalid ID format for this node_type, source pair; we will skip it
            logging.warning(f"Local id {local_id} is invalid for {prefix_lowercase} (known invalid format)")
            return self.known_invalid, ''
        else:
            # This is an unknown invalid ID format; we want to return nothing, which will halt processing after logging
            logging.error(f"Local id {local_id} is invalid for {prefix_lowercase} (UNKNOWN invalid format)")
            return '', ''
    

    def extract_equivalent_identifiers(self, node_type: str, properties: Dict) -> List[str]:
        """Extract equivalent IDs from properties - simplified approach"""
        equivalent_ids = set()
        none_strings = {'null', 'none', 'nan'}
        equiv_id_sources = {'chembl', 'drugbank', 'chebi', 'pubchem', 'kegg', 'mesh', 'ensembl', 'omim'}
        exact_fields = {'standardized_smiles', 'snomedct', 'icd10', 'icd9'}
        
        # Figure out which properties probably contain an identifier
        relevant_properties = set()
        for property_name in properties.keys():
            prop_name_lower = property_name.lower()
            first_word = prop_name_lower.split('_')[0]
            if prop_name_lower in exact_fields or (first_word in equiv_id_sources and ('id' in prop_name_lower or '_list' in prop_name_lower)):
                relevant_properties.add(property_name)

        print(f"     keys are: {properties.keys()}")
        print(f"    relevant properties for equiv ids are: {relevant_properties}\n")
        # if properties.get('SNOMEDCT') == '1.62248710001191e+16':
        #     raise ValueError('At the node I want')

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
    

    def _construct_normalized_curie(self, prefix_key: Union[str, Tuple[str, str]], prefix_lowercase: str, local_id: str) -> str:
        normalized_prefix = self.prefix_lowercase_map[prefix_lowercase]
        # Cache this mapping (source --> normalized prefix), for faster processing later
        self.prefix_lowercase_map[prefix_key] = normalized_prefix
        return f"{normalized_prefix}:{local_id}"

    @staticmethod
    def is_loinc_id(local_id: str) -> bool:
        # LOINC codes: digits followed by dash and check digit (e.g., 27858-0)
        # or LP codes: LP followed by digits and dash-digit (e.g., LP32606-3)
        return bool(re.match(r'^(LP)?\d+-\d$', local_id))

    @staticmethod
    def is_mesh_id(local_id: str) -> bool:
        return bool(re.match(r'^D\d+$', local_id)) or bool(re.match(r'^C\d+$', local_id))

    @staticmethod
    def is_metacyc_pathway_id(local_id: str) -> bool:
        # MetaCyc pathway IDs: either PWY-#### or DESCRIPTIVE-NAME-PWY (but not both)
        return bool(re.match(r'^PWY-\d+$', local_id)) or bool(re.match(r'^[A-Z0-9\-]+-PWY$', local_id))

    @staticmethod
    def is_snomedct_id(local_id: str) -> Optional[bool]:
        # SNOMED CT IDs are numeric strings
        if 'e+' in local_id:  # Known spoke bug where some snomed ct IDs are in scientific notation, like '1.62248710001191e+16'
            return None
        else:
            return local_id.isdigit()

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
        if len(parts) < 1 or len(parts) > 4:
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
    def is_kegg_reaction_id(local_id: str) -> bool:
        # Allows: R followed by exactly 5 digits
        return bool(re.match(r'^R\d{5}$', local_id))

    @staticmethod
    def is_pubchem_compound_id(local_id: str) -> bool:
        # Allows: positive integers (PubChem CIDs are numeric)
        return local_id.isdigit() and int(local_id) > 0

    @staticmethod
    def is_smiles_string(local_id: str) -> bool:
        # Basic SMILES validation: multiple element symbols and valid characters
        uppercase_count = sum(1 for c in local_id if c.isupper())
        valid_chars = set('BCNOPSFHIKLMWUVYXZbcnopslr[]()=#+\\/@.-0123456789')
        has_valid_chars = all(c in valid_chars for c in local_id)
        return uppercase_count >= 2 and has_valid_chars

    @staticmethod
    def is_wikipathways_id(local_id: str) -> bool:
        # Allows: WP followed by digits
        return bool(re.match(r'^WP[0-9]+$', local_id))

    @staticmethod
    def clean_wikipathways_id(local_id: str) -> str:
        # Get rid of version suffix info, like in WP5395_r126912
        return local_id.split('_')[0]

    @staticmethod
    def is_doid_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0070557 from DOID:0070557)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_ncbigene_id(local_id: str) -> bool:
        # Allows: pure digits (Entrez Gene IDs)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_ncbitaxon_id(local_id: str) -> bool:
        # NCBI Taxonomy IDs are positive integers
        return local_id.isdigit() and int(local_id) > 0

    @staticmethod
    def is_omim_id(local_id: str) -> bool:
        # OMIM IDs are 6-digit numbers
        return local_id.isdigit() and len(local_id) == 6

    @staticmethod
    def is_pfam_id(local_id: str) -> bool:
        # Pfam IDs: PF followed by digits
        return bool(re.match(r'^PF\d+$', local_id))

    @staticmethod
    def is_umls_cui(local_id: str) -> bool:
        # UMLS CUI: C followed by 7 digits
        return bool(re.match(r'^C\d{7}$', local_id))

    @staticmethod
    def is_umls_mthu_id(local_id: str) -> bool:
        # UMLS MTHU identifiers: MTHU followed by 6 digits
        return bool(re.match(r'^MTHU\d{6}$', local_id))

    def is_umls_id(self, local_id: str) -> bool:
        # UMLS CUI or MTHU identifiers
        return self.is_umls_cui(local_id) or self.is_umls_mthu_id(local_id)

    @staticmethod
    def is_metacyc_reaction_id(local_id: str) -> bool:
        # Allows: uppercase letters, digits, and hyphens (e.g., R13147, RXN-15029)
        return bool(re.match(r'^R[A-Z0-9-]*$', local_id)) or bool(re.match(r'^RXN-[0-9]+$', local_id))

    @staticmethod
    def is_inchikey_id(local_id: str) -> bool:
        # Allows: standard InChI key format (e.g., AMOFQIUOTAJRKS-UHFFFAOYSA-N)
        return bool(re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', local_id))

    @staticmethod
    def is_icd10_id(local_id: str) -> bool:
        # ICD-10 codes: letter followed by 2 digits, optional dot and alphanumeric
        return bool(re.match(r'^[A-Z]\d{2}(\.[A-Z0-9]+)?$', local_id))

    @staticmethod
    def is_icd9_id(local_id: str) -> bool:
        # ICD-9 codes: single code or range (e.g., 344.81 or 317-319.99)
        if '-' in local_id:
            # Range format: XXX-XXX.XX or XXX.XX-XXX.XX
            parts = local_id.split('-')
            if len(parts) != 2:
                return False
            return all(re.match(r'^\d{3}(\.\d{1,2})?$', part) for part in parts)
        else:
            # Single code format: XXX.XX
            return bool(re.match(r'^\d{3}(\.\d{1,2})?$', local_id))

    @staticmethod
    def is_go_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0004339 from GO:0004339)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_cl_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0000540 from CL:0000540)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_chebi_id(local_id: str) -> bool:
        # ChEBI IDs are positive integers
        return local_id.isdigit() and int(local_id) > 0

    @staticmethod
    def is_chembl_id(local_id: str) -> bool:
        # ChEMBL IDs: CHEMBL followed by digits
        return bool(re.match(r'^CHEMBL\d+$', local_id))

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

    @staticmethod
    def convert_float_to_int_str(local_id: str) -> str:
        if local_id.endswith('.0'):
            return str(int(float(local_id)))
        else:
            return local_id