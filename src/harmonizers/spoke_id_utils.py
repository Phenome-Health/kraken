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


class SpokeIDNormalizer:
    """Simple, robust identifier normalization focused on case handling"""
    
    def __init__(self, biolink_version: str):
        self.biolink_version = biolink_version
        self.normalized_prefixes_to_iris = self._load_prefix_to_iri_map()
        self.prefix_lowercase_map = {prefix.lower(): prefix for prefix in self.normalized_prefixes_to_iris.keys()}
        self.validator_prop = 'validator'
        self.cleaner_prop = 'cleaner'
        self.validator_map = self._load_validator_map()
        self.curie_construction_map = self._load_curie_construction_map()
        self.spoke_underscore_prefixes = {'SNOMED_', 'CLO_', 'ENVO_', 'CHR_', 'HPS_', 'CVCL_', 'BFO_'}
    

    def _load_prefix_to_iri_map(self) -> Dict[str, str]:
        """Load Biolink model prefix map"""
        logging.info(f"Grabbing biolink prefix map for version: {self.biolink_version}")
        url = f"https://raw.githubusercontent.com/biolink/biolink-model/refs/tags/v{self.biolink_version}/project/prefixmap/biolink-model-prefix-map.json"
        prefix_to_iri_map = load_biolink_file(url, self.biolink_version)

        # Remove prefixes as needed
        if 'KEGG' in prefix_to_iri_map:
            del prefix_to_iri_map['KEGG']  # We want to use only KEGG.COMPOUND, KEGG.REACTION, etc.

        # Add prefixes as needed (ones we're making up, that don't exist in biolink)
        prefix_to_iri_map['USZIPCODE'] = "https://www.unitedstateszipcodes.org/"
        prefix_to_iri_map['SMILES'] = "https://pubchem.ncbi.nlm.nih.gov/compound/"
        prefix_to_iri_map['CVCL'] = "https://web.expasy.org/cellosaurus/CVCL_"
        prefix_to_iri_map['VESICLEPEDIA'] = "http://microvesicles.org/exp_summary?exp_id="
        prefix_to_iri_map['NDFRT'] = "http://purl.bioontology.org/ontology/NDFRT/"
        prefix_to_iri_map['BVBRC'] = "https://www.bv-brc.org/view/Genome/"
        prefix_to_iri_map['GeoNames'] = "http://www.geonames.org/search.html?q="  # Note: this doesn't go exactly to page for item, but closest I could find
        prefix_to_iri_map['NHANES'] = "https://dsld.od.nih.gov/label/"  # These IRIs work, but weirdly SPOKE's identifiers for these nodes don't match what they have..
        prefix_to_iri_map['MIRDB'] = "https://mirdb.org/cgi-bin/mature_mir.cgi?name="  # Not sure if it's right to use a 'mature' iri like this for all...
        prefix_to_iri_map['CYTOBAND'] = ""  # Haven't found good iri for these yet..
        prefix_to_iri_map['CHR'] = ""  # Country Health Rankings.. Haven't found good iri for these yet
        prefix_to_iri_map['AHRQ'] = ""  # AHRQ SDOH Database
        prefix_to_iri_map['HPS'] = ""  # Household Pulse Survey
        prefix_to_iri_map['mirbase'] = "https://mirbase.org/hairpin/"  # Biolink has mirbase in here, but their iri doesn't work
        prefix_to_iri_map['metacyc.pathway'] = "https://metacyc.org/pathway?orgid=META&id="  # Biolink has metacyc.reaction, but not pathway
        prefix_to_iri_map['metacyc.ec'] = "https://biocyc.org/META/NEW-IMAGE?type=EC-NUMBER&object=EC-"  # Biolink has metacyc.reaction, but not ec (these are like provisional ec codes, not yet in explorenz)
        prefix_to_iri_map['FIPS.PLACE'] = ""
        prefix_to_iri_map['FIPS.STATE'] = ""
        prefix_to_iri_map['PHARMVAR'] = ""  # This wants a number rather than the symbol..
        prefix_to_iri_map['CDCSVI'] = ""  # CDC Social Vulnerability Index

        # Override prefixes as needed (if Biolink's iri is broken)
        prefix_to_iri_map['OMIM'] = "https://omim.org/entry/"
        prefix_to_iri_map['REACT'] = "https://reactome.org/content/detail/" # Works for Complexes and Pathways (I think)

        return prefix_to_iri_map


    def _load_validator_map(self) -> Dict[str, Dict[str, Callable]]:
        validator = self.validator_prop
        cleaner = self.cleaner_prop
        return {
            'ahrq': {validator: self.is_ahrq_id},
            'bfo': {validator: self.is_bfo_id},
            'bvbrc': {validator: self.is_bvbrc_id},
            'cdcsvi': {validator: self.is_cdcsvi_id},
            'chebi': {validator: self.is_chebi_id},
            'chembl.compound': {validator: self.is_chembl_compound_id},
            'chembl.target': {validator: self.is_chembl_target_id},
            'chr': {validator: self.is_chr_id},
            'cl': {validator: self.is_cl_id},
            'clo': {validator: self.is_clo_id},
            'complexportal': {validator: self.is_complexportal_id},
            'cvcl': {validator: self.is_cellosaurus_id},
            'cytoband': {validator: self.is_cytoband_id},
            'dbsnp': {validator: self.is_dbsnp_id},
            'doid': {validator: self.is_doid_id},
            'drugbank': {validator: self.is_drugbank_id},
            'ec': {validator: self.is_ec_id},
            'envo': {validator: self.is_envo_id},
            'fips.place': {validator: self.is_fips_compound_id},
            'fips.state': {validator: self.is_fips_state_id},
            'geonames': {validator: self.is_geonames_id},
            'go': {validator: self.is_go_id},
            'hps': {validator: self.is_hps_id},
            'icd9': {validator: self.is_icd9_id},
            'icd10': {validator: self.is_icd10_id},
            'inchikey': {validator: self.is_inchikey_id},
            'kegg.compound': {validator: self.is_kegg_compound_id},
            'kegg.drug': {validator: self.is_kegg_drug_id},
            'kegg.reaction': {validator: self.is_kegg_reaction_id},
            'loinc': {validator: self.is_loinc_id},
            'mesh': {validator: self.is_mesh_id},
            'metacyc.ec': {validator: self.is_metacyc_ec_id},
            'metacyc.pathway': {validator: self.is_metacyc_pathway_id},
            'metacyc.reaction': {validator: self.is_metacyc_reaction_id},
            'mirbase': {validator: self.is_mirbase_id},
            'mirdb': {validator: self.is_mirdb_id},
            'ncbigene': {validator: self.is_ncbigene_id},
            'ncbitaxon': {validator: self.is_ncbitaxon_id},
            'ndfrt': {validator: self.is_ndfrt_id},
            'nhanes': {validator: self.is_nhanes_id},
            'omim': {validator: self.is_omim_id},
            'pfam': {validator: self.is_pfam_id},
            'pharmvar': {validator: self.is_pharmvar_id},
            'pubchem.compound': {validator: self.is_pubchem_compound_id},
            'react': {validator: self.is_reactome_id},
            'smiles': {validator: self.is_smiles_string},
            'snomedct': {validator: self.is_snomedct_id, cleaner: self.clean_snomed_id},
            'uberon': {validator: self.is_uberon_id},
            'umls': {validator: self.is_umls_id},
            'uniprotkb': {validator: self.is_uniprot_id},
            'uszipcode': {validator: self.is_uszipcode_id, cleaner: self.clean_zipcode},
            'vesiclepedia': {validator: self.is_vesiclepedia_id},
            'wikipathways': {validator: self.is_wikipathways_id, cleaner: self.clean_wikipathways_id},
        }

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
            ('Compound', 'standardized_smiles'): 'smiles',
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
                    # Grab the proper validation and cleaning functions
                    validate = self.validator_map[prefix][self.validator_prop]
                    clean = self.validator_map[prefix].get(self.cleaner_prop)

                    # Clean up the id as necessary
                    if clean:
                        local_id = clean(local_id)

                    # Stop if we've found a prefix that applies
                    if validate(local_id) in {True, None}:  # It's ok if it's a 'known invalid' format
                        chosen_prefix = prefix
                        break
            else:
                chosen_prefix = prefix_entry
                # Clean up the local ID, if needed
                clean = self.validator_map[chosen_prefix].get(self.cleaner_prop)
                if clean:
                    local_id = clean(local_id)

            if chosen_prefix:
                curie, iri = self.construct_curie(chosen_prefix, local_id)
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


    def construct_curie(self, prefix_lowercase: str, local_id: str) -> Tuple[str, str]:
        validate = self.validator_map[prefix_lowercase][self.validator_prop]
        is_valid_id = validate(local_id)
        if is_valid_id:
            prefix_normalized = self.prefix_lowercase_map[prefix_lowercase]
            iri_root = self.normalized_prefixes_to_iris[prefix_normalized]
            iri = f"{iri_root}{local_id}" if iri_root else ""
            return f"{prefix_normalized}:{local_id}", iri
        elif is_valid_id is None:
            # Indicates this is a known invalid ID format for this node_type, source pair; we will skip it
            logging.warning(f"Local id '{local_id}' is invalid for {prefix_lowercase} (known invalid format)")
            return KNOWN_INVALID, ''
        else:
            # This is an unknown invalid ID format; we want to return nothing, which will halt processing after logging
            logging.error(f"Local id '{local_id}' is invalid for {prefix_lowercase} (UNKNOWN invalid format)")
            return '', ''
    

    def extract_equivalent_identifiers(self, node_type: str, properties: Dict) -> List[str]:
        """Extract equivalent IDs from properties - simplified approach"""
        equivalent_ids = set()
        none_strings = {'null', 'none', 'nan'}
        equiv_id_sources = {'chembl', 'drugbank', 'chebi', 'pubchem', 'kegg', 'mesh', 'ensembl', 'omim'}
        exact_fields = {'standardized_smiles', 'snomedct', 'icd10', 'icd9', 'accession'}
        
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
        # Allows: D, C, or M followed by one or more digits
        return bool(re.match(r'^[DCM]\d+$', local_id))

    @staticmethod
    def is_metacyc_ec_id(local_id: str) -> bool:
        # Allows three digits groups, then a final group of alphanumeric characters
        return bool(re.match(r'^\d+\.\d+\.\d+\.[a-zA-Z0-9]+$', local_id))

    @staticmethod
    def is_metacyc_reaction_id(local_id: str) -> bool:
        # Allows: Hyphen-separated uppercase/numeric or capitalized alpha parts; must contain 'RXN' somewhere
        # e.g., 3.2.1.68-RXN, TRANS-RXN0-593, CYPRIDINA-LUCIFERIN-2-MONOOXYGENASE-RXN, RXN0-5258-Yeast
        has_valid_chars = bool(re.match(r'^[A-Za-z0-9-.+]+$', local_id))
        parts = local_id.split('-')
        has_proper_capitalization = all(part.isupper() or (not any(char.isalpha() for char in part)) or (part.isalpha()) for part in parts)
        return has_valid_chars and has_proper_capitalization and 'RXN' in local_id

    @staticmethod
    def is_metacyc_pathway_id(local_id: str) -> Optional[bool]:
        # MetaCyc pathway IDs: examples: PWY-#### or PWY0-#### or DESCRIPTIVE-NAME-PWY or PWY18C3-9
        has_valid_chars = bool(re.match(r'^[A-Z0-9-+]+$', local_id))
        is_metacyc_id = has_valid_chars and local_id.isupper()
        if is_metacyc_id:
            return True
        elif any(char.isalpha() for char in local_id):
            return None  # Meant to catch english names given as identifiers, like Glycan biosynthesis - 2
        else:
            return False

    @staticmethod
    def is_snomedct_id(local_id: str) -> Optional[bool]:
        # SNOMED CT IDs are numeric strings
        if 'e+' in local_id:  # Known spoke bug where some snomed ct IDs are in scientific notation, like '1.62248710001191e+16'
            return None
        else:
            return local_id.isdigit()

    @staticmethod
    def is_cellosaurus_id(local_id: str) -> bool:
        # Allows: Exactly 4 digits or uppercase letters
        return bool(re.match(r'^[A-Z0-9]{4}$', local_id))

    @staticmethod
    def is_cytoband_id(local_id: str) -> bool:
        # Allows: chromosome (number or X/Y), arm (p/q), band, and optional sub-band e.g., 1p36.33
        return bool(re.match(r'^(\d{1,2}|[XYxy])[pq]\d+(\.\d+)?$', local_id))

    @staticmethod
    def is_mirbase_id(local_id: str) -> bool:
        # Allows: MI or MIMAT followed by exactly 7 digits
        return bool(re.match(r'^(MI|MIMAT)\d{7}$', local_id))

    @staticmethod
    def is_mirdb_id(local_id: str) -> bool:
        # Allows: 3 lowercase letters, an optional miR-, followed by a mix of lowercase letters, digits, and hyphens
        return bool(re.match(r'^[a-z]{3}-(miR-)?[-a-z0-9]+$', local_id))
    @staticmethod
    def is_uberon_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0003233 from UBERON:0003233)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_dbsnp_id(local_id: str) -> Optional[bool]:
        # Allows: rs followed by digits, with an optional version suffix (e.g., .1)
        if local_id == '-':  # Some identifiers are just a hyphen; think it's like a null value
            return None
        elif ',' in local_id:
            # Detect when multiple dbsnp IDs are concatenated into one ID (skip these for now)
            parts = local_id.split(',')
            if bool(re.match(r'^rs[0-9]+(\.\d+)?$', parts[0].strip())):
                return None
        return bool(re.match(r'^rs[0-9]+(\.\d+)?$', local_id))

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
    def is_envo_id(local_id: str) -> bool:
        # Allows: one or more digits
        return bool(re.match(r'^\d+$', local_id))

    @staticmethod
    def is_reactome_id(local_id: str) -> Optional[bool]:
        # Allows: R-HSA-digits (e.g., R-HSA-162582)
        if local_id == 'root':
            return None
        else:
            return bool(re.match(r'^R-[A-Z]{3}-[0-9]+$', local_id))

    @staticmethod
    def is_kegg_reaction_id(local_id: str) -> bool:
        # Allows: R followed by exactly 5 digits
        return bool(re.match(r'^R\d{5}$', local_id))

    @staticmethod
    def is_kegg_drug_id(local_id: str) -> bool:
        # Allows: D followed by exactly 5 digits
        return bool(re.match(r'^D\d{5}$', local_id))

    @staticmethod
    def is_kegg_compound_id(local_id: str) -> bool:
        # Allows: C followed by exactly 5 digits
        return bool(re.match(r'^C\d{5}$', local_id))

    @staticmethod
    def is_pubchem_compound_id(local_id: str) -> bool:
        # Allows: positive integers (PubChem CIDs are numeric)
        return local_id.isdigit() and int(local_id) > 0

    @staticmethod
    def is_smiles_string(local_id: str) -> bool:
        # Basic SMILES validation: multiple element symbols and valid characters
        valid_chars = set('BCNOPSFHIKLMWUVYXZbcnopslr[]()=#%+\\/@.-0123456789')
        has_valid_chars = all(c in valid_chars for c in local_id)
        return has_valid_chars

    @staticmethod
    def is_wikipathways_id(local_id: str) -> bool:
        # Allows: WP followed by digits
        return bool(re.match(r'^WP[0-9]+$', local_id))

    @staticmethod
    def is_vesiclepedia_id(local_id: str) -> bool:
        # Allows: one or more digits
        return bool(re.match(r'^\d+$', local_id))

    @staticmethod
    def clean_wikipathways_id(local_id: str) -> str:
        # Get rid of version suffix info, like in WP5395_r126912
        return local_id.split('_')[0]

    @staticmethod
    def is_doid_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0070557 from DOID:0070557)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_drugbank_id(local_id: str) -> bool:
        # Allows: DB followed by exactly 5 digits
        return bool(re.match(r'^DB\d{5}$', local_id))

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
        # Allows: PF or CL followed by digits
        return bool(re.match(r'^(PF|CL)\d+$', local_id))

    @staticmethod
    def is_pharmvar_id(local_id: str) -> bool:
        # Allows: Gene symbol, asterisk, allele number, and optional sub-allele - e.g., CYP26A1*1.001
        return bool(re.match(r'^[A-Z0-9]+\*\d+(\.\d+)?$', local_id))

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
    def is_uniprot_protein_id(local_id: str) -> bool:
        # Allows: 6 or 10 uppercase alphanumeric characters, requires >=1 letter and >=1 digit
        is_proper_alphanumeric = bool(re.match(r'^([A-Z0-9]{6}|[A-Z0-9]{10})$', local_id))
        has_letter = any(c.isalpha() for c in local_id)
        has_digit = any(c.isdigit() for c in local_id)
        return is_proper_alphanumeric and has_letter and has_digit

    @staticmethod
    def is_uniprot_feature_id(local_id: str) -> bool:
        # Allows: UniProt ID, hyphen, then PRO_ and digits
        pattern = r'^([A-Z0-9]{6}|[A-Z0-9]{10})-PRO_\d+$'
        return bool(re.match(pattern, local_id))

    def is_uniprot_id(self, local_id: str) -> bool:
        # Allows: Regular uniprot protein IDs or the special 'feature' ids
        return self.is_uniprot_protein_id(local_id) or self.is_uniprot_feature_id(local_id)

    @staticmethod
    def is_inchikey_id(local_id: str) -> bool:
        # Allows: standard InChI key format (e.g., AMOFQIUOTAJRKS-UHFFFAOYSA-N)
        return bool(re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', local_id))

    @staticmethod
    def is_icd10_id(local_id: str) -> bool:
        # ICD-10 codes: letter followed by 2 letters or digits, optional dot and alphanumeric
        return bool(re.match(r'^[A-Z][A-Z0-9]{2}(\.[A-Z0-9]+)?$', local_id))

    @staticmethod
    def is_go_id(local_id: str) -> bool:
        # Allows: exactly 7 digits
        return bool(re.match(r'^\d{7}$', local_id))

    @staticmethod
    def is_hps_id(local_id: str) -> bool:
        # Allows: one or more alphabetic characters
        return bool(re.match(r'^[a-zA-Z_]+$', local_id))

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
    def is_chr_id(local_id: str) -> bool:
        # Allows: lowercase letters, digits, and underscores; requires at least one letter
        has_valid_chars = bool(re.match(r'^[a-z0-9_/-]+$', local_id))
        return has_valid_chars and any(char.isalpha() for char in local_id)

    @staticmethod
    def is_cl_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0000540 from CL:0000540)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_clo_id(local_id: str) -> Optional[bool]:
        # Allows: Exactly 7 digits
        if local_id.startswith('http://'):
            return None  # A couple nodes have old URLs as IDs
        else:
            return bool(re.match(r'^\d{7}$', local_id))

    @staticmethod
    def is_complexportal_id(local_id: str) -> bool:
        # Allows: CPX- followed by one or more digits
        return bool(re.match(r'^CPX-\d+$', local_id))

    @staticmethod
    def is_ahrq_id(local_id: str) -> bool:
        # Allows: uppercase letters, digits, and underscores
        return bool(re.match(r'^[A-Z0-9_]+$', local_id))

    @staticmethod
    def is_bfo_id(local_id: str) -> bool:
        # Allows: one or more digits
        return bool(re.match(r'^\d+$', local_id))

    @staticmethod
    def is_bvbrc_id(local_id: str) -> bool:
        # Allows: digits, a period, and more digits
        return bool(re.match(r'^\d+\.\d+$', local_id))

    @staticmethod
    def is_cdcsvi_id(local_id: str) -> bool:
        # Allows: one or more uppercase alphabetic characters
        return bool(re.match(r'^[A-Z]+$', local_id))

    @staticmethod
    def is_chebi_id(local_id: str) -> Optional[bool]:
        # ChEBI IDs are positive integers
        if local_id == 'root':
            return None
        else:
            return local_id.isdigit() and int(local_id) > 0

    @staticmethod
    def is_chembl_compound_id(local_id: str) -> bool:
        # Allows: CHEMBL followed by digits
        return bool(re.match(r'^CHEMBL\d+$', local_id))

    @staticmethod
    def is_chembl_target_id(local_id: str) -> bool:
        # Allows: CHEMBL followed by one or more digits
        return bool(re.match(r'^CHEMBL\d+$', local_id))

    @staticmethod
    def is_hpo_id(local_id: str) -> bool:
        # Allows: digits only (e.g., 0001234 from HP:0001234)
        return bool(re.match(r'^[0-9]+$', local_id))

    @staticmethod
    def is_uszipcode_id(local_id: str) -> bool:
        # Allows: 5-digit US ZIP codes
        return bool(re.match(r'^[0-9]{5}$', local_id)) or local_id == 'US'

    @staticmethod
    def is_fips_compound_id(local_id: str) -> bool:
        # Allows: 6, 7, 11, or 12 digit FIPS-like codes
        return bool(re.match(r'^(\d{6}|\d{7}|\d{11}|\d{12})$', local_id))

    @staticmethod
    def is_fips_state_id(local_id: str) -> bool:
        # Allows: exactly 2 digits
        return bool(re.match(r'^\d{2}$', local_id))

    @staticmethod
    def is_geonames_id(local_id: str) -> bool:
        # Allows: 2-letter country code OR 2 letters followed by one or more dot-separated alphanumeric segments
        return bool(re.match(r'^[A-Z]{2}$', local_id)) or bool(re.match(r'^[A-Z]{2}(\.[A-Z0-9]+)+$', local_id))

    @staticmethod
    def is_ndfrt_id(local_id: str) -> bool:
        # Allows: N followed by exactly 10 digits
        return bool(re.match(r'^N\d{10}$', local_id))

    @staticmethod
    def is_nhanes_id(local_id: str) -> bool:
        # Allows: one or more digits
        return bool(re.match(r'^\d+$', local_id))

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

    def clean_snomed_id(self, local_id: str) -> str:
        return self.convert_float_to_int_str(local_id)

    @staticmethod
    def clean_zipcode(local_id: str) -> str:
        # Strip prefix off of zipcodes like AZ-85039
        return local_id.split('-')[-1]