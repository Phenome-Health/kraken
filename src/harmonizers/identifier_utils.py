"""
Identifier normalization utilities for converting SPOKE identifiers to Biolink CURIEs
"""

import json
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set


class IdentifierNormalizer:
    """Converts SPOKE identifiers to proper Biolink CURIEs"""
    
    def __init__(self):
        self.biolink_prefixes = self._load_biolink_prefixes()
        self.spoke_patterns = self._initialize_spoke_patterns()
    
    def _load_biolink_prefixes(self) -> Dict[str, str]:
        """Load Biolink model prefix map"""
        try:
            # Try to load from local copy first
            biolink_file = Path("/tmp/biolink_prefix_map.json")
            if not biolink_file.exists():
                # Fallback - could download or use embedded copy
                return self._get_embedded_biolink_prefixes()
            
            with open(biolink_file) as f:
                return json.load(f)
        except Exception:
            return self._get_embedded_biolink_prefixes()
    
    def _get_embedded_biolink_prefixes(self) -> Dict[str, str]:
        """Embedded subset of key Biolink prefixes for SPOKE data"""
        return {
            "UBERON": "http://purl.obolibrary.org/obo/UBERON_",
            "CL": "http://purl.obolibrary.org/obo/CL_", 
            "GO": "http://purl.obolibrary.org/obo/GO_",
            "DOID": "http://purl.obolibrary.org/obo/DOID_",
            "HGNC": "http://identifiers.org/hgnc/",
            "NCBIGene": "http://identifiers.org/ncbigene/",
            "ENSEMBL": "http://identifiers.org/ensembl/",
            "UniProtKB": "http://identifiers.org/uniprot/",
            "CHEBI": "http://purl.obolibrary.org/obo/CHEBI_",
            "PUBCHEM.COMPOUND": "http://identifiers.org/pubchem.compound/",
            "CHEMBL.COMPOUND": "http://identifiers.org/chembl.compound/",
            "DRUGBANK": "http://identifiers.org/drugbank/",
            "KEGG.COMPOUND": "http://identifiers.org/kegg.compound/",
            "KEGG.DRUG": "http://identifiers.org/kegg.drug/",
            "MESH": "http://identifiers.org/mesh/",
            "HP": "http://purl.obolibrary.org/obo/HP_",
            "OMIM": "http://identifiers.org/omim/",
            "DBSNP": "http://identifiers.org/dbsnp/",
            "CLINVAR": "http://identifiers.org/clinvar/",
            "EC": "http://identifiers.org/ec-code/",
            "REACTOME": "http://identifiers.org/reactome/",
            "WIKIPATHWAYS": "http://identifiers.org/wikipathways/",
            "PF": "http://identifiers.org/pfam/",
            "CVCL": "http://identifiers.org/cellosaurus/",
            "CLO": "http://purl.obolibrary.org/obo/CLO_",
            "ENVO": "http://purl.obolibrary.org/obo/ENVO_",
            "NCBITaxon": "http://purl.obolibrary.org/obo/NCBITaxon_",
            "GEONAMES": "http://identifiers.org/geonames/"
        }
    
    def _initialize_spoke_patterns(self) -> Dict[str, Dict]:
        """Initialize patterns for normalizing SPOKE identifiers by node type and source"""
        return {
            # Anatomy
            ("Anatomy", "Uberon"): {
                "prefix": "UBERON",
                "pattern": r"^UBERON:(\d+)$",
                "format": "UBERON:{id}"
            },
            
            # Genes
            ("Gene", "Entrez Gene"): {
                "prefix": "NCBIGene", 
                "pattern": r"^(\d+)$",
                "format": "NCBIGene:{id}"
            },
            ("Gene", "miRBase"): {
                "prefix": "MIRBASE",
                "pattern": r"^(\w+)$", 
                "format": "MIRBASE:{id}"
            },
            
            # Proteins  
            ("Protein", "UniProt"): {
                "prefix": "UniProtKB",
                "pattern": r"^([A-Z0-9]+)$",
                "format": "UniProtKB:{id}"
            },
            
            # Variants
            ("Variant", "dbSNP"): {
                "prefix": "DBSNP",
                "pattern": r"^(rs\d+)$", 
                "format": "DBSNP:{id}"
            },
            ("Variant", "unknown"): {
                "prefix": "DBSNP",
                "pattern": r"^(rs\d+)$",
                "format": "DBSNP:{id}"
            },
            
            # Diseases
            ("Disease", "Disease Ontology"): {
                "prefix": "DOID",
                "pattern": r"^DOID:(\d+)$",
                "format": "DOID:{id}"  
            },
            
            # Cell types
            ("CellType", "Cell Ontology"): {
                "prefix": "CL",
                "pattern": r"^CL:(\d+)$",
                "format": "CL:{id}"
            },
            
            # Cell lines
            ("CellLine", "Cellosaurus"): {
                "prefix": "CVCL", 
                "pattern": r"^(CLO_\d+)$",
                "format": "CVCL:{id}"
            },
            ("CellLine", "CellLineOntology"): {
                "prefix": "CLO",
                "pattern": r"^(CLO_\d+)$", 
                "format": "CLO:{id}"
            },
            
            # Biological processes, molecular functions, cellular components
            ("BiologicalProcess", "Gene Ontology"): {
                "prefix": "GO",
                "pattern": r"^GO:(\d+)$",
                "format": "GO:{id}"
            },
            ("MolecularFunction", "Gene Ontology"): {
                "prefix": "GO", 
                "pattern": r"^GO:(\d+)$",
                "format": "GO:{id}"
            },
            ("CellularComponent", "Gene Ontology"): {
                "prefix": "GO",
                "pattern": r"^GO:(\d+)$", 
                "format": "GO:{id}"
            },
            
            # Compounds
            ("Compound", "ChEBI"): {
                "prefix": "CHEBI",
                "extract_ids": ["CHEBI_ids"],
                "inchikey_fallback": True
            },
            ("Compound", "PubChem"): {
                "prefix": "PUBCHEM.COMPOUND",
                "extract_ids": ["pubchem_compound_ids"],
                "inchikey_fallback": True
            },
            ("Compound", "ChEMBL"): {
                "prefix": "CHEMBL.COMPOUND", 
                "extract_ids": ["chembl_ids"],
                "inchikey_fallback": True
            },
            ("Compound", "DrugBank"): {
                "prefix": "DRUGBANK",
                "extract_ids": ["drugbank_ids"], 
                "inchikey_fallback": True
            },
            ("Compound", "KEGG"): {
                "prefix": "KEGG.COMPOUND",
                "extract_ids": ["kegg_drug_ids"],
                "inchikey_fallback": True
            },
            
            # EC numbers
            ("EC", "ExplorEnz"): {
                "prefix": "EC",
                "pattern": r"^(\d+\.\d+\.\d+\.\d+)$",
                "format": "EC:{id}"
            },
            ("EC", "metacyc"): {
                "prefix": "EC", 
                "pattern": r"^(\d+\.\d+\.\d+\.[\dM]+)$",
                "format": "EC:{id}"
            },
            
            # Pathways
            ("Pathway", "unknown"): {
                "custom_handler": self._handle_reactome_pathway
            },
            ("Pathway", "WikiPathways"): {
                "prefix": "WIKIPATHWAYS",
                "pattern": r"^(WP\d+_r\d+)$",
                "format": "WIKIPATHWAYS:{id}"
            },
            
            # Protein domains
            ("ProteinDomain", "Pfam"): {
                "prefix": "PF",
                "pattern": r"^(PF\d+)$", 
                "format": "PF:{id}"
            },
            
            # MicroRNA
            ("MiRNA", "miRDB"): {
                "prefix": "MIRBASE",
                "pattern": r"^(hsa-miR-[\w-]+)$",
                "format": "MIRBASE:{id}"
            },
            
            # Organisms
            ("Organism", "ncbi-taxonomy"): {
                "prefix": "NCBITaxon",
                "pattern": r"^(\d+)$",
                "format": "NCBITaxon:{id}"
            },
            
            # Symptoms  
            ("Symptom", "MeSH"): {
                "prefix": "MESH",
                "pattern": r"^([DC]\d+)$",
                "format": "MESH:{id}" 
            },
            ("Symptom", "HPO"): {
                "prefix": "HP",
                "pattern": r"^([DC]\d+)$",  # May need adjustment
                "format": "HP:{id}"
            },
            
            # Environment
            ("Environment", "Environment Ontology"): {
                "prefix": "ENVO",
                "pattern": r"^ENVO_(\d+)$", 
                "format": "ENVO:{id}"
            }
        }
    
    def _handle_reactome_pathway(self, identifier: str, properties: Dict) -> str:
        """Custom handler for Reactome pathway identifiers"""
        if identifier.startswith("reactome:"):
            # Extract R-HSA-XXXXX part
            match = re.search(r"reactome:(R-HSA-\d+)", identifier)
            if match:
                return f"REACTOME:{match.group(1)}"
        return f"REACTOME:{identifier}"
    
    def normalize_spoke_identifier(self, node_type: str, source: str, 
                                   identifier: str, properties: Dict) -> Tuple[str, List[str]]:
        """
        Normalize SPOKE identifier to Biolink CURIE format
        
        Returns:
            Tuple of (primary_curie, equivalent_curies_list)
        """
        key = (node_type, source)
        
        # Check if we have a specific pattern for this type-source combination
        if key in self.spoke_patterns:
            pattern_config = self.spoke_patterns[key]
            
            # Handle custom processors
            if "custom_handler" in pattern_config:
                primary_id = pattern_config["custom_handler"](identifier, properties)
                return primary_id, [primary_id]
            
            # Handle compound nodes with extractable IDs
            if "extract_ids" in pattern_config:
                return self._handle_compound_identifiers(identifier, properties, pattern_config)
            
            # Handle standard pattern matching
            if "pattern" in pattern_config and "format" in pattern_config:
                match = re.match(pattern_config["pattern"], identifier)
                if match:
                    primary_id = pattern_config["format"].format(id=match.group(1))
                    return primary_id, [primary_id]
        
        # Fallback: try to detect common patterns
        return self._fallback_normalize(node_type, identifier, properties)
    
    def _handle_compound_identifiers(self, identifier: str, properties: Dict, 
                                     pattern_config: Dict) -> Tuple[str, List[str]]:
        """Handle compound nodes with multiple possible ID sources"""
        equivalent_ids = []
        primary_id = None
        
        # Extract IDs from specified properties
        for id_property in pattern_config.get("extract_ids", []):
            if id_property in properties and properties[id_property]:
                ids = properties[id_property]
                if isinstance(ids, list):
                    for id_val in ids:
                        if id_val:
                            curie = f"{pattern_config['prefix']}:{id_val}"
                            equivalent_ids.append(curie)
                            if not primary_id:
                                primary_id = curie
                elif ids:
                    curie = f"{pattern_config['prefix']}:{ids}"
                    equivalent_ids.append(curie)
                    if not primary_id:
                        primary_id = curie
        
        # Use InChI key as fallback if enabled and no other IDs found
        if not primary_id and pattern_config.get("inchikey_fallback"):
            if identifier.startswith("inchikey:"):
                inchi_key = identifier.replace("inchikey:", "")
                primary_id = f"INCHIKEY:{inchi_key}"
                equivalent_ids.append(primary_id)
        
        # Final fallback to original identifier with prefix
        if not primary_id:
            primary_id = f"{pattern_config['prefix']}:{identifier}"
            equivalent_ids.append(primary_id)
        
        return primary_id, equivalent_ids
    
    def _fallback_normalize(self, node_type: str, identifier: str, 
                            properties: Dict) -> Tuple[str, List[str]]:
        """Fallback normalization for unmatched patterns"""
        
        # If already a proper CURIE, return as-is
        if ":" in identifier and not identifier.startswith("inchikey:"):
            # Check if prefix is in Biolink model
            prefix = identifier.split(":")[0]
            if prefix.upper() in self.biolink_prefixes:
                return identifier, [identifier]
        
        # Generate a reasonable fallback based on node type
        fallback_prefixes = {
            "Gene": "NCBIGene",
            "Protein": "UniProtKB", 
            "Variant": "DBSNP",
            "Disease": "DOID",
            "Compound": "CHEBI",
            "Anatomy": "UBERON",
            "CellType": "CL",
            "BiologicalProcess": "GO",
            "MolecularFunction": "GO",
            "CellularComponent": "GO"
        }
        
        prefix = fallback_prefixes.get(node_type, "UNKNOWN")
        fallback_id = f"{prefix}:{identifier}"
        
        return fallback_id, [fallback_id]
    
    def extract_equivalent_identifiers(self, properties: Dict) -> List[str]:
        """Extract additional equivalent identifiers from node properties"""
        equivalent_ids = []
        
        # Common property names that contain equivalent IDs
        equivalent_id_properties = [
            "mesh_id", "mesh_ids", "mesh_list",
            "omim_list", "omim_ids", 
            "chembl_id", "chembl_ids",
            "drugbank_ids", "CHEBI_ids",
            "pubchem_compound_ids", "kegg_drug_ids",
            "ensembl", "refseq", "seqids",
            "accession", "accession_id", "xrefs",
            "ICD10", "SNOMEDCT", "equivalent_ids"
        ]
        
        for prop_name in equivalent_id_properties:
            if prop_name in properties and properties[prop_name]:
                prop_value = properties[prop_name]
                
                # Handle different formats
                if isinstance(prop_value, list):
                    for item in prop_value:
                        if item and str(item).strip():
                            equivalent_ids.extend(self._normalize_equivalent_id(prop_name, str(item)))
                elif isinstance(prop_value, str) and prop_value.strip():
                    equivalent_ids.extend(self._normalize_equivalent_id(prop_name, prop_value))
        
        return list(set(equivalent_ids))  # Remove duplicates
    
    def _normalize_equivalent_id(self, property_name: str, value: str) -> List[str]:
        """Normalize a single equivalent identifier based on its property name"""
        value = value.strip()
        if not value or value.lower() in ["null", "", "none"]:
            return []
        
        # Map property names to prefixes
        property_prefix_map = {
            "mesh_id": "MESH",
            "mesh_ids": "MESH",
            "omim_list": "OMIM", 
            "omim_ids": "OMIM",
            "chembl_id": "CHEMBL.COMPOUND",
            "chembl_ids": "CHEMBL.COMPOUND",
            "drugbank_ids": "DRUGBANK", 
            "CHEBI_ids": "CHEBI",
            "pubchem_compound_ids": "PUBCHEM.COMPOUND",
            "kegg_drug_ids": "KEGG.DRUG",
            "ensembl": "ENSEMBL",
            "ICD10": "ICD10",
            "SNOMEDCT": "SNOMEDCT"
        }
        
        prefix = property_prefix_map.get(property_name)
        if prefix:
            # Clean up common prefixes if already present
            clean_value = re.sub(r"^(CHEBI|MESH|OMIM|DB):", "", value)
            return [f"{prefix}:{clean_value}"]
        
        # For xrefs and other complex fields, try to parse
        if property_name == "xrefs":
            return self._parse_xref_value(value)
        
        return []
    
    def _parse_xref_value(self, xref: str) -> List[str]:
        """Parse cross-reference values that may contain multiple IDs"""
        ids = []
        
        # Split on common delimiters
        parts = re.split(r"[;,\|]", xref)
        for part in parts:
            part = part.strip()
            if ":" in part:
                # Looks like it's already a CURIE
                ids.append(part)
            
        return ids