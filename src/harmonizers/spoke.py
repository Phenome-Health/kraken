"""
SPOKE harmonizer - converts SPOKE format to unified Biolink schema
"""

from pathlib import Path
import sys
from typing import List, Optional, Tuple
import jsonlines
import logging
from ..utils.metagraph import generate_metagraph_for_source
from ..utils.constants import *
from .spoke_id_utils import SpokeIDNormalizer


def harmonize_spoke(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str, build_metagraph: bool):
    """Harmonize SPOKE mixed JSONL to unified Biolink schema using streaming"""
    logging.info(f"Harmonizing SPOKE: {input_file} -> {nodes_output}, {edges_output}")

    node_count = 0
    edge_count = 0
    
    # Initialize identifier normalizer
    id_norm = SpokeIDNormalizer(biolink_version=biolink_version)
    
    # Keep track of normalized node IDs for edge mapping
    spoke_to_normalized_id = {}

    with jsonlines.open(input_file, 'r') as reader, \
         jsonlines.open(nodes_output, 'w') as nodes_writer, \
         jsonlines.open(edges_output, 'w') as edges_writer:
        
        for item in reader:
            item_type = item.get('type')
            
            if item_type == 'node':
                harmonized_node = harmonize_spoke_node(item, id_norm)

                if harmonized_node:  # Occasionally we skip nodes (if invalid identifier, etc...)
                    # Store mapping for edge processing
                    spoke_to_normalized_id[item['id']] = harmonized_node[ID]
                    
                    nodes_writer.write(harmonized_node)

                    node_count += 1
                    if node_count % 1000000 == 0:
                        logging.info(f"Processed {node_count} SPOKE nodes")

            elif item_type == 'relationship':
                harmonized_edge = harmonize_spoke_edge(item, spoke_to_normalized_id)
                if harmonized_edge:
                    edges_writer.write(harmonized_edge)

                    edge_count += 1
                    if edge_count % 1000000 == 0:
                        logging.info(f"Processed {edge_count} SPOKE edges")

    logging.info(f"SPOKE harmonization complete: {node_count} nodes, {edge_count} edges")

    if build_metagraph:
        # Generate metagraph for harmonized output
        # Store metagraphs in artifacts/metagraphs/harmonized/source_name/
        artifacts_root = Path("artifacts")
        metagraph_dir = artifacts_root / "metagraphs" / "harmonized" / "spoke"
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, "spoke")
        logging.info("SPOKE metagraph generated")


def harmonize_spoke_node(node_item: dict, id_norm: SpokeIDNormalizer) -> dict:
    """Harmonize a single SPOKE node"""
    properties = node_item.get('properties', {})
    labels = node_item.get('labels', [])
    if not labels:
        raise ValueError(f"SPOKE node is missing labels: {node_item}")
    
    # Get node type and source for identifier normalization
    node_type = labels[0]  # Primary label
    
    # Get source(s) - handle both 'source' and 'sources' 
    primary_source, secondary_sources = get_all_sources(node_item)
    
    original_identifier = properties.get('identifier', node_item['id'])

    # Skip invalid/unhelpful nodes
    if primary_source == 'CellLineOntology' and original_identifier.startswith('http'):
        return None  # Example of such a node's 'identifier': 'http://www.ebi.ac.uk/cellline#cancer_cell_line'
    if primary_source == 'CDC/ATSDR Social Vulnerability Index':
        return None  # There's only one node from this source in SPOKE v6, and it doesn't exactly map to the right mesh term

    # Normalize the identifier
    normalized_id = id_norm.normalize_spoke_identifier(node_type, primary_source, original_identifier, properties)
    
    # Extract additional equivalent identifiers from properties
    additional_equivalent_ids = id_norm.extract_equivalent_identifiers(node_type, properties)
    all_equivalent_ids = list(set([normalized_id] + additional_equivalent_ids))
    
    harmonized_node = {
        ID: normalized_id,
        CATEGORIES: map_spoke_labels_to_biolink(labels, primary_source),
        PROVIDED_BY: [SPOKE_INFORES],
        EQUIVALENT_IDS: all_equivalent_ids,
        'spoke_nodes': [node_item]  # Make this a list, because we want it merged during entity resolution (e.g., SPOKE maps separate nodes to the same KEGG.REACTION identifier)
    }
    if properties.get(NAME):
        harmonized_node[NAME] = properties[NAME]
        harmonized_node[SYNONYMS] = [harmonized_node[NAME]]  # TODO: down the road see about extracting SPOKE synonyms
    
    return harmonized_node


def harmonize_spoke_edge(edge_item: dict, spoke_to_normalized_id: dict) -> Optional[dict]:
    """Harmonize a single SPOKE edge"""
    edge_type = edge_item.get('label')
    if not edge_type:
        raise ValueError(f"SPOKE edge is missing type: {edge_item}")

    spoke_subject_id = edge_item['start']['id']
    spoke_object_id = edge_item['end']['id']
    try:
        predicate, subject_id, object_id, qual_predicate, qual_direction, qual_aspect = map_spoke_edge_type_to_biolink(edge_type, 
                                                                                                                       spoke_subject_id, 
                                                                                                                       spoke_object_id)
    except:
        logging.error(f"Failed to find biolink type for spoke edge type '{edge_type}': {edge_item}")
        sys.exit(1)
    
    # Get source(s) - handle both 'source' and 'sources' 
    primary_source, secondary_sources = get_all_sources(edge_item)

    # Map SPOKE internal IDs to normalized CURIEs
    if subject_id in spoke_to_normalized_id and object_id in spoke_to_normalized_id:
        normalized_subject_id = spoke_to_normalized_id[subject_id]
        normalized_object_id = spoke_to_normalized_id[object_id]
    else:
        logging.warning(f"No normalized IDs available for {subject_id} and/or {object_id}. Skipping this edge.")
        return None
    
    # Remove the full start/end node objects (replace with their SPOKE IDs instead - saves a lot of space)
    edge_item['start'] = subject_id
    edge_item['end'] = object_id

    harmonized_edge = {
        SUBJECT: normalized_subject_id,
        OBJECT: normalized_object_id,
        PREDICATE: predicate,
        PRIMARY_KS: primary_source,  # TODO: Convert to infores curies where possible...
        AGGREGATOR_KS: SPOKE_INFORES,
        'spoke_edges': [edge_item]  # Make this a list, because we want it merged during entity resolution (e.g., SPOKE maps separate nodes to the same KEGG.REACTION identifier, creating duplicate edges)
    }
    # Tack on any additional sources
    if secondary_sources:
        harmonized_edge[SUPPORTING_SOURCES] = secondary_sources  # TODO: Convert to infores curies where possible...
    # Tack on any qualifiers
    if qual_predicate:
        harmonized_edge[QUALIFIED_PREDICATE] = qual_predicate
    if qual_direction:
        harmonized_edge[QUALIFIED_DIRECTION] = qual_direction
    if qual_aspect:
        harmonized_edge[QUALIFIED_ASPECT] = qual_aspect
    
    return harmonized_edge


def map_spoke_labels_to_biolink(labels: List[str], source: str) -> List[str]:
    """Map SPOKE node labels to Biolink categories"""

    # Simple mapping - extend as needed
    label_mapping = {
        'Compound': 'ChemicalEntity',  # This more accurate than SmallMolecule?
        'Variant': 'SequenceVariant',
        'Organism--ncbi-taxonomy': 'OrganismTaxon',  # Need to consider source for these
        'Organism--BV-BRC': 'OrganismalEntity',  # Need to consider source for these
        'Protein': 'Protein',
        'Location': 'GeographicLocation',
        'ClinicalLab': 'ClinicalFinding',
        'Reaction': 'MolecularActivity',
        'Gene': 'Gene',
        'ProteinDomain': 'ProteinDomain',
        'DietarySupplement': 'Food',  # Doesn't seem to be a good Biolink type for supplements... invent one? add 'supplement' flag?
        'Anatomy': 'AnatomicalEntity',
        'BiologicalProcess': 'BiologicalProcess',
        'CellLine': 'CellLine',
        'Disease': 'Disease',
        'EC': 'BiologicalEntity',  # TODO: Collapse these nodes into annotations on Protein nodes (see issue)
        'PwGroup': 'MacromolecularComplex',  # Think these are protein working groups?
        'Pathway': 'Pathway',
        'SideEffect': 'DiseaseOrPhenotypicFeature',  # The fact that this is a side effect would be implied by the edge from the drug
        'Blend': 'ChemicalMixture',
        'MolecularFunction': 'MolecularActivity',
        'CellType': 'Cell',
        'MiRNA': 'MicroRNA',
        'Complex': 'MacromolecularComplex',
        'Symptom': 'PhenotypicFeature',
        'Haplotype': 'Haplotype',
        'CellularComponent': 'CellularComponent',
        'ExtracellularParticle': 'CellularComponent',  # Definition says "in or around" the cell.. ok?
        'SDoH': 'SocioeconomicExposure',
        'Cytoband': 'GenomicEntity',  # This reasonable?
        'ProteinFamily': 'ProteinFamily',
        'PharmacologicClass': 'Drug',  # Need to flag these as *classes* of drugs somehow?
        'Environment': 'EnvironmentalFeature'
    }
    biolink_types = set()
    for label in labels:
        if label in label_mapping:
            node_type = label_mapping[label]
        else:
            node_type = label_mapping[f"{label}--{source}"]
        biolink_types.add(f"{BIOLINK_PREFIX}:{node_type}")
    return list(biolink_types)


def map_spoke_edge_type_to_biolink(edge_type: str, 
                                   original_subject_id: str, 
                                   original_object_id: str) -> Tuple[str, str, str, Optional[str], Optional[str], Optional[str]]:
    """Map SPOKE edge types to Biolink predicates; flips edges as necessary to use canonical predicates."""
    core_edge_type = '_'.join(edge_type.split('_')[:-1])  # Gets rid of suffix indicating node categories, like _GiP

    # Simple mapping - extend as needed
    predicate = 'type'
    flip = 'flip'
    type_map = {
        'PREVALENCEIN': {predicate: 'associated_with'},  # SPOKE only uses for SocioeconomicExposure-->GeographicLocation edges TODO: would occurs_in be better?
        'INTERACTS': {predicate: 'interacts_with'},
        'INTERACTS_as_LR': {predicate: 'physically_interacts_with'},
        'MAPS': {predicate: 'is_sequence_variant_of'},  # SPOKE only uses for SequenceVariant-->Gene edges
        'TARGETS': {predicate: 'regulates'},  # SPOKE only uses for MiRNA-->Gene edges
        'EXPRESSES': {predicate: 'expressed_in', flip: True},  # SPOKE only uses for Anatomy-->Gene edges
        'DOWNREGULATES': {predicate: 'regulates', QUALIFIED_DIRECTION: 'downregulated'},  # Anatomy-->Gene, Compound-->Gene, Gene-->Gene, Variant-->Gene
        'BINDS': {predicate: 'binds'},  # Compound-->Protein/ProteinDomain
        'UPREGULATES': {predicate: 'regulates', QUALIFIED_DIRECTION: 'upregulated'},  # Protein-->Gene
        'REGULATES': {predicate: 'regulates'},
        'EXPRESSEDIN': {predicate: 'expressed_in'},  # Gene/Protein-->CellType
        'EXPRESSEDIN_GeiD': {predicate: 'gene_associated_with_condition'},  # Gene-->Disease
        'PARTICIPATES': {predicate: 'has_participant', flip: True},  # Gene-->BioProcess, Gene-->CellComponent, etc. TODO: these violate domain/range! 
        'BELONGS': {predicate: 'related_to'},  # Variant-->Gene, Haplotype-->Gene, Variant-->Haplotype  TODO: need to improve this.. has_part?
        'ASSOCIATES': {predicate: 'associated_with'},
        'ISA': {predicate: 'subclass_of'},
        'PARTOF': {predicate: 'has_part',  flip: True},
        'ENCODES_GeP': {predicate: 'gene_product_of', flip: True},  # Gene-->Protein
        'ENCODES_GeM': {predicate: 'gene_product_of', flip: True},  # Gene-->MicroRNA
        'ENCODES_OeP': {predicate: 'produces'},  # Organism-->Protein  TODO: Get a check on this?
        'HAS': {predicate: 'related_to'},  # Protein-->EC, CellLine-->Variant  TODO: improve this?
        'PREVALENCE': {predicate: 'associated_with'},  # Disease-->GeographicLocation... TODO: would occurs_in be better?
        'ISOLATEDIN': {predicate: 'located_in'},  # Organism-->GeographicLocation
        'RESPONDS_TO': {predicate: 'affects', flip: True},  # Organism-->Compound
        'CAUSES': {predicate: 'causes'},
        'CONTAINS': {predicate: 'has_part'},
        'FOUNDIN': {predicate: 'located_in'},
        'PRESENTS': {predicate: 'has_phenotype'},  # Disease-->Symptom
        'RESEMBLES': {predicate: 'similar_to'},  # Disease-->Disease
        'MEASURESIN': {predicate: 'related_to'},  # ClinicalLab-->Anatomy  TODO: Flag these edges as being measured_in?
        'LOCALIZES': {predicate: 'disease_has_location'},  # Disease-->Anatomy
        'PRODUCES': {predicate: 'produces'},  # Reaction-->Compound
        'MEASURES': {predicate: 'assesses'},
        'TREATS': {predicate: 'treats'},
        'CONSUMES': {predicate: 'consumes'},
        'HASROLE': {predicate: 'has_chemical_role'},  # Compound-->Compound
        'MENTIONED_CLINICAL_TRIALS_FOR': {predicate: 'mentioned_in_clinical_trials_for'},  # TODO: ask Gwenlyn what she does with these mentions..
        'IN_CLINICAL_TRIALS_FOR': {predicate: 'in_clinical_trials_for'},
        'CATALYZES': {predicate: 'catalyzes'},  # EC-->Reaction
        'MARKER_NEG': {predicate: 'exacerbates_condition'},  # Gene-->Disease. technically gene is not in the domain, but hard to find a better one..
        'MARKER_POS': {predicate: 'ameliorates_condition'},  # Gene-->Disease. technically gene is not in the domain, but hard to find a better one..
        'AFFECTS': {predicate: 'affects'},
        'INCREASEDIN': {predicate: 'has_increased_amount', flip: True},  # Protein-->Disease
        'CONTRAINDICATES': {predicate: 'contraindicated_in'},  # Compound-->Disease
        'MORTALITY': {predicate: 'occurs_in'},  # Disease-->GeographicLocation
        'MEMBEROF': {predicate: 'has_member', flip: True},  # ProteinDomain-->ProteinFamily
        'INCLUDES': {predicate: 'has_member'},  # PharmacologicClass-->Compound
        'AFFECT': {predicate: 'affects_response_to'},  # Variant-->Compound
        'TRANSPORTS': {predicate: 'affects', QUALIFIED_ASPECT: 'transport'},  # Protein-->Compound
        'DERIVES_FROM': {predicate: 'derives_from'},  # CellLine-->Disease
        'CLEAVESTO': {predicate: 'affects', QUALIFIED_ASPECT: 'cleavage'},  # Protein-->Protein
        'RESPONSE_TO': {predicate: 'affects_response_to'},  # Gene-->Compound
        'SAME': {predicate: 'same_as'},  # CellLine-->CellLine TODO: maybe check this one.. (look at example edges)
        'RESISTANT_TO': {predicate: 'associated_with_resistance_to'},  # Gene-->Compound
        'DECREASEDIN': {predicate: 'has_decreased_amount', flip: True},  # Protein-->Disease
        'ADVRESPONSE_TO': {predicate: 'contraindicated_in', flip: True},  # Gene-->Compound
        'REDUCES_SEN': {predicate: 'decreases_response_to'}  # Gene-->Compound
    }
    type_mapping = type_map[edge_type] if edge_type in type_map else type_map[core_edge_type]
    subject_id = original_object_id if type_mapping.get(flip) else original_subject_id
    object_id = original_subject_id if type_mapping.get(flip) else original_object_id
    
    return (f"{BIOLINK_PREFIX}:{type_mapping[predicate]}", subject_id, object_id,
            type_mapping.get(QUALIFIED_PREDICATE), type_mapping.get(QUALIFIED_DIRECTION), type_mapping.get(QUALIFIED_ASPECT))


def get_all_sources(item: dict) -> Tuple[str, List[str]]:
    properties = item['properties']
    sources = []
    if 'source' in properties and properties['source']:
        sources.append(properties['source'])
    if 'sources' in properties and properties['sources']:
        if isinstance(properties['sources'], list):
            sources.extend(properties['sources'])
        else:
            sources.append(str(properties['sources']))

    primary_source = sources[0] if sources else 'unknown'
    secondary_sources = sources[1:] if len(sources) > 1 else []

    return primary_source, secondary_sources
