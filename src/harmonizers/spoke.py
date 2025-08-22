"""
SPOKE harmonizer - converts SPOKE format to unified Biolink schema
"""

from pathlib import Path
from typing import List, Optional, Tuple
import jsonlines
import logging
from ..utils.metagraph import generate_metagraph_for_source
from .simple_identifier_utils import SimpleIdentifierNormalizer


def harmonize_spoke(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str, rules: dict):
    """Harmonize SPOKE mixed JSONL to unified Biolink schema using streaming"""
    logging.info(f"Harmonizing SPOKE: {input_file} -> {nodes_output}, {edges_output}")

    node_count = 0
    edge_count = 0
    
    # Initialize identifier normalizer
    id_norm = SimpleIdentifierNormalizer(biolink_version=biolink_version)
    
    # Keep track of normalized node IDs for edge mapping
    spoke_to_normalized_id = {}

    with jsonlines.open(input_file, 'r') as reader, \
         jsonlines.open(nodes_output, 'w') as nodes_writer, \
         jsonlines.open(edges_output, 'w') as edges_writer:
        
        for line_num, item in enumerate(reader, 1):
            item_type = item.get('type')
            
            if item_type == 'node':
                harmonized_node = harmonize_spoke_node(item, id_norm)
                
                # Store mapping for edge processing
                spoke_to_normalized_id[item['id']] = harmonized_node['id']
                
                nodes_writer.write(harmonized_node)
                node_count += 1
                
                if node_count % 500000 == 0:
                    logging.info(f"Processed {node_count} SPOKE nodes")
            
            elif item_type == 'relationship':
                harmonized_edge = harmonize_spoke_edge(item, spoke_to_normalized_id)
                edges_writer.write(harmonized_edge)
                edge_count += 1
                
                if edge_count % 1000000 == 0:
                    logging.info(f"Processed {edge_count} SPOKE edges")
    
    logging.info(f"SPOKE harmonization complete: {node_count} nodes, {edge_count} edges")
    
    # Generate metagraph for harmonized output
    if rules.get('generate_metagraph', True):
        # Store metagraphs in artifacts/metagraphs/harmonized/source_name/
        artifacts_root = Path("artifacts")
        metagraph_dir = artifacts_root / "metagraphs" / "harmonized" / "spoke"
        
        metagraph_config = rules.get('metagraph_config', {
            'generate_summaries': True,
            'generate_cytoscape': True,
            'generate_html_viewer': True,
            'cytoscape_thresholds': [1, 5, 10]
        })
        
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, "spoke", metagraph_config)
        logging.info("SPOKE metagraph generated")


def harmonize_spoke_node(node_item: dict, id_norm: SimpleIdentifierNormalizer) -> dict:
    """Harmonize a single SPOKE node"""
    properties = node_item.get('properties', {})
    labels = node_item.get('labels', [])
    if not labels:
        raise ValueError(f"SPOKE node is missing labels: {node_item}")
    
    # Get node type and source for identifier normalization
    node_type = labels[0]  # Primary label
    
    # Get source(s) - handle both 'source' and 'sources' 
    primary_source, secondary_sources = get_all_sources(node_item)
    
    # Normalize the identifier
    original_identifier = properties.get('identifier', node_item['id'])
    normalized_id = id_norm.normalize_spoke_identifier(node_type, primary_source, original_identifier)
    
    # Extract additional equivalent identifiers from properties
    additional_equivalent_ids = id_norm.extract_equivalent_identifiers(node_type, properties)
    all_equivalent_ids = list(set([normalized_id] + additional_equivalent_ids))
    
    harmonized_node = {
        'id': normalized_id,
        'categories': map_spoke_labels_to_biolink(labels),
        'provided_by': ['infores:spoke'],
        'equivalent_ids': all_equivalent_ids,
        'spoke_node': node_item
    }
    if properties.get('name'):
        harmonized_node['name'] = properties['name']
        harmonized_node['synonyms'] = [harmonized_node['name']]  # TODO: down the road see about extracting SPOKE synonyms
    
    return harmonized_node


def harmonize_spoke_edge(edge_item: dict, spoke_to_normalized_id: dict) -> dict:
    """Harmonize a single SPOKE edge"""
    edge_type = edge_item.get('label')
    if not edge_type:
        raise ValueError(f"SPOKE edge is missing type: {edge_item}")

    spoke_subject_id = edge_item['start']['id']
    spoke_object_id = edge_item['end']['id']
    predicate, subject_id, object_id, qual_predicate, qual_direction, qual_aspect = map_spoke_edge_type_to_biolink(edge_type, spoke_subject_id, spoke_object_id)
    
    # Get source(s) - handle both 'source' and 'sources' 
    primary_source, secondary_sources = get_all_sources(edge_item)

    # Map SPOKE internal IDs to normalized CURIEs
    normalized_subject_id = spoke_to_normalized_id.get(subject_id, subject_id)
    normalized_object_id = spoke_to_normalized_id.get(object_id, object_id)
    
    # Remove the full start/end node objects (replace with their SPOKE IDs instead - saves a lot of space)
    edge_item['start'] = subject_id
    edge_item['end'] = object_id

    harmonized_edge = {
        'subject': normalized_subject_id,
        'object': normalized_object_id,
        'predicate': predicate,
        'primary_knowledge_source': primary_source,  # TODO: Convert to infores curies where possible...
        'aggregator_knowledge_source': 'infores:spoke',
        'spoke_edge': edge_item
    }
    # Tack on any additional sources
    if secondary_sources:
        harmonized_edge['supporting_data_sources'] = secondary_sources  # TODO: Convert to infores curies where possible...
    # Tack on any qualifiers
    if qual_predicate:
        harmonized_edge['qualified_predicate'] = qual_predicate
    if qual_direction:
        harmonized_edge['qualified_direction'] = qual_direction
    if qual_aspect:
        harmonized_edge['qualified_aspect'] = qual_aspect
    
    return harmonized_edge


def map_spoke_labels_to_biolink(labels: List[str]) -> List[str]:
    """Map SPOKE node labels to Biolink categories"""

    # Simple mapping - extend as needed
    label_mapping = {
        'Compound': 'ChemicalEntity',  # This more accurate than SmallMolecule?
        'Variant': 'SequenceVariant',
        'Organism': 'OrganismTaxon',
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
        'EC': 'EC',  # TODO: Collapse these nodes into annotations on Protein nodes (see issue)
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
    
    return [f"biolink:{label_mapping[label]}" for label in labels]


def map_spoke_edge_type_to_biolink(edge_type: str, 
                                   original_subject_id: str, 
                                   original_object_id: str) -> Tuple[str, str, str, Optional[str], Optional[str], Optional[str]]:
    """Map SPOKE edge types to Biolink predicates; flips edges as necessary to use canonical predicates."""
    core_edge_type = '_'.join(edge_type.split('_')[:-1])  # Gets rid of suffix indicating node categories, like _GiP

    # Simple mapping - extend as needed
    predicate = 'type'
    flip = 'flip'
    qual_predicate = 'qualified_predicate'
    qual_direction = 'qualified_direction'
    qual_aspect = 'qualified_aspect'
    type_map = {
        'PREVALENCEIN': {predicate: 'associated_with'},  # SPOKE only uses for SocioeconomicExposure-->GeographicLocation edges TODO: would occurs_in be better?
        'INTERACTS': {predicate: 'interacts_with'},
        'MAPS': {predicate: 'is_sequence_variant_of'},  # SPOKE only uses for SequenceVariant-->Gene edges
        'TARGETS': {predicate: 'regulates'},  # SPOKE only uses for MiRNA-->Gene edges
        'EXPRESSES': {predicate: 'expressed_in', flip: True},  # SPOKE only uses for Anatomy-->Gene edges
        'DOWNREGULATES': {predicate: 'regulates', qual_direction: 'downregulated'},  # Anatomy-->Gene, Compound-->Gene, Gene-->Gene, Variant-->Gene
        'BINDS': {predicate: 'binds'},  # Compound-->Protein/ProteinDomain
        'UPREGULATES': {predicate: 'regulates', qual_direction: 'upregulated'},  # Protein-->Gene
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
        'RESPONDS': {predicate: 'affects', flip: True},  # Organism-->Compoound
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
        'MARKER': {predicate: 'biomarker_for'},  # Gene-->Disease
        'AFFECTS': {predicate: 'affects'},
        'INCREASEDIN': {predicate: 'has_increased_amount', flip: True},  # Protein-->Disease
        'CONTRAINDICATES': {predicate: 'contraindicated_in'},  # Compound-->Disease
        'MORTALITY': {predicate: 'occurs_in'},  # Disease-->GeographicLocation
        'MEMBEROF': {predicate: 'has_member', flip: True},  # ProteinDomain-->ProteinFamily
        'INCLUDES': {predicate: 'has_member'},  # PharmacologicClass-->Compound
        'AFFECT': {predicate: 'affects_response_to'},  # Variant-->Compound
        'TRANSPORTS': {predicate: 'affects', qual_aspect: 'transport'},  # Protein-->Compound
        'DERIVES': {predicate: 'derives_from'},  # CellLine-->Disease
        'CLEAVESTO': {predicate: 'affects', qual_aspect: 'cleavage'},  # Protein-->Protein
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
    
    return (f"biolink:{type_mapping[predicate]}", subject_id, object_id, 
            type_mapping.get(qual_predicate), type_mapping.get(qual_direction), type_mapping.get(qual_aspect))


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
