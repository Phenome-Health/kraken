"""
SPOKE harmonizer - converts SPOKE format to unified Biolink schema
"""

import json
from pathlib import Path
import sys
from typing import List, Optional, Tuple
import jsonlines
import logging

from ..utils.general import load_biolink_file, create_node, create_edge
from ..utils.constants import *
from ..utils.kg_io import stream_mixed_jsonl
from .spoke_id_utils import SpokeIDNormalizer


def harmonize_spoke(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    """Harmonize SPOKE mixed JSONL to unified Biolink schema using streaming"""
    logging.info(f"Harmonizing SPOKE: {input_file} -> {nodes_output}, {edges_output}")

    # Initialize identifier normalizer and other biolink info needed for conversion
    id_norm = SpokeIDNormalizer(biolink_version=biolink_version)
    infores_url = "https://raw.githubusercontent.com/biolink/information-resource-registry/refs/heads/main/infores_catalog.yaml"
    infores_info = load_biolink_file(infores_url, biolink_version)
    klat_map = {item['id']: {'knowledge_level': item.get('knowledge_level', 'not_provided'),
                             'agent_type': item.get('agent_type', 'not_provided')}
                for item in infores_info['information_resources']}
    
    # Keep track of normalized node IDs for edge mapping
    spoke_to_normalized_id = {}

    node_count = 0
    edge_count = 0
    with jsonlines.open(nodes_output, 'w') as nodes_writer, \
         jsonlines.open(edges_output, 'w') as edges_writer:
        
        for item in stream_mixed_jsonl(input_file):
            item_type = item.get('type')
            
            if item_type == 'node':
                harmonized_node = harmonize_spoke_node(item, id_norm)

                if harmonized_node:  # Occasionally we skip nodes (if invalid identifier, etc...)
                    # Store mapping for edge processing
                    spoke_to_normalized_id[item['id']] = harmonized_node[ID]
                    
                    nodes_writer.write(harmonized_node)
                    node_count += 1

            elif item_type == 'relationship':
                harmonized_edge = harmonize_spoke_edge(item, spoke_to_normalized_id, klat_map)
                if harmonized_edge:
                    edges_writer.write(harmonized_edge)
                    edge_count += 1

    logging.info(f"SPOKE harmonization complete: {node_count} nodes, {edge_count} edges")


def harmonize_spoke_node(node_item: dict, id_norm: SpokeIDNormalizer) -> Optional[dict]:
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

    # Handle special case in SPOKE nodes
    if primary_source == 'Complex Portal' and node_type == 'Complex':
        original_identifier = properties['complex_portal']  # The 'identifier' for these is meaningless; true standard id is tucked away here

    # Normalize the identifier
    normalized_id, iri = id_norm.normalize_spoke_identifier(node_type, primary_source, original_identifier, properties)

    if normalized_id == KNOWN_INVALID:
        # logging.warning(f"Skipping node as curie extraction failed (known failure). {node_item}")
        return None
    elif normalized_id:
        # Extract additional equivalent identifiers from properties
        additional_equivalent_ids = id_norm.extract_equivalent_identifiers(node_type, properties)
        all_equivalent_ids = list(set([normalized_id] + additional_equivalent_ids))

        harmonized_node = create_node(curie=normalized_id,
                                      categories=map_spoke_labels_to_biolink(labels, primary_source, normalized_id),
                                      provided_by=[SPOKE_INFORES],
                                      equivalent_ids=all_equivalent_ids,
                                      name=properties.get('name'),
                                      synonyms=[properties['name']] if properties.get('name') else None,
                                      iri=iri,
                                      attributes={SPOKE_INFORES: {'id': node_item['id']}})

        return harmonized_node
    else:
        logging.error(f"Failed to convert SPOKE 'identifier' to a proper curie. {node_item}")
        sys.exit(1)


def harmonize_spoke_edge(edge_item: dict, spoke_to_normalized_id: dict, klat_map: dict) -> Optional[dict]:
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

    normalized_primary_ks = normalize_source(primary_source, edge_item)

    harmonized_edge = create_edge(subject_id=normalized_subject_id,
                                  object_id=normalized_object_id,
                                  predicate=predicate,
                                  primary_ks=normalized_primary_ks,
                                  knowledge_level=klat_map.get(normalized_primary_ks, dict()).get('knowledge_level', 'not_provided'),
                                  agent_type=klat_map.get(normalized_primary_ks, dict()).get('agent_type', 'not_provided'),
                                  aggregator_ks=SPOKE_INFORES,
                                  supporting_sources=list({normalize_source(s, edge_item) for s in secondary_sources}),
                                  qualified_predicate=qual_predicate,
                                  qualified_direction=qual_direction,
                                  qualified_aspect=qual_aspect,
                                  attributes={SPOKE_INFORES: {'id': edge_item['id']}})

    return harmonized_edge


def map_spoke_labels_to_biolink(labels: List[str], source: str, standardized_id: str) -> List[str]:
    """Map SPOKE node labels to Biolink categories"""

    # Simple mapping - extend as needed
    label_mapping = {
        'Compound': 'ChemicalEntity',  # This more accurate than SmallMolecule?
        'Variant': 'SequenceVariant',
        'Organism--ncbi-taxonomy': 'OrganismTaxon',  # Need to consider source for these
        'Organism--BV-BRC': 'OrganismalEntity',  # Need to consider source for these
        'Protein': 'Polypeptide' if 'PRO_' in standardized_id else 'Protein',  # Some are really protein features
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


def normalize_source(spoke_source: str, spoke_edge: dict) -> str:
   spoke_source_cleaned = spoke_source.lower().replace(' ', '')
   mappings = {
       'ahrqsdohdatabase': 'ahrq-sdoh',
       'bgee': f'{INFORES_PREFIX}:bgee',
       'mirdb': 'mirdb',
       'opentargets': f'{INFORES_PREFIX}:open-targets',
       'string': f'{INFORES_PREFIX}:string',
       'stitch': f'{INFORES_PREFIX}:stitch',
       'humanproteinatlas': f'{INFORES_PREFIX}:hpa',
       'ncbigene2go': 'ncbi-gene2go',
       'bindingdb': f'{INFORES_PREFIX}:bindingdb',
       'cmap/lincscompound(trt_cp)': f'{INFORES_PREFIX}:lincs',  # Check?
       'bv-brc': 'bv-brc',
       'uniprot': f'{INFORES_PREFIX}:uniprot',
       'intact': f'{INFORES_PREFIX}:intact',
       'cmap/lincsknockdown(trt_xrt,trt_sh)': f'{INFORES_PREFIX}:lincs',  # Check?
       'clinvar': f'{INFORES_PREFIX}:clinvar',
       'interpro': f'{INFORES_PREFIX}:interpro',
       'chebi': f'{INFORES_PREFIX}:chebi',
       'hpo': f'{INFORES_PREFIX}:hpo',
       'places': 'cdc-places',
       'cancercelllineencyclopedia': 'ccle',
       'nhanes': 'nhanes',
       'countyhealthrankings': 'chr-r',
       'loinc': f'{INFORES_PREFIX}:loinc',
       'ncbipubmed': f'{INFORES_PREFIX}:pubmed',
       'sider4.1': f'{INFORES_PREFIX}:sider',
       'unitedstateszipcode_database': 'us-zipcode',
       'protcid': 'prot-cid',
       'metacyc': f'{INFORES_PREFIX}:metacyc',
       'cmap/lincsoverexpression(trt_oe)': f'{INFORES_PREFIX}:lincs',  # Check?
       'wikipathways': f'{INFORES_PREFIX}:wikipathways',
       'ctkp': f'{INFORES_PREFIX}:multiomics-clinicaltrials',
       'diseases': f'{INFORES_PREFIX}:diseases',
       'chembl': f'{INFORES_PREFIX}:chembl',
       'kegg': f'{INFORES_PREFIX}:kegg',
       'drugcentral': f'{INFORES_PREFIX}:drugcentral',
       'celltaxonomy': 'cell-taxonomy',
       'bioplex(pharos)': f'{INFORES_PREFIX}:pharos',
       'tflink': 'tf-link',
       'tri': 'tri',
       'ucmr5': 'ucmr',
       'superfund': 'epa-superfund',
       'reactome': f'{INFORES_PREFIX}:reactome',
       'uberon': f'{INFORES_PREFIX}:uberon',
       'gwascatalog': f'{INFORES_PREFIX}:gwas-catalog',
       'ensemblhg38': f'{INFORES_PREFIX}:ensembl-gene',
       'cmap/lincsligand(trt_lig)': f'{INFORES_PREFIX}:lincs',
       'ucmr4': 'ucmr',
       'diseaseontology': f'{INFORES_PREFIX}:disease-ontology',
       'cancerrx': f'{INFORES_PREFIX}:gdsc',
       'ncbi-taxonomy': f'{INFORES_PREFIX}:ncbi-taxonomy',
       'worldhealthorganization': 'who',
       'worldhealhorganization': 'who',  # Some SPOKE edges have this typo
       'celllineontology': 'clo',
       'pfam': f'{INFORES_PREFIX}:pfam',
       'gwas': f'{INFORES_PREFIX}:gwas-catalog',
       'pharmvar': f'pharmvar',
       'explorenz': 'explor-enz',
       'complexportal': f'{INFORES_PREFIX}:complex-portal',
       '2020nationalemissionsinventory(nei)data': 'epa-nei',
       'cellontology': f'{INFORES_PREFIX}:cl',
       'whoambientairqualitydatabase': 'who-air-quality',
       'geonames': 'geonames',
       'cdc/atsdrsocialvulnerabilityindex': 'svi',
       'pharmgkb': f'{INFORES_PREFIX}:pharmgkb',
       'cellosaurus': 'cellosaurus',
       'tcdb': 'tcdb',
       'mirbase': f'{INFORES_PREFIX}:mirbase',
       'eqtlcatalogue': 'eqtl-catalogue',
       'cellphonedb': 'cellphone-db',
       'civic': f'{INFORES_PREFIX}:civic',
       'nationalcenterforhealthstatistics.u.s.censusbureau,householdpulsesurvey,2024.lackofsocialconnection4.2': 'pulse-survey',
       'pathophenodb': f'{INFORES_PREFIX}:path-pheno-db',
       'airqualitystatisticsreport': 'epa-air-quality-stats',
       'https://github.com/hadlock_lab/recover/': 'hadlock-recover',
       'https://github.com/hadlock_lab/incov/': 'hadlock-incov',
       'unknown': SPOKE_INFORES  # If SPOKE doesn't give a source for the edge, just list SPOKE as the source.. (better than nothing)
   }
   if spoke_source_cleaned in mappings:
       return mappings[spoke_source_cleaned]
   else:
       logging.error(f'Encountered an unmapped edge source: {spoke_source} ({spoke_source_cleaned}). Edge: {spoke_edge}')
       sys.exit(1)


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
