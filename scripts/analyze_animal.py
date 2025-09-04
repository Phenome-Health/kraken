

from collections import defaultdict
import json
import os
import sys
import jsonlines

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl
from utils.logging_config import setup_logging


setup_logging()

nodes_path = '/Users/amyglen/phenome-kg/artifacts/integrated/kraken_nodes.jsonl'
categories_map = {node['id']: node['categories'] for node in stream_nodes_from_jsonl(nodes_path)}

taxon_meta = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {'count': 0, 'example': None})))

for edge in stream_edges_from_jsonl('/Users/amyglen/phenome-kg/artifacts/integrated/kraken_edges.jsonl'):
    subj_categories = categories_map[edge['subject']]
    obj_categories = categories_map[edge['object']]
    predicate = edge['predicate']
    primary_ks = edge['primary_knowledge_source']
    aggregator_ks = edge.get('aggregator_knowledge_source', 'none')

    if predicate != 'biolink:close_match':
        for subj_category in subj_categories:
            for obj_category in obj_categories:
                if subj_category == 'biolink:OrganismTaxon' or obj_category == 'biolink:OrganismTaxon':
                    triple = f"{subj_category}--{predicate}--{obj_category}"
                    taxon_meta[triple][primary_ks][aggregator_ks]['count'] += 1
                    
                    # Save an example edge if we haven't seen one yet
                    if not taxon_meta[triple][primary_ks][aggregator_ks]['example']:
                        taxon_meta[triple][primary_ks][aggregator_ks]['example'] = edge


# Dump the counts and example edges
taxon_metainfo_path = '/Users/amyglen/phenome-kg/scripts/animal_meta.json'
with open(taxon_metainfo_path, 'w+') as output_file:
    json.dump(taxon_meta, output_file, indent=2)

