import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from utils.constants import PRIMARY_KS, SUPPORTING_SOURCES, SUBJECT, OBJECT, ID
from utils.kg_io import get_harmonized_file_paths, PROJECT_ROOT, stream_edges_from_jsonl, stream_nodes_from_jsonl
from utils.metagraph import generate_metagraph_for_source

import jsonlines

target_kg_name = "spoke"
target_knowledge_sources = {"us-zipcode", "ucmr", "tri", "epa-superfund", "who", "epa-nei",
                            "who-air-quality", "geonames", "epa-air-quality-stats"}

nodes_path, edges_path = get_harmonized_file_paths(target_kg_name)


subgraph_dir = PROJECT_ROOT / "artifacts" / "subgraphs"
subgraph_dir.mkdir(parents=True, exist_ok=True)
subgraph_nodes_path = subgraph_dir / "spoke_env_geo_nodes.jsonl"
subgraph_edges_path = subgraph_dir / "spoke_env_geo_edges.jsonl"


# First extract the edges, recording which nodes they involve
print(f"Extracting edges..")
subgraph_node_ids = set()
num_subgraph_edges = 0
with jsonlines.open(subgraph_edges_path, "w") as subgraph_writer:
    for edge in stream_edges_from_jsonl(edges_path):
        primary_ks = edge[PRIMARY_KS]
        supporting_sources = set(edge.get(SUPPORTING_SOURCES, []))
        from_target_source = primary_ks in target_knowledge_sources or supporting_sources.intersection(target_knowledge_sources)
        if (target_knowledge_sources and from_target_source) or not target_knowledge_sources:
            subgraph_node_ids.add(edge[SUBJECT])
            subgraph_node_ids.add(edge[OBJECT])
            subgraph_writer.write(edge)
            num_subgraph_edges += 1
print(f"Finished extracting {num_subgraph_edges} subgraph edges.")

# Then go through and extract the nodes the subgraph edges involve
print(f"Extracting nodes..")
num_subgraph_nodes = 0
with jsonlines.open(subgraph_nodes_path, "w") as subgraph_writer:
    for node in stream_nodes_from_jsonl(nodes_path):
        if node[ID] in subgraph_node_ids:
            subgraph_writer.write(node)
            num_subgraph_nodes += 1
print(f"Finished extracting {num_subgraph_nodes} subgraph nodes.")

generate_metagraph_for_source(subgraph_nodes_path, subgraph_edges_path, subgraph_dir, "spoke_env_geo_subgraph")
