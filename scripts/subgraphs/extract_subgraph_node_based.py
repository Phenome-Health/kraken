import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from utils.constants import PRIMARY_KS, SUPPORTING_SOURCES, SUBJECT, OBJECT, ID, EQUIVALENT_IDS
from utils.kg_io import get_harmonized_file_paths, PROJECT_ROOT, stream_edges_from_jsonl, stream_nodes_from_jsonl
from utils.metagraph import generate_metagraph_for_source

import jsonlines

target_kg_name = "spoke"
target_prefixes = {"ENVO"}

nodes_path, edges_path = get_harmonized_file_paths(target_kg_name)


subgraph_dir = PROJECT_ROOT / "artifacts" / "subgraphs"
subgraph_dir.mkdir(parents=True, exist_ok=True)
subgraph_nodes_path = subgraph_dir / "spoke_envo_nodes.jsonl"


print(f"Extracting nodes..")
num_subgraph_nodes = 0
target_node_ids = set()
with jsonlines.open(subgraph_nodes_path, "w") as subgraph_writer:
    for node in stream_nodes_from_jsonl(nodes_path):
        equiv_ids = node[EQUIVALENT_IDS]
        prefixes = {equiv_id.split(":")[0] for equiv_id in equiv_ids}
        if prefixes.intersection(target_prefixes):
            target_node_ids.add(node[ID])
            subgraph_writer.write(node)
            num_subgraph_nodes += 1
print(f"Finished extracting {num_subgraph_nodes} subgraph nodes.")
