import jsonlines

from kraken.orchestrator import KrakenBuildOrchestrator
from kraken.utils.constants import EQUIVALENT_IDS, ID, PROJECT_ROOT
from kraken.utils.kg_io import stream_nodes_from_jsonl

target_kg_name = "spoke"
target_prefixes = {"ENVO"}

orchestrator = KrakenBuildOrchestrator()
nodes_path, edges_path = orchestrator.config.all_harmonized_paths_resolved(target_kg_name)


subgraph_dir = PROJECT_ROOT / "artifacts" / "subgraphs"
subgraph_dir.mkdir(parents=True, exist_ok=True)
subgraph_nodes_path = subgraph_dir / "spoke_envo_nodes.jsonl"


print("Extracting nodes..")
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
