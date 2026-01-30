import jsonlines

from kraken.metagraph import generate_metagraph_for_source
from kraken.orchestrator import KrakenBuildOrchestrator
from kraken.utils.constants import ID, OBJECT, PRIMARY_KS, PROJECT_ROOT, SUBJECT, SUPPORTING_SOURCES
from kraken.utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl

target_kg_name = "kg2"
target_knowledge_sources = {"infores:icd10pcs-umls"}

orchestrator = KrakenBuildOrchestrator()
nodes_path, edges_path = orchestrator.config.all_harmonized_paths_resolved[target_kg_name]

subgraph_name = "kg2_icd10pcsumls"
subgraph_dir = PROJECT_ROOT / "artifacts" / "subgraphs"
subgraph_dir.mkdir(parents=True, exist_ok=True)
subgraph_nodes_path = subgraph_dir / f"{subgraph_name}_nodes.jsonl"
subgraph_edges_path = subgraph_dir / f"{subgraph_name}_edges.jsonl"


# First extract the edges, recording which nodes they involve
print("Extracting edges..")
subgraph_node_ids = set()
num_subgraph_edges = 0
with jsonlines.open(subgraph_edges_path, "w") as subgraph_writer:
    for edge in stream_edges_from_jsonl(edges_path):
        primary_ks = edge[PRIMARY_KS]
        supporting_sources = set(edge.get(SUPPORTING_SOURCES, []))
        from_target_source = primary_ks in target_knowledge_sources or supporting_sources.intersection(
            target_knowledge_sources
        )
        if (target_knowledge_sources and from_target_source) or not target_knowledge_sources:
            subgraph_node_ids.add(edge[SUBJECT])
            subgraph_node_ids.add(edge[OBJECT])
            subgraph_writer.write(edge)
            num_subgraph_edges += 1
print(f"Finished extracting {num_subgraph_edges} subgraph edges.")

# Then go through and extract the nodes the subgraph edges involve
print("Extracting nodes..")
num_subgraph_nodes = 0
with jsonlines.open(subgraph_nodes_path, "w") as subgraph_writer:
    for node in stream_nodes_from_jsonl(nodes_path):
        if node[ID] in subgraph_node_ids:
            subgraph_writer.write(node)
            num_subgraph_nodes += 1
print(f"Finished extracting {num_subgraph_nodes} subgraph nodes.")

generate_metagraph_for_source(subgraph_nodes_path, subgraph_edges_path, subgraph_dir, subgraph_name)
