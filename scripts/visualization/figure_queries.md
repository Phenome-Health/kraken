# Figure queries and commands (KRAKEN v2.3.0)

Record of exactly how each manuscript figure was generated. Update this whenever a figure is finalized.

UI screenshots were captured from a local copy of the web interface (`http://localhost:5173`) pointed at the
production backend (`https://kestrel.krakenkg.com`). The links below use the public site; the part after `#` is
identical on localhost. Screenshot settings are the page width and height in CSS pixels and the pixel density.

## Figure 1: example node and edge

- Node panel: look up Creatinine by id, which opens its Node Details panel directly.
  `https://app.krakenkg.com/#/search?ids=RM:0028430`
- Edge panel: run the query, click the "correlated with" edge in the first result, then "Show attributes".
  `https://app.krakenkg.com/#/path?start=RM:0028430&end=MONDO:0005575&hops=1-2&ranking=established&traversal=deep_dive`
- Panels captured at full height, 4x pixel density, then composed:
  `uv run --group viz python scripts/visualization/plot_node_edge_panels.py scripts/visualization/ex_node_full.png scripts/visualization/ex_edge_full.png -o fig1.pdf --dpi 1289 --gap 0.1 --pad 0.03`

## Figure 2: multi-hop results for cholesterol

- Click the end node of the first result to open Node Details.
  `https://app.krakenkg.com/#/path?start=CHEBI:16113&end_category=Pathway,BiologicalProcessOrActivity&hops=1-2&ranking=established&traversal=deep_dive&intermediate_node_category.in=biolink:Protein`
- Screenshot: 1300 x 900, 4x.

## Figure 3: query builder

- Same link as Figure 2, then click "Edit Query". Cropped to the builder (start node through Execute Query).
- Screenshot: 950 wide, 5x.

## Figure 4: node search

- Click the first result (atrial fibrillation, `MONDO:0004981`).
  `https://app.krakenkg.com/#/search?text=atrial%20fibrillation`
- Screenshot: 1500 x 880, 4x.

## Figure 5: chord diagram

`uv run --group viz python scripts/visualization/plot_metagraph_chord.py /Volumes/AmySSD/kraken-data/artifacts/metagraphs/kraken_metagraph_2.3.0.json --dpi 500`

Uses `chord_sqrt.png`.

## Figure 6: source network

`uv run --group viz python scripts/visualization/plot_source_network.py /Volumes/AmySSD/kraken-data/artifacts/integrated/kraken_edges_2.3.0.jsonl --nodes /Volumes/AmySSD/kraken-data/artifacts/integrated/kraken_nodes_2.3.0.jsonl --cache scan_2.3.0.pkl --exclude "ubergraph,robokop-kg,rtx-kg2,translator-kg-open,sri-node-normalizer" --dpi 400 --legend-loc "lower left" -o fig6.pdf`

Result: 240 links between 98 sources.

## Figure 7: knowledge level by agent type

`uv run --group viz python scripts/visualization/plot_klat.py /Volumes/AmySSD/kraken-data/artifacts/metagraphs/kraken_metagraph_2.3.0.json -o fig7.pdf --dpi 600`

Uses the log-scaled heatmap (`fig7_heatmap_log`, copied to `fig7`).

## Figure 8: enrichment of the 88 signature proteins

- Default permutations (1,000). Row "inflammatory response" selected. 400 targets at FDR < 0.05 (seen with `limit=500`).
  `https://app.krakenkg.com/#/enrichment?nodes=NCBIGene:147945,NCBIGene:3952,NCBIGene:133,NCBIGene:3484,NCBIGene:2167,NCBIGene:177,NCBIGene:1435,NCBIGene:3603,NCBIGene:92737,NCBIGene:5653,NCBIGene:56477,NCBIGene:10673,NCBIGene:6441,NCBIGene:54,NCBIGene:181,NCBIGene:1956,NCBIGene:117156,NCBIGene:4319,NCBIGene:4973,NCBIGene:4151,NCBIGene:5806,NCBIGene:8788,NCBIGene:2277,NCBIGene:9965,NCBIGene:3627,NCBIGene:3977,NCBIGene:1476,NCBIGene:8797,NCBIGene:920,NCBIGene:3558,NCBIGene:1473,NCBIGene:6401,NCBIGene:58191,NCBIGene:923,NCBIGene:27189,NCBIGene:4854,NCBIGene:8685,NCBIGene:2658,NCBIGene:6373,NCBIGene:7037,NCBIGene:1357,NCBIGene:11093,NCBIGene:3958,NCBIGene:3596,NCBIGene:50604,NCBIGene:763,NCBIGene:5328,NCBIGene:7039,NCBIGene:4353,NCBIGene:9048,NCBIGene:6504,NCBIGene:3566,NCBIGene:4313,NCBIGene:26291,NCBIGene:83886,NCBIGene:7079,NCBIGene:56729,NCBIGene:127482540,NCBIGene:142,NCBIGene:355,NCBIGene:6370,NCBIGene:3593,NCBIGene:3976,NCBIGene:145264,NCBIGene:3339,NCBIGene:1522,NCBIGene:8600,NCBIGene:3425,NCBIGene:8743,NCBIGene:2250,NCBIGene:100,NCBIGene:9332,NCBIGene:53832,NCBIGene:5327,NCBIGene:3605,NCBIGene:90865,NCBIGene:23584,NCBIGene:126014,NCBIGene:57823,NCBIGene:8993,NCBIGene:133395145,NCBIGene:4803,NCBIGene:6347,NCBIGene:654,NCBIGene:23765,NCBIGene:51744,NCBIGene:6359,NCBIGene:3559&target=Pathway,BiologicalProcessOrActivity&predicate=has_participant,participates_in&limit=30`
- Screenshot: 1300 x 1000, 4x.

## Figure 9: immune module subgraph (41 seed proteins)

- Row "cytokine-mediated signaling pathway" selected in the Intermediate Nodes table. 89 nodes, 283 collapsed (735 total) edges.
  `https://app.krakenkg.com/#/subgraph?nodes=NCBIGene:9048,NCBIGene:6359,NCBIGene:6347,NCBIGene:6370,NCBIGene:56477,NCBIGene:9332,NCBIGene:51744,NCBIGene:920,NCBIGene:923,NCBIGene:1401,NCBIGene:1435,NCBIGene:3627,NCBIGene:6373,NCBIGene:58191,NCBIGene:355,NCBIGene:3593,NCBIGene:3596,NCBIGene:3603,NCBIGene:3605,NCBIGene:27189,NCBIGene:23765,NCBIGene:3558,NCBIGene:50604,NCBIGene:53832,NCBIGene:3559,NCBIGene:90865,NCBIGene:3566,NCBIGene:3958,NCBIGene:3976,NCBIGene:3977,NCBIGene:8685,NCBIGene:4353,NCBIGene:4973,NCBIGene:8993,NCBIGene:5806,NCBIGene:6504,NCBIGene:57823,NCBIGene:8797,NCBIGene:8743,NCBIGene:8600,NCBIGene:10673&hops=2&limit=50&ranking=deep_dive&traversal=deep_dive&intermediate_node_category.in=biolink:Protein,biolink:Gene,biolink:Pathway&prefix.in=NCBIGene,HGNC,GO,REACT&agent_type.not_in=text_mining_agent`
- Screenshot: page 1500 x 1000 at 4x, cropped to the graph and table (top navigation bar and right-hand details panel excluded).

## Figure 10: asthma PRS proteins to compounds (three panels)

Query-graph queries with the protein, drug and asthma pinned, text-mined edges excluded on every edge and expansion,
and equivalent terms allowed on the drug node (1 hop; categories listed so equivalents of another type still match).
Each panel is the result's graph only (page 1400 x 1000 at 4x, cropped to the graph); the `layout` parameter places the protein node closer to the drug node.

- Panel A, IL5RA, mepolizumab, asthma (no expansion on asthma):
  `https://app.krakenkg.com/#/querygraph?n0=ATC:R03DX09&n0.category=Drug,Protein&n0.expand=exact_match,close_match,same_as&n0.expand.hops=1&n0.expand.agent_type.not_in=text_mining_agent&n1=HGNC:6017&n2=MONDO:0004979&e0=n0-n1&e0.agent_type.not_in=text_mining_agent&e1=n0-n2&e1.predicate=treats_or_applied_or_studied_to_treat&e1.agent_type.not_in=text_mining_agent&ranking=established&traversal=deep_dive&layout=n0.757.91,n1.557.91,n2.1057.91`
- Panel B, MMP10, marimastat, asthma (asthma expanded 2 hops to subtypes, parent types and equivalent terms; 1 path per pair):
  `https://app.krakenkg.com/#/querygraph?n0=RM:0224430&n0.category=Drug,SmallMolecule&n0.expand=exact_match,close_match,same_as&n0.expand.hops=1&n0.expand.agent_type.not_in=text_mining_agent&n1=HGNC:7156&n2=MONDO:0004979&n2.expand=subclass_of,superclass_of,exact_match,close_match,same_as&n2.expand.hops=2&n2.expand.agent_type.not_in=text_mining_agent&e0=n0-n1&e0.agent_type.not_in=text_mining_agent&e1=n0-n2&e1.predicate=treats_or_applied_or_studied_to_treat&e1.agent_type.not_in=text_mining_agent&paths_per_pair=1&ranking=established&traversal=deep_dive&layout=n0.757.91,n1.557.91,n2.1057.91`
- Panel C, TNFRSF8, brentuximab vedotin (node `ATC:L01FX05`, named "Adcetris"), asthma (same asthma expansion; 2 paths per pair):
  `https://app.krakenkg.com/#/querygraph?n0=ATC:L01FX05&n0.category=Drug,SmallMolecule,Protein&n0.expand=exact_match,close_match,same_as&n0.expand.hops=1&n0.expand.agent_type.not_in=text_mining_agent&n1=HGNC:11923&n2=MONDO:0004979&n2.expand=subclass_of,superclass_of,exact_match,close_match,same_as&n2.expand.hops=2&n2.expand.agent_type.not_in=text_mining_agent&e0=n0-n1&e0.agent_type.not_in=text_mining_agent&e1=n0-n2&e1.predicate=treats_or_applied_or_studied_to_treat&e1.agent_type.not_in=text_mining_agent&paths_per_pair=2&ranking=established&traversal=deep_dive&layout=n0.757.91,n1.557.91,n2.1057.91`
- Composed with:
  `uv run --group viz python scripts/visualization/plot_two_panels.py scripts/visualization/fig10_panel_a.png scripts/visualization/fig10_panel_b.png scripts/visualization/fig10_panel_c.png -o fig10.pdf --size 7 --dpi 526`

## Use case text: full 260-input subgraph (not a figure)

Rerun on v2.3.0 for the numbers quoted in the use case paragraph. POST to `https://kestrel.krakenkg.com/api/subgraph`
with the request saved in `one_off_kraken_usecase/combi_subgraph_v2.3.0_request.json` (result in
`combi_subgraph_v2.3.0.json`):

- 260 input ids (from `combi_subgraph_input_ids.py`), `max_path_length` 2, `limit` 400, default ranking.
- Connectors excluded: `biolink:DiseaseOrPhenotypicFeature`, `biolink:AnatomicalEntity`,
  `biolink:InformationContentEntity`, `biolink:OrganismTaxon`, `biolink:Human`.

Result: 660 nodes, 22,917 edges; 245 genes and proteins, 201 drug-annotated chemicals, 171 chemicals without drug
annotation, 35 processes, 5 clinical findings, 3 cell types. Fourteen connectors are analytes measured in the
source study but discarded by the model.

## Supplementary Data (PDF and xlsx)

`uv run --group viz python scripts/visualization/make_supplementary_tables.py /Volumes/AmySSD/kraken-data/artifacts/metagraphs/kraken_metagraph_2.3.0.json -o supplementary_data.xlsx --example-node scripts/visualization/supplementary_example_node.json --example-edge scripts/visualization/supplementary_example_edge.json --figure-queries scripts/visualization/supplementary_figure_queries.json`

- S1 to S4 come from the metagraph. S5 uses the Figure 1 creatinine node (`RM:0028430`) and edge, copied from the
  v2.3.0 NDJSON files into the two `supplementary_example_*.json` files. S6 is rendered from
  `supplementary_figure_queries.json`, which holds the public-site versions of the links above; update it when a
  figure's query changes.
