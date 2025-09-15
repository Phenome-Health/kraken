# kraken-build

Build system for Phenome Health's unified general-purpose knowledge graph, **KRAKEN** _(Knowledge Research & Analysis Kit for Evidence Networks)_.

## Overview

KRAKEN combines knowledge from RTX-KG2, SPOKE, UMLS, LIPID MAPS, and RefMet into a unified biomedical knowledge graph. The system harmonizes data formats, resolves entities across sources, and exports to multiple flat-file formats (ready for import into downstream tools).

## Knowledge Sources

- **KG2**: RTX-KG2 biomedical knowledge graph
- **SPOKE**: Scalable Precision medicine Open Knowledge Engine
- **UMLS**: Unified Medical Language System
- **LIPID MAPS**: Lipid structure database
- **RefMet**: Reference list of metabolite names


## Usage
First configure sources and build steps in `config/build_config.yaml` as desired. Then run:

```bash
uv run python build.py
```

## Build Process

The build process has three main phases:

1. **Harmonize**: Convert each source to common Biolink format/schema
2. **Integrate**: Merge sources using entity resolution, leveraging equivalency mappings provided by each source
3. **Post-process**: Export flat files for ArangoDB and Biomapper

## Output

- **Unified KG**: `artifacts/integrated/kraken_*.jsonl`
- **ArangoDB**: `artifacts/export/arango/`
- **BiomapperKG**: `artifacts/export/biomapper/`
- **Metagraphs**: `artifacts/metagraphs/`

## Requirements

Python 3.10+ with dependencies managed by `uv`.