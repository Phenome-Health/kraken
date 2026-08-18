# kraken-build

[![CI](https://github.com/Phenome-Health/kraken/actions/workflows/ci.yml/badge.svg)](https://github.com/Phenome-Health/kraken/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Build system for Phenome Health's unified general-purpose knowledge graph, **KRAKEN** _(Knowledge Research & Analysis Kit for Evidence Networks)_.

## Overview

KRAKEN integrates biomedical knowledge from many sources into a single [Biolink Model](https://biolink.github.io/biolink-model/)-compliant knowledge graph, tailored for multiomic and wellness research. Release 2.1.0 spans roughly **15.3M nodes** and **112.7M edges** drawn from **114 primary knowledge sources**. The build system harmonizes each source to a common schema, resolves entities across sources using provided equivalency mappings, and exports flat files ready for downstream tools.

## Knowledge sources

KRAKEN integrates:

- **Aggregator knowledge graphs** — RTX-KG2, ROBOKOP, and Translator KG Open
- **Ontologies & vocabularies** — UMLS, LOINC, LIPID MAPS, RefMet
- **Specialized resources** — the NIH CDE Repository, PGS Catalog, and published biological-age and biological-BMI datasets

Each aggregator contributes many upstream primary knowledge sources; see `artifacts/metagraphs/` for the full per-source breakdown.

## Usage

Configure sources and build steps in `config/build_config.yaml`, then run:

```bash
uv run python build.py
```

## Build process

1. **Harmonize** — convert each source to a common Biolink format/schema
2. **Integrate** — merge sources, performing entity resolution using each source's equivalency mappings
3. **Post-process** — generate a coherent, small test version of the graph

Validation and metagraph creation are done both on each harmonized artifact and on the final integrated graph.

## Output

- **Unified KG**: `artifacts/integrated/kraken_*.jsonl`
- **Metagraphs**: `artifacts/metagraphs/`

## Requirements

Python 3.10+ with dependencies managed by [`uv`](https://docs.astral.sh/uv/):

```bash
uv sync
```

## Citation

KRAKEN is described in a forthcoming manuscript; please cite that publication. Each release is also archived on Zenodo — concept DOI [10.5281/zenodo.21940866](https://doi.org/10.5281/zenodo.21940866), which always resolves to the latest version.

## License

Released under the [MIT License](LICENSE).

## Acknowledgments

Development of this codebase was assisted by [Claude Code](https://claude.com/claude-code) (Anthropic).
