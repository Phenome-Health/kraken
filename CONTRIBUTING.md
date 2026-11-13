# Contributing to KRAKEN

## Setup
Requires [uv](https://docs.astral.sh/uv/) and Python ≥3.10. The repo pins 3.12
(`.python-version`) and CI runs on it, so 3.12 is the safe default.
```bash
uv sync --all-groups   # runtime + dev + viz deps
```

## Building a release

### Getting the input data

The input files for the current release are archived in a Phenome Health S3 bucket,
`kraken-kg-archive` (us-west-2), which needs team AWS access. With the AWS CLI, `zstd`
(`brew install zstd`) and an AWS CLI profile named `kraken-admin` (log in with
`aws sso login --profile kraken-admin`):

```bash
scripts/restore_inputs.sh /path/to/kraken-data
```

This downloads `s3://kraken-kg-archive/input_data` (~32 GB, compressed) into
`/path/to/kraken-data/input_data/` and unpacks it into the files the build reads (~250 GB). It
can be re-run: anything already downloaded or unpacked is skipped, so if your AWS login expires
partway, log in again and re-run it.

Then set `base_path` in `config/build_config.yaml` to `/path/to/kraken-data` and run
`uv run python -m kraken.build`. Outputs go to `<base_path>/artifacts/`.

The same bucket also holds the 2.3.0 output (`artifacts/2.3.0/`, zstd-compressed node and edge
files plus `build_info.json`) and the Kestrel data tarballs built from it (`kestrel-tarballs/`).

### Updating sources

Sources are downloaded by hand into `input_data/<source>/` (usually in a version subfolder), with
the version recorded in `config/build_config.yaml` (`sources.<source>.version`). A finished build lists the
versions it used in its `build_info.json`. Where each one came from:

| Source | Where to get it                                                                                                                                                                                                                                                                                                                   |
| --- |-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `babel` | [stars.renci.org/var/babel_outputs](https://stars.renci.org/var/babel_outputs/) (`input_data/babel/<release>/download.sh` fetches a release)                                                                                                                                                                                      |
| `kg2` | RTX-KG2 KG2c: [KG-Registry entry](https://kghub.org/kg-registry/resource/rtx-kg2/rtx-kg2.html)                                                                                                                                                                                                                                    |
| `robokop` | [KG-Registry entry](https://kghub.org/kg-registry/resource/robokop/robokop.html)                                                                                                                                                                                                                                                  |
| `translator-kg-open` | [kgx-storage.ci.transltr.io/releases/translator_kg_open](https://kgx-storage.ci.transltr.io/releases/translator_kg_open/)                                                                                                                                                                                                         |
| `multiomics-kg` | Multiomics KG: download URLs are in the [multiomics_kp manifest](https://github.com/multiomicsKP/multiomics_kp/blob/main/manifest.json)                                                                                                                                                                                           |
| `microbiome-kg` | Microbiome KG: download URLs are in the [microbiome_kp manifest](https://github.com/multiomicsKP/microbiome_kp/blob/main/manifest.json)                                                                                                                                                                                           |
| `ctkg` | Clinical Trials KG: download URLs are in the [clinical_trials_kp manifest](https://github.com/multiomicsKP/clinical_trials_kp/blob/main/manifest.json)                                                                                                                                                                            |
| `dakg` | Drug Approvals KG: download URLs are in the [drug_approvals_kp manifest](https://github.com/multiomicsKP/drug_approvals_kp/blob/main/manifest.json)                                                                                                                                                                               |
| `ncbigene` | [ftp.ncbi.nlm.nih.gov/gene/DATA](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) (`gene_info.gz`, `gene_summary.gz`), plus `taxdump.tar.gz` from [ftp.ncbi.nlm.nih.gov/pub/taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)                                                                                                         |
| `loinc` | [loinc.org/downloads](https://loinc.org/downloads/) (free account)                                                                                                                                                                                                                                                                |
| `umls` | `LOINC_PART_UMLS.tsv`, a LOINC-part-to-UMLS mapping that Gwenlyn made from a [UMLS](https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html) download (from `MRCONSO` file). UMLS needs a personal UTS licence, which can take a few days to approve. Planned to be replaced by a full direct UMLS ingest |
| `lipidmaps` | [LIPID MAPS LMSD](https://www.lipidmaps.org/databases/lmsd/overview) (`LMSD.sdf.zip`)                                                                                                                                                                                                                                             |
| `refmet` | [Metabolomics Workbench RefMet](https://www.metabolomicsworkbench.org/databases/refmet/index.php) (`refmet.csv`)                                                                                                                                                                                                                  |
| `cdes` | [NIH CDE Repository](https://cde.nlm.nih.gov/) search export (`SearchExport.json`)                                                                                                                                                                                                                                                |
| `pgs-catalog` | [PGS Catalog](https://www.pgscatalog.org/) metadata (`pgs_all_metadata.xlsx`). The build downloads each score's harmonized file from EBI's FTP and the Ensembl 112 gene model itself                                                                                                                                              |
| `long-covid` | Long COVID tables from Qi Wei at the Institute for Systems Biology (RECOVER-funded), sent by email. Qi Wei is contact for any new versions.                                                                                                                                                                                       |
| `icd10cm` | CDC [ICD-10-CM order file](https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/ICD10CM/2027/icd10cm-code-descriptions-2027.zip) (`icd10cm-order-2027.txt`), Mondo [SSSOM mappings](https://purl.obolibrary.org/obo/mondo/mappings/mondo.sssom.tsv) saved with the download date, and the [Disease Ontology](https://github.com/DiseaseOntology/HumanDiseaseOntology/releases) release (`doid.obo`), saved with its release date|
| `icd10` | WHO [ICD-10](https://icdcdn.who.int/icd10/index.html), 2019 edition, plain text tabular format (`icd102019enMeta.zip`; the build reads its `icd102019syst_codes.txt`). WHO holds the copyright; its download page shows no licence text. Harmonized after `babel`: it reads the harmonized output of `babel` and of `icd10cm`, so harmonize it again when either changes (the build says so) |
| `bio-bmi` | Watanabe et al., *Nat Med* 2023 ([doi:10.1038/s41591-023-02248-0](https://doi.org/10.1038/s41591-023-02248-0)). Fixed paper; Lance or Noa know the work                                                                                                                                                                           |
| `bio-age` | Earls et al., *J Gerontol A* 2019 ([doi:10.1093/gerona/glz220](https://doi.org/10.1093/gerona/glz220)). Fixed paper; Lance or Noa know the work                                                                                                                                                                                   |

The two bio-bmi and bio-age papers' tables were taken from the publications: one from a supplementary table, the
other typed in by hand from an image of a table (might be possible to get spreadsheet form from authors).

### Hardware

The build trades memory for disk: it streams through large intermediate files instead of holding
the graph in memory, so it runs on a laptop, but **disk speed matters a lot**. Use a fast NVMe
SSD, not a spinning disk or a slow USB drive.

Releases up to 2.3.0 were built on:

- **Laptop:** MacBook Pro, Apple M4 Pro (14 cores), 48 GB memory
- **Data drive:** external Samsung 990 PRO 2 TB NVMe SSD in an OWC Express 1M2 enclosure
  (Thunderbolt / USB4, 40 Gb/s), set as `base_path`

Space needed for one release is roughly 400 GB: ~206 GB of unpacked inputs, ~93 GB of harmonized
intermediates (`artifacts/harmonized/`) and ~82 GB of final node and edge files. A 1 TB drive is
the practical minimum; 2 TB leaves room to keep the previous release alongside. On that setup,
harmonizing all sources takes about an hour and integrating about 2.3 hours.

## Before pushing
```bash
scripts/fix.sh     # auto-format (black) + auto-fix (ruff)
scripts/check.sh   # ruff + black --check + pytest — same as CI
```
Green `check.sh` ⇒ green CI.

## Adding a source
Each source is a `BaseHarmonizer` subclass in `src/kraken/harmonizers/`. Set its
property mappings (or override `harmonize`) and build nodes/edges via
`create_node`/`create_edge` — CURIEs are formed through biomapper2. Then register
it in `orchestrator.py` and add a source entry to `config/build_config.yaml`.

## Commits
We use [Conventional Commits](https://www.conventionalcommits.org/). The type
drives releases (release-please defaults), so it matters:
`feat:` → minor · `fix:`/`docs:`/`perf:` → patch · `feat!:`/`BREAKING CHANGE:` → major ·
`chore:`/`refactor:`/`test:`/`ci:` → no release.

## Releases
Automated via [release-please](https://github.com/googleapis/release-please):
merging to `main` opens/updates a release PR (version bump + changelog);
merging *that* PR cuts the release. Never bump the version by hand.

**Two version numbers.** release-please versions the *code* (`pyproject.toml`). The *graph's*
version is separate: `kraken_version` in `config/build_config.yaml`, set by hand before a build.
They are kept in sync by convention (both 2.3.0 so far), but nothing enforces it, and they can
briefly differ: `build_info.json` records both, and the 2.3.0 graph was built with package
version 2.2.0 because the code release was cut after the build. When releasing a new graph,
set `kraken_version` to the version release-please will give the code, and cut the code release
along with it.

Each GitHub release is archived on Zenodo automatically (concept DOI
[10.5281/zenodo.21940866](https://doi.org/10.5281/zenodo.21940866)), through the Zenodo GitHub
integration on Amy Glen's Zenodo account. Record managers: Amy Glen and Drew Witherington.

## The biomapper2 dependency
Pinned to a [biomapper2](https://github.com/Phenome-Health/biomapper2) release tag in `pyproject.toml`. `uv.lock` holds
the exact commit, and a build that runs harmonization writes it to `build_info.json` (`biomapper2_commit`). To take a newer release, change
the tag, then:
```bash
uv lock --upgrade-package biomapper2 && uv sync
```
