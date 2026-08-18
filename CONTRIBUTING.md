# Contributing to KRAKEN

## Setup
Requires [uv](https://docs.astral.sh/uv/) and Python 3.12 (CI's version).
```bash
uv sync --all-groups   # runtime + dev + viz deps
```

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

## The biomapper2 dependency
Pinned as a git dep (`main`) in `pyproject.toml`. To pick up [biomapper2](https://github.com/Phenome-Health/biomapper2)
changes, merge them there first, then:
```bash
uv lock --upgrade-package biomapper2 && uv sync
```
