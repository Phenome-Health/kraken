"""
Configuration models for KRAKEN build system
"""

from pathlib import Path
from typing import Self

from pydantic import BaseModel, Field, model_validator

from kraken.utils.constants import PROJECT_ROOT
from kraken.utils.general import to_list


class HarmonizationConfig(BaseModel):
    output_directory: str
    zip_inputs_after: bool = False


class IntegrationConfig(BaseModel):
    output_directory: str


class MetagraphConfig(BaseModel):
    output_directory: str


class EntityResolutionConfig(BaseModel):
    """Clustering-based entity resolution (opt-in; replaces integrate.py's node
    merge). Default OFF — the legacy merge remains the build default until this is
    validated. See src/kraken/entity_resolution/ and docs/entity_resolution_plan.md.

    Curated inputs (branch families, prefix ranking, weights, ground truth) live at
    fixed paths under config/entity_resolution/ and are loaded directly by the ER
    modules — they are not build-config knobs.
    """

    enabled: bool = False


class StepsConfig(BaseModel):
    harmonize: bool
    integrate: bool
    postprocess: bool


class OptionsConfig(BaseModel):
    include_sources: str | list[str] | None = None
    exclude_sources: str | list[str] | None = None
    metagraph_creation: bool = False
    validate_output: bool = True
    validation_only: bool = False
    # Ask for interactive confirmation of the source versions before building. Automatically skipped when
    # stdin is not a TTY (e.g. cron/CI), so it never blocks non-interactive runs.
    confirm_source_versions: bool = True


class TestExportConfig(BaseModel):
    output_directory: str
    num_edges: int


class PostProcessingConfig(BaseModel):
    test_export: TestExportConfig | None = None


class SourceConfig(BaseModel):
    version: str | None = None  # version/release of the source that was ingested (e.g. "2.10.2", "june2025")
    input_file: str | None = None
    nodes_input: str | None = None
    edges_input: str | None = None

    # Computed (set by KrakenConfig validator)
    input_file_resolved: Path | None = Field(default=None, init=False)
    nodes_input_resolved: Path | None = Field(default=None, init=False)
    edges_input_resolved: Path | None = Field(default=None, init=False)

    def resolve(self, base_path: Path) -> None:
        """Resolve paths against base_path (mutates in place)"""
        if self.input_file:
            self.input_file_resolved = base_path / self.input_file
        if self.nodes_input:
            self.nodes_input_resolved = base_path / self.nodes_input
        if self.edges_input:
            self.edges_input_resolved = base_path / self.edges_input


class KrakenConfig(BaseModel):
    biolink_version: str
    kraken_version: str
    kg_label: str | None = None  # human-readable build name, e.g. "kraken-lite"
    log_level: str = "INFO"
    base_path: str | None = None
    harmonization: HarmonizationConfig
    integration: IntegrationConfig
    metagraph: MetagraphConfig
    steps: StepsConfig
    options: OptionsConfig
    post_processing: PostProcessingConfig | None = None
    entity_resolution: EntityResolutionConfig = Field(default_factory=EntityResolutionConfig)
    sources: dict[str, SourceConfig]

    # Computed field (not from yaml)
    sources_to_use: set[str] = Field(default_factory=set, init=False)

    @model_validator(mode="after")
    def validate_and_resolve_sources(self) -> Self:
        all_sources = set(self.sources.keys())

        include_sources = set(to_list(self.options.include_sources))
        exclude_sources = set(to_list(self.options.exclude_sources))

        if not include_sources.issubset(all_sources) or not exclude_sources.issubset(all_sources):
            raise ValueError(
                f"Source names in 'include_sources' and 'exclude_sources' must exist under 'sources': {all_sources}"
            )

        overlap = include_sources.intersection(exclude_sources)
        if overlap:
            raise ValueError(f"Sources cannot be both included and excluded: {overlap}")

        source_pool = include_sources if include_sources else all_sources
        self.sources_to_use = source_pool - exclude_sources

        # Every source actually included in the build must declare a version (recorded in build_info.json
        # and the metagraphs for provenance). Only in-use sources are checked, so an excluded or unused
        # source with an unknown version never blocks the build.
        missing_versions = sorted(name for name in self.sources_to_use if not self.sources[name].version)
        if missing_versions:
            raise ValueError(
                f"These sources are included in the build but have no 'version' set under 'sources' in the "
                f"build config: {missing_versions}. Set a version for each (a release, a download date, or an "
                f"explicit value like 'unknown'). Only sources actually included in the build are checked."
            )

        # Resolve paths for each source
        for source_name in self.sources_to_use:
            self.sources[source_name].resolve(self.base_path_resolved)

        return self

    # Computed path properties

    @property
    def base_path_resolved(self) -> Path:
        return Path(self.base_path) if self.base_path else PROJECT_ROOT

    @property
    def harmonized_dir(self) -> Path:
        return self.base_path_resolved / self.harmonization.output_directory

    @property
    def metagraph_dir(self) -> Path:
        return self.base_path_resolved / self.metagraph.output_directory

    @property
    def integrated_dir(self) -> Path:
        return self.base_path_resolved / self.integration.output_directory

    @property
    def test_export_dir(self) -> Path:
        return self.base_path_resolved / self.post_processing.test_export.output_directory

    @property
    def integrated_debug_dir(self) -> Path:
        return self.integrated_dir / "debug"

    # Entity-resolution outputs (curated inputs load from fixed config paths in
    # the ER modules themselves). These resolve against base_path.
    @property
    def er_nodenorm_cache_path(self) -> Path:
        return self.base_path_resolved / "artifacts" / "entity_resolution" / "nodenorm_cache.sqlite"

    @property
    def er_membership_path(self) -> Path:
        return self.integrated_dir / f"kraken_membership_{self.kraken_version}.jsonl"

    @property
    def integrated_nodes_path(self) -> Path:
        return self.integrated_dir / f"kraken_nodes_{self.kraken_version}.jsonl"

    @property
    def integrated_edges_path(self) -> Path:
        return self.integrated_dir / f"kraken_edges_{self.kraken_version}.jsonl"

    @property
    def source_versions(self) -> dict[str, str | None]:
        """Map each source in use to its configured version (None if unknown/unspecified).

        Recorded in build_info.json so a KG release can be traced back to the exact
        source versions that went into it.
        """
        return {source: self.sources[source].version for source in sorted(self.sources_to_use)}

    @property
    def create_metagraphs(self) -> bool:
        return self.options.metagraph_creation

    @property
    def zip_inputs_after(self) -> bool:
        return self.harmonization.zip_inputs_after

    @property
    def all_harmonized_paths_resolved(self) -> dict[str, tuple[Path, Path]]:
        """Get harmonized paths for all sources to use"""
        return {
            source: (
                self.harmonized_dir / source / "nodes.jsonl",
                self.harmonized_dir / source / "edges.jsonl",
            )
            for source in self.sources_to_use
        }

    @property
    def all_source_input_paths_resolved(self) -> dict[str, list[Path]]:
        """Get resolved input paths for all sources to use"""
        return {
            source: [
                p
                for p in [
                    self.sources[source].input_file_resolved,
                    self.sources[source].nodes_input_resolved,
                    self.sources[source].edges_input_resolved,
                ]
                if p
            ]
            for source in self.sources_to_use
        }
