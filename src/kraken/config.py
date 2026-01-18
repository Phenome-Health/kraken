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
    primary_source: str


class MetagraphConfig(BaseModel):
    output_directory: str


class StepsConfig(BaseModel):
    harmonize: bool
    integrate: bool
    postprocess: bool


class OptionsConfig(BaseModel):
    include_sources: str | list[str] | None = None
    exclude_sources: str | list[str] | None = None
    metagraph_creation: bool = False


class TestExportConfig(BaseModel):
    output_directory: str
    num_edges: int


class PostProcessingConfig(BaseModel):
    test_export: TestExportConfig | None = None


class SourceConfig(BaseModel):
    input_file: str | None = None
    nodes_input: str | None = None
    edges_input: str | None = None
    can_merge_existing_nodes: bool


class KrakenConfig(BaseModel):
    biolink_version: str
    kraken_version: str
    log_level: str = "INFO"
    base_path: str | None = None
    harmonization: HarmonizationConfig
    integration: IntegrationConfig
    metagraph: MetagraphConfig
    steps: StepsConfig
    options: OptionsConfig
    post_processing: PostProcessingConfig | None = None
    sources: dict[str, SourceConfig]

    # Computed field (not from yaml)
    sources_to_use: set[str] = Field(default_factory=set, init=False)

    @model_validator(mode="after")
    def validate_and_resolve_sources(self) -> Self:
        all_sources = set(self.sources.keys())

        # Validate primary_source exists
        if self.integration.primary_source not in all_sources:
            raise ValueError(
                f"primary_source '{self.integration.primary_source}' must exist under 'sources': {all_sources}"
            )

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
    def integrated_debug_dir(self) -> Path:
        return self.integrated_dir / "debug"

    @property
    def integrated_nodes_path(self) -> Path:
        return self.integrated_dir / f"kraken_nodes_{self.kraken_version}.jsonl"

    @property
    def integrated_edges_path(self) -> Path:
        return self.integrated_dir / f"kraken_edges_{self.kraken_version}.jsonl"

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
    def all_source_input_paths_resolved(self) -> dict[str, dict[str, Path]]:
        """Get full input file paths for all sources to use"""
        result = {}
        for source_name in self.sources_to_use:
            source_config = self.sources[source_name]
            paths = {}
            if source_config.input_file:
                paths["input_file"] = self.base_path_resolved / source_config.input_file
            if source_config.nodes_input:
                paths["nodes_input"] = self.base_path_resolved / source_config.nodes_input
            if source_config.edges_input:
                paths["edges_input"] = self.base_path_resolved / source_config.edges_input
            result[source_name] = paths
        return result
