# clingen.py
import json
import logging
import re
import sys
import time
from pathlib import Path
from typing import Any

import requests
from biomapper2.core.normalizer import Normalizer
from bs4 import BeautifulSoup

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.biolink_client import BiolinkClient
from kraken.utils.constants import CLINGEN_INFORES, ID
from kraken.utils.general import create_edge_key
from kraken.utils.kg_io import save_to_jsonl

CLINGEN_ACMG_API_URL = "https://actionability.clinicalgenome.org/ac/api/summ/brief"


class ClinGenHarmonizer(BaseHarmonizer):
    """
    Harmonizer for ClinicalGenome ACMG actionability data.

    Processes gene-disease and variant-disease associations from the
    ACMG Clinical Genome Resource (ClinGen) Actionability Working Group.
    """

    source_name = "clingen"
    source_infores = CLINGEN_INFORES

    def __init__(self, biolink_client: BiolinkClient):
        super().__init__(biolink_client)
        self.normalizer = Normalizer(biolink_version=biolink_client.version)

    def harmonize(
        self,
        nodes_output: Path,
        edges_output: Path,
        *,
        input_file: Path | None = None,
        nodes_input: Path | None = None,
        edges_input: Path | None = None,
    ):
        """
        Harmonize ClinGen ACMG actionability data.
        """
        if not input_file:
            raise ValueError(f"{self.source_name} requires input_file")

        logging.info(f"Harmonizing {self.source_name}: {input_file} -> {nodes_output}, {edges_output}")

        nodes: dict[str, dict[str, Any]] = {}
        edges: dict[str, dict[str, Any]] = {}

        # Load the data (either from pre-downloaded file or fetch from API)
        data = self._load_or_fetch_data(input_file)

        # Process each record
        for record in data:
            self._process_record(record, nodes, edges)

        logging.info(f"Saving {len(nodes)} ClinGen ACMG nodes and {len(edges)} edges")
        save_to_jsonl(nodes.values(), nodes_output, mode="w")
        save_to_jsonl(edges.values(), edges_output, mode="w")

        logging.info(f"{self.source_name} harmonization complete: {len(nodes)} nodes, {len(edges)} edges")

    def _load_or_fetch_data(self, input_file: Path) -> list[dict[str, Any]]:
        """Load data from file or fetch from API if file doesn't exist."""
        if input_file.exists():
            logging.info(f"Loading ClinGen ACMG data from {input_file}")
            with open(input_file) as f:
                return json.load(f)
        else:
            logging.info(f"Fetching ClinGen ACMG data from {CLINGEN_ACMG_API_URL}")
            response = requests.get(CLINGEN_ACMG_API_URL)
            response.raise_for_status()
            data = response.json()

            # Save the fetched data for future use
            input_file.parent.mkdir(parents=True, exist_ok=True)
            with open(input_file, "w") as f:
                json.dump(data, f, indent=2)
            logging.info(f"Saved ClinGen ACMG data to {input_file}")

            return data

    def _process_record(
        self, record: dict[str, Any], nodes: dict[str, dict[str, Any]], edges: dict[str, dict[str, Any]]
    ) -> None:
        """Process a single ClinGen record."""
        doc_id = record.get("docId")
        modes_of_inheritance = record.get("modesOfInheritance", [])
        context = record["context"]
        context_type = "Adult"  # SKIP pediatric for now...
        context_data = context.get(context_type)

        # Skip any records with incomplete status (they generally don't have html pages..)
        if not context_data or context_data["status"].get("stg2", "").lower() == "incomplete":
            return

        logging.info(f"On record {doc_id}")

        # Scrape actionability scores once per record
        source_iri = f"https://actionability.clinicalgenome.org/ac/{context_type}/ui/stg2SummaryRpt?doc={doc_id}"
        assertions_by_gene, outcome_intervention_pairs = self._scrape_actionability_scores(source_iri, doc_id)

        # Process gene-disease associations
        genes = context_data.get("genes", [])
        for gene_info in genes:
            gene_symbol = gene_info["gene"]
            gene_omim = gene_info["geneOmim"]

            # Create gene node
            gene_node = self._process_gene(gene_symbol, gene_omim)
            nodes[gene_node[ID]] = gene_node

            # Process diseases associated with this gene
            diseases = gene_info.get("diseases", [])
            for disease_info in diseases:
                # Create disease node from omim and preferredMondo identifiers
                disease_node = self._process_disease(disease_info, gene_omim)
                if disease_node:
                    nodes[disease_node[ID]] = disease_node

                    # Create gene-disease edge with assertions specific to this gene-condition pair
                    edge = self._create_gene_disease_edge(
                        gene_node[ID],
                        disease_node[ID],
                        gene_symbol,
                        disease_info.get("omim"),
                        doc_id,
                        context_type,
                        modes_of_inheritance,
                        outcome_intervention_pairs,
                        assertions_by_gene,
                    )
                    if edge:
                        edge_key = create_edge_key(edge)
                        edges[edge_key] = edge

    def _process_disease(self, disease_info: dict[str, Any], gene_omim: str) -> dict[str, Any] | None:
        """
        Process a disease/condition into a harmonized node.
        Uses the omim and preferredMondo identifiers from the disease_info.
        """
        disease_label = disease_info.get("label", "")
        disease_omim = disease_info.get("omim")
        disease_mondo = disease_info.get("preferredMondo")

        # Build a dict of identifiers to normalize
        id_dict = {}
        if disease_mondo:
            id_dict["mondo"] = disease_mondo
        if disease_omim and disease_omim != gene_omim:  # Sometimes they incorrectly give the gene OMIM on the disease
            id_dict["omim"] = disease_omim

        if not id_dict:
            logging.error(f"No disease identifiers found for: {disease_label}")
            sys.exit(1)

        # Normalize disease identifiers to standard curies
        disease_curies_dict, _ = self.normalizer.get_curies(id_dict, stop_on_invalid_id=True)

        if disease_curies_dict:
            disease_curie = list(sorted(disease_curies_dict.keys(), reverse=True))[
                0
            ]  # OMIM identifiers seem more accurate
            disease_iri = disease_curies_dict[disease_curie]
            equivalent_ids = list(disease_curies_dict.keys())
        else:
            # Fallback: if normalization fails, use the MONDO or OMIM directly
            logging.error(
                f"Could not normalize disease: {disease_label} with IDs {id_dict}. full disease item is: {disease_info}"
            )
            sys.exit(1)

        return self.create_node(
            curie=disease_curie,
            categories=["biolink:Disease"],
            equivalent_ids=equivalent_ids,
            provided_by=self.source_infores,
            name=disease_label,
            urls=disease_iri,
            synonyms=[disease_label] if disease_label else None,
        )

    def _process_gene(self, gene_symbol: str, gene_omim: str) -> dict[str, Any]:
        """Process a gene into a harmonized node."""
        # Build a dict of identifiers to normalize
        id_dict = {"omim": gene_omim}

        # Normalize to standard gene identifiers (HGNC, NCBIGene, etc.)
        gene_curies_dict, _ = self.normalizer.get_curies(id_dict, stop_on_invalid_id=True)

        if gene_curies_dict:
            gene_curie = list(gene_curies_dict.keys())[0]
            gene_iri = gene_curies_dict[gene_curie]
            equivalent_ids = list(gene_curies_dict.keys())
        else:
            # Fallback: use HGNC symbol as identifier
            logging.error(f"Could not normalize gene: {gene_symbol} - gene omim: {gene_omim}")
            sys.exit(1)

        return self.create_node(
            curie=gene_curie,
            categories=["biolink:Gene"],
            equivalent_ids=equivalent_ids,
            provided_by=self.source_infores,
            name=gene_symbol,
            urls=gene_iri,
            synonyms=[gene_symbol],
        )

    def _create_gene_disease_edge(
        self,
        gene_id: str,
        disease_id: str,
        gene_symbol: str,
        disease_local_ids: list[str],
        doc_id: str,
        context_type: str,
        modes_of_inheritance: list[str],
        outcome_intervention_pairs: list[dict[str, Any]],
        assertions_by_gene: dict[str, list[dict[str, str]]],
    ) -> dict[str, Any] | None:
        """
        Create an edge representing a gene-disease association.

        Args:
            gene_id: Curie for the gene
            disease_id: Curie for the disease
            gene_symbol: Gene symbol (for matching assertions)
            disease_local_ids: OMIM and/or MONDO local IDs for the disease (for matching assertions)
            doc_id: ClinGen document ID
            context_type: Context (Adult/Pediatric)
            modes_of_inheritance: List of inheritance modes
            outcome_intervention_pairs: List of outcome-intervention pair dictionaries
            assertions_by_gene: Dictionary mapping gene symbols to their assertions for conditions
        """
        source_iri = f"https://actionability.clinicalgenome.org/ac/{context_type}/ui/stg2SummaryRpt?doc={doc_id}"

        # Find the assertion specific to this gene-condition pair
        matching_assertion = self._find_matching_assertion(gene_symbol, disease_local_ids, assertions_by_gene)

        # Build attributes
        attributes = {
            "doc_id": doc_id,
            "modes_of_inheritance": modes_of_inheritance,
            "source_iri": source_iri,
            "source_iri_json": f"https://actionability.clinicalgenome.org/ac/{context_type}/api/sepio/doc/{doc_id}",
        }

        # Add gene-condition specific assertion if found
        if matching_assertion:
            attributes["gene_condition_assertion"] = matching_assertion

        # Add outcome-intervention pairs (same for all edges from this record)
        if outcome_intervention_pairs:
            attributes["outcome_intervention_pairs"] = outcome_intervention_pairs

        return self.create_edge(
            subject_id=gene_id,
            object_id=disease_id,
            predicate="biolink:contributes_to",
            primary_ks=self.source_infores,
            knowledge_level="knowledge_assertion",
            agent_type="manual_agent",
            qualifiers={"context_qualifier": context_type.lower()},
            attributes=attributes,
        )

    def _scrape_actionability_scores(
        self, url: str, doc_id: str, debug: bool = False
    ) -> tuple[dict[str, list[dict[str, str]]], list[dict[str, Any]]]:
        """
        Scrape actionability scores and final assertions from a ClinGen HTML page.

        Args:
            url: Full URL to the ClinGen actionability page
                 e.g., "https://actionability.clinicalgenome.org/ac/Adult/ui/stg2SummaryRpt?doc=AC102"
            doc_id: Document ID in ClinGen for the given record (used in logging)
            debug: If True, print the HTML for inspection

        Returns:
            Tuple of (assertions_by_gene, outcome_intervention_pairs) where:
                - assertions_by_gene: Dict mapping gene symbols to lists of assertion dicts
                - outcome_intervention_pairs: List of outcome-intervention pair dicts
            Returns ({}, []) if scraping fails
        """
        logging.info(f"Scraping actionability assertions/scores from html at {url}")
        start_time = time.time()

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            soup = BeautifulSoup(response.content, "html.parser")

            final_assertions: list[dict[str, str]] = []
            outcome_intervention_pairs: list[dict[str, Any]] = []

            # Scrape Final Assertions
            page_text = soup.get_text()

            # Look for the "Actionability Assertion" section
            if "Actionability Assertion" in page_text:
                # Find the section in the soup
                assertion_marker = soup.find(string=lambda text: text and "Gene Condition Pairs(s)" in text)

                if assertion_marker:
                    # Navigate to find the parent container
                    container = assertion_marker.find_parent()
                    if container:
                        current = container.find_next_sibling()

                        gene_condition_pattern: list[str] = []
                        assertion_pattern: list[str] = []

                        # Walk through siblings collecting data
                        while current and "Final Consensus Scores" not in current.get_text():
                            text = current.get_text(strip=True)

                            if debug:
                                print(f"Examining: {text[:100]}")

                            # Check if this looks like a gene-condition pair (contains ⇔)
                            if "⇔" in text:
                                gene_condition_pattern.append(text)
                                if debug:
                                    print(f"  -> Gene-condition: {text}")
                            # Check if this looks like an assertion (contains "Actionability")
                            elif "Actionability" in text and len(text) < 50:
                                assertion_pattern.append(text)
                                if debug:
                                    print(f"  -> Assertion: {text}")

                            current = current.find_next_sibling()

                        # Pair them up (they should alternate)
                        for i in range(min(len(gene_condition_pattern), len(assertion_pattern))):
                            final_assertions.append(
                                {"gene_condition_pair": gene_condition_pattern[i], "assertion": assertion_pattern[i]}
                            )

                # Alternative approach: use regex to find patterns in the full text
                if not final_assertions:
                    # Extract the assertion section
                    assertion_section_match = re.search(
                        r"Actionability Assertion.*?(?=Actionability Rationale|Final Consensus Scores)",
                        page_text,
                        re.DOTALL,
                    )

                    if assertion_section_match:
                        assertion_text = assertion_section_match.group(0)

                        if debug:
                            print("Found assertion section:")
                            print(assertion_text[:500])

                        # Find all gene-condition pairs (contain ⇔ and a number)
                        gene_cond_matches = re.findall(r"([A-Z0-9]+⇔\d+[^\n]*)", assertion_text)

                        # Find all actionability levels
                        assertion_matches = re.findall(r"(\w+(?:\s+\w+)?\s+Actionability)", assertion_text)

                        if debug:
                            print(f"\nFound {len(gene_cond_matches)} gene-condition pairs")
                            print(f"Found {len(assertion_matches)} assertions")
                            print(f"Gene-condition pairs: {gene_cond_matches}")
                            print(f"Assertions: {assertion_matches}")

                        # Pair them up
                        for i in range(min(len(gene_cond_matches), len(assertion_matches))):
                            final_assertions.append(
                                {
                                    "gene_condition_pair": gene_cond_matches[i].strip(),
                                    "assertion": assertion_matches[i].strip(),
                                }
                            )

            if not final_assertions and debug:
                print("WARNING: No final assertions found")

            # Organize assertions by gene symbol
            assertions_by_gene: dict[str, list[dict[str, str]]] = {}
            for assertion in final_assertions:
                gene_condition_pair = assertion.get("gene_condition_pair", "")

                # Parse the gene symbol from the pair (e.g., "BRCA1⇔114480 (...)")
                if "⇔" in gene_condition_pair:
                    gene_symbol = gene_condition_pair.split("⇔")[0].strip()

                    if gene_symbol not in assertions_by_gene:
                        assertions_by_gene[gene_symbol] = []

                    assertions_by_gene[gene_symbol].append(
                        {"gene_condition_pair": gene_condition_pair, "assertion": assertion.get("assertion", "")}
                    )

            # Now find the "Final Consensus Scores" section for outcome/intervention pairs
            scr_table = soup.find("div", class_="scrTable")

            if debug:
                print("=" * 80)
                print("SCORING TABLE:")
                print("=" * 80)
                if scr_table:
                    print(scr_table.prettify())
                else:
                    print("No scrTable found!")
                print("=" * 80)

            if not scr_table:
                elapsed_time = time.time() - start_time
                logging.warning(f"No scoring table found for doc {doc_id} (took {elapsed_time:.2f}s)")
            else:
                # Find all data rows (rows with class "data")
                data_rows = scr_table.find_all("div", class_="data row")

                if debug:
                    print(f"\nFound {len(data_rows)} data rows")

                for row in data_rows:
                    oi_pair = row.find("div", class_=lambda x: x and "oiPair" in x)
                    severity = row.find("div", class_=lambda x: x and "severity" in x and "scrData" in x)
                    likelihood = row.find("div", class_=lambda x: x and "likelihood" in x and "scrData" in x)
                    effectiveness_div = row.find_all("div", class_=lambda x: x and "scrData" in x)
                    # The third scrData div is typically effectiveness
                    effectiveness = effectiveness_div[2] if len(effectiveness_div) > 2 else None
                    noi = row.find("div", class_=lambda x: x and "noi" in x and "scrData" in x)
                    total_score = row.find("div", class_=lambda x: x and "totalScore" in x and "scrData" in x)

                    if debug:
                        print("\nProcessing row:")
                        print(f"  oi_pair: {oi_pair.get_text(strip=True) if oi_pair else 'None'}")
                        print(f"  severity: {severity.get_text(strip=True) if severity else 'None'}")
                        print(f"  likelihood: {likelihood.get_text(strip=True) if likelihood else 'None'}")
                        print(f"  effectiveness: {effectiveness.get_text(strip=True) if effectiveness else 'None'}")
                        print(f"  noi: {noi.get_text(strip=True) if noi else 'None'}")
                        print(f"  total: {total_score.get_text(strip=True) if total_score else 'None'}")

                    if oi_pair:
                        pair_text = oi_pair.get_text(strip=True)

                        # Split by " / " to separate outcome and intervention
                        if " / " in pair_text:
                            parts = pair_text.split(" / ", 1)
                            outcome = parts[0].strip()
                            intervention = parts[1].strip()

                            pair = {
                                "outcome": outcome,
                                "intervention": intervention,
                                "scores": {
                                    "severity": severity.get_text(strip=True) if severity else None,
                                    "likelihood": likelihood.get_text(strip=True) if likelihood else None,
                                    "effectiveness": effectiveness.get_text(strip=True) if effectiveness else None,
                                    "nature_of_intervention": noi.get_text(strip=True) if noi else None,
                                    "total": total_score.get_text(strip=True) if total_score else None,
                                },
                            }
                            outcome_intervention_pairs.append(pair)

            elapsed_time = time.time() - start_time

            # Log warnings if we didn't find expected data
            if not assertions_by_gene:
                logging.warning(f"No final assertions found for doc {doc_id}")
            if not outcome_intervention_pairs:
                logging.warning(f"No outcome/intervention pairs found for doc {doc_id}")

            # Return empty structures if we found neither assertions nor pairs
            if not assertions_by_gene and not outcome_intervention_pairs:
                logging.warning(f"No data found for doc {doc_id} (took {elapsed_time:.2f}s)")
                return {}, []

            logging.info(
                f"Successfully scraped {sum(len(v) for v in assertions_by_gene.values())} gene-condition assertions "
                f"and {len(outcome_intervention_pairs)} outcome-intervention pairs from doc "
                f"{doc_id} (took {elapsed_time:.2f}s)"
            )
            return assertions_by_gene, outcome_intervention_pairs

        except Exception as e:
            elapsed_time = time.time() - start_time
            logging.error(f"Error scraping doc {doc_id} after {elapsed_time:.2f}s: {e}")
            return {}, []

    def _find_matching_assertion(
        self, gene_symbol: str, disease_ids: list[str], assertions_by_gene: dict[str, list[dict[str, str]]]
    ) -> dict[str, str] | None:
        """
        Find the assertion that matches the current gene-condition edge.

        Args:
            gene_symbol: Gene symbol (e.g., "BRCA1")
            disease_ids: OMIM and/or MONDO local IDs for the disease (e.g., ["114480"])
            assertions_by_gene: Dictionary of assertions organized by gene symbol

        Returns:
            The matching assertion dictionary, or None if no match found
        """
        if gene_symbol not in assertions_by_gene:
            return None

        gene_assertions = assertions_by_gene[gene_symbol]

        # Look for the matching disease
        for assertion in gene_assertions:
            gene_condition_pair = assertion.get("gene_condition_pair", "")
            if any(disease_id in gene_condition_pair for disease_id in disease_ids):
                return assertion

        # If we couldn't match by disease ID, log a warning and return nothing
        logging.warning(
            f"Multiple assertions found for gene {gene_symbol} but couldn't match "
            f"to disease ID(s) {disease_ids}. Skipping final assertion for this gene-condition pair."
        )
        return None
