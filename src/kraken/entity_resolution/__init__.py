"""Clustering-based entity resolution for KRAKEN.

Replaces the order-dependent node merge in ``integrate.py`` (primary_source +
can_merge_existing_nodes) with an order-independent pipeline whose output is a
function of the match graph alone:

    match_graph  -> weighted CURIE-pair evidence (equivalency cliques, source
                    match predicates, name similarity), accumulated across
                    sources and pre-filtered at tau
    clustering   -> per connected component: Babel cliques start whole, then
                    the strongest links merge first and equal ones are decided
                    together, on exact integer weights (no seed, no id order)
    guardrails   -> hereditary "one branch / one taxon / one structural id"
                    rules that no merge may break, checked again on the output
    materialize  -> deterministic representative selection + order-independent
                    property reconciliation

The single artifact the rest of the build consumes is a membership map
``curie -> cluster_id`` (see ``resolve``), which drops into the existing
order-independent edge integration (``integrate.integrate_edges``) in place of
the old ``equivalency_index``.

See ``docs/entity_resolution_plan.md`` for the requirements and rationale.
"""
