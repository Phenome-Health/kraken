# Entity Resolution — the process (north star)

This is the intended ER process, in order. Everything in the implementation should
serve this; if the code and this doc disagree, this doc is the north star.

1. Take every id from the equivalent-ids list on the harmonized source files and make a
   separate node in the match graph. Babel is one of those sources — it's ingested with one
   node per id, and its cliques and gene/protein conflations come in as `same_as` edges,
   which count as equivalence evidence just like an equiv-ids list. Where an aggregator
   (kg2 / ROBOKOP / Translator) asserts two ids are equivalent — in an equiv-ids list or a
   `same_as` / `exact_match` edge — and Babel knows both ids, Babel decides; the aggregator's
   claim only counts for ids Babel doesn't know — except for a Babel clique OUTLIER: an id whose
   name matches none of its clique's (which agree among themselves) but two or more of another
   clique's. For those, Babel's evidence is dropped and the aggregators' claims count again
   (see `babel_outliers`). Ids whose names match also get linked:
   ignoring spacing, hyphens and possessives, and plurals outside chemistry, but never
   stereo signs, charges or primes (see `name_sim.name_keys`). A name counts only when it is
   attributable to the id itself: a native source names its own id, while a canonicalized
   aggregator's node carries the clique's preferred label, so its name is used only for an id
   Babel has never heard of — where no Babel node will ever name it, and without it the id
   carries no name into the match graph at all and can never match anything.

2. Set those nodes' id / category / name / taxon based on Babel's node for that exact id
   (its own label, etc.), backing up to deriving it from the source if possible.

3. Convert each of those individual nodes' categories to a "family" per our curated
   mappings.

4. Remove all conflicting edges from the match graph (per ALL of our guardrails).

5. Now we have a clean match graph that's ready to go. Cluster it using label propagation
   to assign a canonical (cluster) id to every individual node in the match graph.

6. Split any resulting clusters that still violate a guardrail because of transitive
   relations (conflicts the pairwise edge-prune in step 4 couldn't see). The split regrows the
   cluster along its strongest edges, starting from whole Babel cliques rather than from single
   ids — so a cluster holding two structures is cut BETWEEN them, instead of shedding whichever
   id is most weakly attached (which is usually the structure id itself, since it often carries
   nothing but clique edges). A clique that breaks a guardrail on its own is still cut: Babel
   cliques do mix two species, and keeping one whole must never override a guardrail.

   Then put back any Babel gene/protein clique that step 5 tore apart: label propagation
   finds communities, and a gene's clique conflated with its protein's is two dense halves
   at one uniform weight, so it splits whenever only one half has outside evidence (ACE).
   A gene and the protein it encodes in one species are one entity, so those clusters are
   rejoined — as long as the guardrails allow the union, which is what keeps the cliques
   that mix two species' orthologs apart. Only this family: elsewhere a torn Babel clique
   is usually label propagation correctly separating a class from its members.

7. Now we're ready to merge nodes. The "canonical" id for the merged node is chosen based
   on our prefix ranking. The chosen name should preferably be the label for that exact
   chosen id per Babel, backing up as appropriate. The description preferably comes
   from the same id as the name, backing up as appropriate. Synonyms should be retained
   from sources IF their equiv ids are a subset of the merged node's equiv ids. Retain the
   taxon from Babel as well — basically never throw away anything from Babel or from
   sources.

8. Now nodes are properly merged. Go through and remap edges appropriately: for aggregator
   sources, look at their ORIGINAL subj/obj, map that to what is now canonical, and
   duplicate edges if necessary (like for kg2, which can have multiple original subj/obj
   per edge).

9. Retain the equivalence signal as edges. For every asserted equivalence (an equiv-ids
   list, a Babel clique, or any source's `same_as` / `exact_match` edge) whose two ids
   ended up in DIFFERENT clusters, keep a `biolink:close_match` edge between their
   representatives — so e.g. TP53 protein-isoforms that didn't merge into the main TP53
   node stay linked to it. (Not `same_as` — we decided they aren't the same thing.)
   Equiv-list assertions carry the source as the primary knowledge source; Babel's cliques
   carry Babel (`infores:sri-node-normalizer`) as the primary knowledge source. Do
   NOT emit name-similarity edges. Globally drop self-edges (an edge whose endpoints merged
   into the same node).
