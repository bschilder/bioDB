"""biodb: biomedical knowledge graph helpers, ported from ``AoU.phenome``.

Module map:

* :mod:`biodb.opentargets` -- Open Targets Platform bulk downloaders
  and parsers (disease/drug/PGx/expression/essentiality/pathways) +
  gene-association matrix builders.
* :mod:`biodb.opentargets_graphql` -- Targeted Open Targets GraphQL
  queries (``query_target``, ``query_disease``, ``query_drug``,
  ``query_variant``).
* :mod:`biodb.monarch` -- Monarch Initiative association readers
  (causal gene-to-disease and friends).
* :mod:`biodb.ontology` -- OBO / OWL / Mondo loaders, N-hop keyword
  set expansion, hierarchical keyword set generation, attention
  analysis, gene-phenotype matrix construction.
* :mod:`biodb.ontology_owl` -- Generic owlready2-based primitives that
  work for any OBO Foundry OWL file: SO, HPO, EFO, GO, ChEBI, … (load,
  descendants/ancestors walk, most-recent common ancestor).
* :mod:`biodb.uniprot` -- UniProt REST client (protein sequences,
  features, cross-references).
* :mod:`biodb.harmonizome` -- Maayan-Lab Harmonizome client (~114
  curated gene-attribute datasets; ``list_datasets``, ``download_datasets``,
  ``get_gmt``, ``load_gene_attribute_matrix``, ``get_dataset_metadata``).
* :mod:`biodb.clinvar` -- ClinVar VCF download + parsing, CLNSIG
  long-tail → 6-class simplification, BED / sites format converters.
* :mod:`biodb.aou_allxall` -- All of Us *All-by-All* PheWAS atlas
  (~3,600 phenotypes × ~414k WGS participants); public-API client for
  per-phenotype gene-burden tables (pLoF/missenseLC/synonymous) with
  concurrent bulk download + consolidated Parquet output.
* :mod:`biodb.orphanet` -- Orphanet rare-disease epidemiology
  (Orphadata ``en_product9_prev``, CC-BY-4.0): per-ORPHAcode prevalence
  band → per-individual K for liability-threshold effect modeling.
* :mod:`biodb.string` -- STRING database physical PPI edges with
  continuous combined-score weights (``download_physical_links``,
  ``load_physical_links``, ``physical_ppi_edges``). The physical
  sub-network is direct binding evidence only — closer to "PPI"
  than the full functional-coupling network in ``protein.links``.

Shared utilities live in :mod:`biodb.utils` (random seeding,
similarity helpers, token counting, ``create_gene_association_matrix``,
``filter_adaptive``).
"""

from biodb import (
    aou_allxall,
    clinvar,
    gprofiler,
    gtr,
    gwas_atlas,
    harmonizome,
    mapping,
    monarch,
    msigdb,
    ols,
    omicspred,
    ontology,
    ontology_owl,
    opentargets,
    opentargets_graphql,
    orphanet,
    pubmed,
    snomed,
    string,
    transform,
    uniprot,
    utils,
)

# A small slice of high-frequency public symbols is re-exported at the
# top level for convenience. The full APIs live on the submodules.
from biodb.clinvar import (
    bed_to_sites,
    df_to_bed,
    df_to_sites,
    download_vcf,
    simplify_annotations,
    vcf_to_df,
)
from biodb.gtr import (
    aggregate_gene_sets as gtr_aggregate_gene_sets,
)
from biodb.gtr import (
    gene_sets as gtr_gene_sets,
)
from biodb.gtr import (
    panel_text as gtr_panel_text,
)
from biodb.gtr import (
    query_test as gtr_query_test,
)
from biodb.gtr import (
    search_tests as gtr_search_tests,
)
from biodb.gtr import (
    to_gmt as gtr_to_gmt,
)
from biodb.harmonizome import (
    download_datasets as harmonizome_download_datasets,
)
from biodb.harmonizome import (
    get_dataset_metadata as harmonizome_get_dataset_metadata,
)
from biodb.harmonizome import (
    get_gmt,
    load_gene_attribute_matrix,
)
from biodb.harmonizome import (
    list_datasets as harmonizome_list_datasets,
)
from biodb.mapping import map_gene_ids
from biodb.monarch import (
    get_gene_associations as monarch_get_gene_associations,
)
from biodb.monarch import (
    read_causal_gene_to_disease_association,
)
from biodb.ontology import (
    get_ancestors,
    get_descendants,
    get_mrca,
    get_ontology,
)
from biodb.opentargets import (
    ensure_cached_shards,
    get_dataset,
    get_gene_associations,
    get_pathways,
    get_targets,
    list_available_versions,
    list_datasets,
)
from biodb.opentargets_graphql import (
    query_disease,
    query_drug,
    query_target,
    query_variant,
)
from biodb.string import (
    physical_ppi_edges as string_physical_ppi_edges,
)
from biodb.uniprot import (
    get_dbxrefs,
    get_features,
    get_sequences,
    query_protein,
)
from biodb.utils import (
    RANDOM_SEED,
    cosine_similarity,
    count_tokens,
    create_gene_association_matrix,
    dot_product_similarity,
    euclidean_similarity,
    filter_adaptive,
    l2_normalize,
    set_random_seed,
)

__version__ = "0.1.0"

__all__ = [
    "RANDOM_SEED",
    "aou_allxall",
    "bed_to_sites",
    "clinvar",
    "cosine_similarity",
    "count_tokens",
    "create_gene_association_matrix",
    "df_to_bed",
    "df_to_sites",
    "dot_product_similarity",
    "download_vcf",
    "ensure_cached_shards",
    "euclidean_similarity",
    "filter_adaptive",
    "get_ancestors",
    "get_dataset",
    "get_dbxrefs",
    "get_descendants",
    "get_features",
    "get_gene_associations",
    "get_gmt",
    "get_mrca",
    "get_ontology",
    "get_pathways",
    "get_sequences",
    "get_targets",
    "gprofiler",
    "gtr",
    "gtr_aggregate_gene_sets",
    "gtr_gene_sets",
    "gtr_panel_text",
    "gtr_query_test",
    "gtr_search_tests",
    "gtr_to_gmt",
    "gwas_atlas",
    "harmonizome",
    "harmonizome_download_datasets",
    "harmonizome_get_dataset_metadata",
    "harmonizome_list_datasets",
    "l2_normalize",
    "list_available_versions",
    "list_datasets",
    "load_gene_attribute_matrix",
    "map_gene_ids",
    "mapping",
    "monarch",
    "monarch_get_gene_associations",
    "msigdb",
    "ols",
    "omicspred",
    "ontology",
    "ontology_owl",
    "opentargets",
    "opentargets_graphql",
    "orphanet",
    "pubmed",
    "query_disease",
    "query_drug",
    "query_protein",
    "query_target",
    "query_variant",
    "read_causal_gene_to_disease_association",
    "set_random_seed",
    "simplify_annotations",
    "snomed",
    "string",
    "string_physical_ppi_edges",
    "transform",
    "uniprot",
    "utils",
    "vcf_to_df",
]
