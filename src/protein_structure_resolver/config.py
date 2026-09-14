"""외부 서비스 주소와 파이프라인 기본 설정."""

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
RCSB_FILE_URL = "https://files.rcsb.org/download"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction"

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
MIN_RCSB_SEQUENCE_LENGTH = 25
PAGE_SIZE = 100
GRAPHQL_BATCH_SIZE = 100
REQUEST_TIMEOUT = 30
DOWNLOAD_TIMEOUT = 120
MAX_ATTEMPTS = 4
RETRYABLE_STATUS = frozenset({408, 429, 500, 502, 503, 504})
RESOLVER_VERSION = "0.4.0"
USER_AGENT = f"protein-structure-resolver/{RESOLVER_VERSION}"
METADATA_SCHEMA_VERSION = "1.0"
CACHE_FORMAT_VERSION = "2"
RANKING_RULE_VERSION = "1.0"
ALPHAFOLD_RANKING_RULE_VERSION = "1.0"
ESMFOLD2_MODEL = "esmfold2-fast-2026-05"
RESOLUTION_ORDER = (
    "rcsb_experimental",
    "alphafold_db",
    "esmfold2",
)

EXPERIMENTAL_ENTITY_QUERY = """
query($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    entity_poly { pdbx_seq_one_letter_code_can }
    rcsb_polymer_entity_container_identifiers {
      asym_ids
      auth_asym_ids
    }
    entry {
      rcsb_id
      rcsb_entry_info {
        experimental_method
        resolution_combined
      }
      rcsb_accession_info { initial_release_date }
      exptl { method }
      refine { ls_R_factor_R_free }
      pdbx_vrpt_summary_geometry {
        clashscore
        percent_ramachandran_outliers
      }
    }
    polymer_entity_instances {
      rcsb_id
      rcsb_polymer_instance_info { modeled_residue_count }
      rcsb_polymer_entity_instance_container_identifiers {
        asym_id
        auth_asym_id
      }
    }
  }
}
"""

COMPUTATIONAL_ENTITY_QUERY = """
query($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    entity_poly { pdbx_seq_one_letter_code_can }
    rcsb_polymer_entity_container_identifiers {
      reference_sequence_identifiers {
        database_accession
        database_name
      }
    }
    entry {
      rcsb_id
      rcsb_comp_model_provenance {
        entry_id
        source_db
      }
    }
  }
}
"""

SELECTED_ENTITY_FEATURE_QUERY = """
query($entry_id: String!, $entity_id: String!) {
  polymer_entity(entry_id: $entry_id, entity_id: $entity_id) {
    rcsb_polymer_entity_feature {
      type
      feature_id
      name
      description
      feature_positions {
        beg_seq_id
        end_seq_id
      }
    }
    polymer_entity_instances {
      rcsb_id
      rcsb_polymer_entity_instance_container_identifiers {
        asym_id
        auth_asym_id
      }
      rcsb_polymer_instance_feature {
        type
        feature_positions {
          beg_seq_id
          end_seq_id
        }
      }
    }
  }
}
"""
