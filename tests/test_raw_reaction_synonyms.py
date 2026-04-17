import json

from scripts.ops.raw_reactions.synonym_enrichment import (
    RawReactionSynonymEnricher,
    parse_synonym_jsonl,
)
from tests.conftest import make_chem


class FakeLLMClient:
    def __init__(self, responses):
        self.responses = list(responses)
        self.completion_tokens_total = 0

    def fetch_answer_str(self, *args, **kwargs):
        return self.responses.pop(0)

    def get_future_result(self, future, executor):
        return future.result()


def make_enricher(ops_context, responses=None):
    llm_client = FakeLLMClient(responses or [])
    return RawReactionSynonymEnricher(
        str(ops_context.data_dir),
        ops_context.compounds,
        ops_context.store,
        ops_context.logger,
        llm_client,
        ops_context.reaction_llm,
    )


def test_collect_unmapped_names_counts_participant_names(ops_context):
    enricher = make_enricher(ops_context)
    raw_fn = ops_context.raw_reactions.layout.raw("wiki_crc_rp")
    ops_context.store.write_jsonl(
        [
            {
                "cid": ops_context.water["cid"],
                "reactions": [
                    {
                        "reagents": [
                            {"name": "Water", "phase": "l"},
                            {"name": "Grain spirit", "phase": "l"},
                        ],
                        "products": [{"name": "Methanol", "phase": "l"}],
                    },
                    {
                        "reagents": [{"name": "grain spirit", "phase": "l"}],
                        "products": [{"name": "Methanol", "phase": "l"}],
                    },
                ],
            }
        ],
        raw_fn,
        backup=False,
    )

    entries = enricher.collect_unmapped_names("wiki_crc_rp")

    assert entries[0]["norm_name"] == "grainspirit"
    assert entries[0]["count"] == 2
    assert entries[0]["raw_variants"] == [
        {"name": "Grain spirit", "count": 1},
        {"name": "grain spirit", "count": 1},
    ]


def test_parse_synonym_jsonl_groups_by_normalized_input_name(ops_context):
    parsed = parse_synonym_jsonl(
        "\n".join([
            json.dumps({"name": "Grain spirit", "synonyms": ["Ethyl alcohol"]}),
            json.dumps({"name": "unknown", "synonyms": ["Water"]}),
            "not json",
        ]),
        {"grainspirit"},
        ops_context.compounds.normalize_chem_name,
    )

    assert parsed == {"grainspirit": ["Ethyl alcohol"]}


def test_mine_batch_counts_votes_by_normalized_synonym(ops_context):
    responses = [
        json.dumps({"name": "Grain spirit", "synonyms": ["Ethyl alcohol", "ethyl alcohol"]}),
        json.dumps({"name": "Grain spirit", "synonyms": ["ethyl  alcohol"]}),
        json.dumps({"name": "Grain spirit", "synonyms": ["Methanol"]}),
    ]
    enricher = make_enricher(ops_context, responses)
    entry = {
        "name": "Grain spirit",
        "norm_name": "grainspirit",
        "count": 3,
        "raw_variants": [{"name": "Grain spirit", "count": 3}],
        "examples": [],
    }

    audit = enricher._mine_batch([entry], model="test-model", rounds=3, min_votes=2)

    assert audit[0]["status"] == "ready"
    assert audit[0]["approved"] is True
    assert audit[0]["target_cid"] == ops_context.ethanol["cid"]
    assert audit[0]["add_synonym"] == "Grain spirit"
    assert audit[0]["candidate_synonyms"][0]["norm_synonym"] == "ethylalcohol"
    assert audit[0]["candidate_synonyms"][0]["votes"] == 2


def test_apply_adds_approved_original_name_to_matched_compound(ops_context, tmp_path):
    enricher = make_enricher(ops_context)
    audit_fn = tmp_path / "audit.jsonl"
    ops_context.store.write_jsonl(
        [
            {
                "name": "Grain spirit",
                "status": "ready",
                "approved": True,
                "target_cid": ops_context.ethanol["cid"],
                "add_synonym": "Grain spirit",
            }
        ],
        str(audit_fn),
        backup=False,
    )

    summary = enricher.apply(str(audit_fn))

    assert summary["applied"] == 1
    assert "Grain spirit" in ops_context.compounds.cid_chem_map[
        ops_context.ethanol["cid"]
    ]["cmpdsynonym"]


def test_apply_skips_synonym_that_would_create_conflict(ops_context, tmp_path):
    shared = make_chem(99, "Shared Target", "CCN", synonyms=["Shared Name"])
    ops_context.compounds.update_chems(list(ops_context.compounds.chems) + [shared])
    enricher = make_enricher(ops_context)
    audit_fn = tmp_path / "audit.jsonl"
    ops_context.store.write_jsonl(
        [
            {
                "name": "Shared Name",
                "status": "ready",
                "approved": True,
                "target_cid": ops_context.ethanol["cid"],
                "add_synonym": "Shared Name",
            }
        ],
        str(audit_fn),
        backup=False,
    )

    summary = enricher.apply(str(audit_fn))

    assert summary["applied"] == 0
    assert summary["skipped_by_reason"] == {"would_create_conflict": 1}
    assert "Shared Name" not in ops_context.compounds.cid_chem_map[
        ops_context.ethanol["cid"]
    ]["cmpdsynonym"]
