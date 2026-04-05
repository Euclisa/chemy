from tests.conftest import make_reaction


def test_parsed_reactions_and_details_prefer_ord_over_llm(compiler_context):
    compiler = compiler_context.compiler
    store = compiler_context.store

    llm_reaction = make_reaction(
        compiler_context.reactions,
        [compiler_context.water],
        [compiler_context.ethanol],
        source="openai/gpt-oss-120b",
        confidence=0.4,
    )
    ord_reaction = {**llm_reaction, "source": "ord", "confidence": 0.9}
    extra_llm = make_reaction(
        compiler_context.reactions,
        [compiler_context.ethanol],
        [compiler_context.methanol],
        source="qwen/qwen3-235b-a22b",
        confidence=0.6,
    )

    llm_details = {
        "rid": llm_reaction["rid"],
        "source": llm_reaction["source"],
        "description": "llm detail",
    }
    ord_details = {
        "rid": llm_reaction["rid"],
        "source": "ord",
        "description": "ord detail",
    }

    store.write_jsonl([llm_reaction, extra_llm], compiler.reactions_parsed_llm_fn, backup=False)
    store.write_jsonl([ord_reaction], compiler.reactions_parsed_ord_fn, backup=False)
    store.write_jsonl([llm_details], compiler.reactions_details_llm_fn, backup=False)
    store.write_jsonl([ord_details], compiler.reactions_details_ord_fn, backup=False)

    parsed_map = {entry["rid"]: entry for entry in compiler.parsed_reactions}
    details_map = {entry["rid"]: entry for entry in compiler.reactions_details}

    assert parsed_map[llm_reaction["rid"]]["source"] == "ord"
    assert parsed_map[extra_llm["rid"]]["source"] == "qwen/qwen3-235b-a22b"
    assert details_map[llm_reaction["rid"]]["description"] == "ord detail"


def test_write_parsed_reactions_splits_files_and_clears_cache(compiler_context):
    compiler = compiler_context.compiler
    store = compiler_context.store

    compiler.__dict__["parsed_reactions"] = ("cached",)

    reactions = [
        make_reaction(
            compiler_context.reactions,
            [compiler_context.water],
            [compiler_context.ethanol],
            source="ord",
        ),
        make_reaction(
            compiler_context.reactions,
            [compiler_context.ethanol],
            [compiler_context.methanol],
            source="openai/gpt-oss-120b",
        ),
    ]

    compiler.write_parsed_reactions(reactions)

    assert "parsed_reactions" not in compiler.__dict__
    assert [entry["source"] for entry in store.load_jsonl(compiler.reactions_parsed_ord_fn)] == ["ord"]
    assert [
        entry["source"] for entry in store.load_jsonl(compiler.reactions_parsed_llm_fn)
    ] == ["openai/gpt-oss-120b"]


def test_merge_parsed_files_uses_source_priority(compiler_context, tmp_path):
    compiler = compiler_context.compiler
    store = compiler_context.store

    low_priority = make_reaction(
        compiler_context.reactions,
        [compiler_context.water],
        [compiler_context.ethanol],
        source="x-ai/grok-4-fast",
    )
    high_priority = {**low_priority, "source": "ord", "confidence": 0.95}

    low_fn = tmp_path / "low.jsonl"
    high_fn = tmp_path / "high.jsonl"
    out_fn = tmp_path / "merged.jsonl"
    store.write_jsonl([low_priority], str(low_fn), backup=False)
    store.write_jsonl([high_priority], str(high_fn), backup=False)

    compiler.merge_parsed_files(str(out_fn), str(low_fn), str(high_fn))

    merged = store.load_jsonl(str(out_fn))
    assert merged == [high_priority]


def test_get_chems_reactions_occurrence_and_generate_edges(compiler_context):
    compiler = compiler_context.compiler
    store = compiler_context.store

    reaction_one = make_reaction(
        compiler_context.reactions,
        [compiler_context.water],
        [compiler_context.ethanol],
        source="ord",
    )
    reaction_two = make_reaction(
        compiler_context.reactions,
        [compiler_context.water],
        [compiler_context.methanol],
        source="openai/gpt-oss-120b",
    )

    compiler.write_parsed_reactions([reaction_one, reaction_two])

    occurrence = compiler.get_chems_reactions_occurence()

    assert occurrence == {
        compiler_context.water["cid"]: 2,
        compiler_context.ethanol["cid"]: 1,
        compiler_context.methanol["cid"]: 1,
        compiler_context.unknown["cid"]: 0,
    }

    compiler.generate_edges()
    edges = store.load_jsonl(compiler_context.compounds.chems_edges_fn)

    assert {(edge["source"], edge["target"]) for edge in edges} == {
        (compiler_context.water["cid"], compiler_context.ethanol["cid"]),
        (compiler_context.water["cid"], compiler_context.methanol["cid"]),
    }
