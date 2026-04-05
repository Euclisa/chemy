from concurrent.futures import Future, ThreadPoolExecutor
from types import SimpleNamespace

import pytest

from scripts.infra.llm_client import LLMClient


def _build_completion(content, prompt_tokens=3, completion_tokens=7):
    return SimpleNamespace(
        usage=SimpleNamespace(
            prompt_tokens=prompt_tokens,
            completion_tokens=completion_tokens,
        ),
        choices=[SimpleNamespace(message=SimpleNamespace(content=content))],
    )


def test_constructor_tolerates_missing_openai(monkeypatch):
    monkeypatch.setattr("scripts.infra.llm_client.OpenAI", None)

    client = LLMClient(api_key="test")

    assert client.client is None


def test_restrict_ctt_sets_and_resets_limit():
    client = LLMClient(logger=None)

    assert client._max_ctt is None
    with client.restrict_ctt(12):
        assert client._max_ctt == 12
    assert client._max_ctt is None


def test_fetch_answer_updates_token_counters():
    client = LLMClient(logger=None)
    client.client = SimpleNamespace(
        chat=SimpleNamespace(
            completions=SimpleNamespace(
                create=lambda **kwargs: _build_completion("hello", 11, 5)
            )
        )
    )

    result = client.fetch_answer_str("hi", "model-x")

    assert result == "hello"
    assert client.input_tokens_total == 11
    assert client.completion_tokens_total == 5


def test_fetch_answer_raises_completion_limit_before_request():
    client = LLMClient(logger=None)
    client.client = object()
    client.completion_tokens_total = 10

    with client.restrict_ctt(10):
        with pytest.raises(LLMClient.CompletionTokensLimitReached):
            client._fetch_answer([], "model-x")


def test_get_future_result_returns_none_on_generic_exception():
    errors = []
    client = LLMClient(logger=SimpleNamespace(log_err=errors.append))
    future = Future()
    future.set_exception(RuntimeError("boom"))

    with ThreadPoolExecutor(max_workers=1) as executor:
        result = client.get_future_result(future, executor)

    assert result is None
    assert errors == ["Exception occured: boom"]


def test_get_future_result_reraises_completion_limit():
    client = LLMClient(logger=None)
    future = Future()
    future.set_exception(LLMClient.CompletionTokensLimitReached())

    with ThreadPoolExecutor(max_workers=1) as executor:
        with pytest.raises(LLMClient.CompletionTokensLimitReached):
            client.get_future_result(future, executor)


def test_submit_entries_to_llm_writes_single_items_lists_and_batches(tmp_path):
    client = LLMClient(logger=None)
    out_fn = tmp_path / "llm.jsonl"
    logs = []
    logger = SimpleNamespace(log=logs.append, track=lambda it, *args, **kwargs: it)

    client.submit_entries_to_llm(
        str(out_fn),
        [1, 2, 3],
        2,
        lambda entry: {"value": entry},
        logger,
    )

    written_values = sorted(int(line.split(": ")[1].rstrip("}")) for line in out_fn.read_text().strip().splitlines())
    assert written_values == [1, 2, 3]
    assert logs == ["Generation: submitted 3 entries"]

    client.submit_entries_to_llm(
        str(out_fn),
        [4, 5, 6, 7],
        2,
        lambda batch: [{"value": item} for item in batch],
        logger,
        batch_size=2,
    )

    final_values = sorted(int(line.split(": ")[1].rstrip("}")) for line in out_fn.read_text().strip().splitlines())
    assert final_values == [1, 2, 3, 4, 5, 6, 7]
