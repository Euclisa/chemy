from pathlib import Path

from scripts.infra import JsonlStore


def test_write_and_load_jsonl_round_trip(tmp_path):
    store = JsonlStore()
    path = tmp_path / "entries.jsonl"
    entries = [{"cid": 2, "name": "ethanol"}, {"cid": 1, "name": "water"}]

    store.write_jsonl(entries, str(path), backup=False, no_sort=True)

    assert store.load_jsonl(str(path)) == entries
    assert store.entries_to_map(entries, "cid", "name") == {2: "ethanol", 1: "water"}
    assert store.load_jsonl_map(str(path), "cid", "name") == {2: "ethanol", 1: "water"}


def test_write_jsonl_sorts_and_creates_backup(tmp_path):
    store = JsonlStore()
    path = tmp_path / "sorted.jsonl"
    initial = [{"cid": 3}, {"cid": 1}]
    updated = [{"cid": 2}, {"cid": 1}, {"cid": 3}]
    store.register_sorting(str(path), "cid")

    store.write_jsonl(initial, str(path), backup=False)
    store.write_jsonl(updated, str(path))

    backup_path = Path(f"{path}.backup")
    assert backup_path.exists()
    assert store.load_jsonl(str(backup_path)) == [{"cid": 1}, {"cid": 3}]
    assert store.load_jsonl(str(path)) == [{"cid": 1}, {"cid": 2}, {"cid": 3}]


def test_vault_shards_entries_and_loads_them_back(tmp_path):
    store = JsonlStore()
    store._dir_vault_entries_per_file = 2

    vault_dir = tmp_path / "vault"
    store.register_sorting(str(vault_dir), "cid")
    store.register_vault(str(vault_dir), "chunk_")

    entries = [{"cid": 5}, {"cid": 2}, {"cid": 4}, {"cid": 1}, {"cid": 3}]
    store.write_jsonl(entries, str(vault_dir), backup=False)

    shard_names = sorted(p.name for p in vault_dir.glob("*.jsonl"))
    assert shard_names == ["chunk_1.jsonl", "chunk_2.jsonl", "chunk_3.jsonl"]
    assert store.load_jsonl(str(vault_dir)) == [
        {"cid": 1},
        {"cid": 2},
        {"cid": 3},
        {"cid": 4},
        {"cid": 5},
    ]
