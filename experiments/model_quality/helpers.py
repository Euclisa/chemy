import argparse
import json
from copy import deepcopy
from pathlib import Path

from scripts.infra.store import JsonlStore


EXPERIMENT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = EXPERIMENT_DIR / "config.yaml"
DEFAULT_DATA_DIR = EXPERIMENT_DIR / "data"

DEFAULT_CONFIG = {
    "seed": 1729,
    "acceptance_threshold": 0.4,
    "batch_size": 5,
    "max_workers": 2,
    "chemy_data_dir": "data",
    "strong_model": "google/gemini-2.5-pro",
    "candidate_models": [
        "openai/gpt-oss-120b",
        "qwen/qwen3-235b-a22b",
        "x-ai/grok-4-fast",
    ],
    "bakeoff_presets": [
        "wiki_crc_rp",
        "simple_inorganic_wiki_rp",
        "simple_organic_wiki_rp",
        "oxidizers_wiki_rp",
        "corrosive_wiki_rp",
        "top_rare_rp",
        "wiki_uncommon_p",
    ],
    "tasks_per_band": 5,
    "complexity_bands": ["simple", "medium", "hard"],
    "generation_self_review": True,
    "gold_rounds": 3,
}


def _parse_scalar(value):
    value = value.strip()
    if value in {"", "null", "None", "~"}:
        return None
    if value.startswith("[") and value.endswith("]"):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            return [
                item.strip().strip('"').strip("'")
                for item in value[1:-1].split(",")
                if item.strip()
            ]
    if value in {"true", "True"}:
        return True
    if value in {"false", "False"}:
        return False
    if (
        (value.startswith('"') and value.endswith('"'))
        or (value.startswith("'") and value.endswith("'"))
    ):
        return value[1:-1]
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _load_simple_yaml(path):
    """Load the small YAML subset used by config.yaml without adding PyYAML."""
    root = {}
    stack = [(-1, root)]
    with open(path) as f:
        for raw_line in f:
            if not raw_line.strip() or raw_line.lstrip().startswith("#"):
                continue
            indent = len(raw_line) - len(raw_line.lstrip(" "))
            line = raw_line.strip()
            if ":" not in line:
                raise ValueError(f"Unsupported config line: {raw_line.rstrip()}")
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()

            while indent <= stack[-1][0]:
                stack.pop()

            parent = stack[-1][1]
            if value == "":
                node = {}
                parent[key] = node
                stack.append((indent, node))
            else:
                parent[key] = _parse_scalar(value)
    return root


def deep_merge(base, override):
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(path=DEFAULT_CONFIG_PATH):
    if path is None:
        return deepcopy(DEFAULT_CONFIG)
    path = Path(path)
    if not path.exists():
        return deepcopy(DEFAULT_CONFIG)
    loaded = _load_simple_yaml(path)
    return deep_merge(DEFAULT_CONFIG, loaded)


def resolve_path(path, *, root=None):
    path = Path(path)
    if path.is_absolute():
        return path
    root = Path.cwd() if root is None else Path(root)
    return root / path


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_jsonl(path):
    return JsonlStore().load_jsonl(str(path))


def write_jsonl(entries, path):
    ensure_parent(path)
    with open(path, "w") as f:
        for entry in entries:
            f.write(json.dumps(entry) + "\n")


def append_jsonl(entry, path):
    ensure_parent(path)
    with open(path, "a") as f:
        f.write(json.dumps(entry) + "\n")


def add_common_args(parser):
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG_PATH),
        help="Path to model-quality experiment config.",
    )
    parser.add_argument(
        "--data-dir",
        default=str(DEFAULT_DATA_DIR),
        help="Directory for frozen tasks and experiment outputs.",
    )
    return parser


def build_arg_parser(description):
    return add_common_args(argparse.ArgumentParser(description=description))
