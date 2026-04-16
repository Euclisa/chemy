import json
from pathlib import Path

from .helpers import (
    build_arg_parser,
    dedupe_by_sample_id,
    deterministic_sample,
    file_sha256,
    filter_by_confidence,
    load_config,
    read_jsonl,
    resolve_path,
    write_jsonl,
)


def sample_band(entries, band_config, *, seed, salt):
    pool = filter_by_confidence(
        entries,
        band_config["confidence_min"],
        band_config.get("confidence_max"),
    )
    sample = deterministic_sample(pool, band_config["n"], seed, salt)
    return pool, sample


def build_samples(config, *, data_dir, root=None):
    data_dir = Path(data_dir)
    root = Path.cwd() if root is None else Path(root)
    verdict_path = resolve_path(config["verdict_path"], root=root)
    verdicts = read_jsonl(verdict_path)
    deduped = dedupe_by_sample_id(verdicts)
    seed = config["seed"]

    samples_dir = data_dir / "samples"
    manifests_dir = data_dir / "manifests"
    samples_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "seed": seed,
        "verdict_path": str(verdict_path),
        "verdict_sha256": file_sha256(verdict_path),
        "input_rows": len(verdicts),
        "deduped_rows": len(deduped),
        "samples": {},
    }

    outputs = {}
    for name, band_config in config["samples"].items():
        pool, sample = sample_band(deduped, band_config, seed=seed, salt=f"audit:{name}")
        out_path = samples_dir / f"audit_{name}.jsonl"
        write_jsonl(sample, out_path)
        outputs[f"audit_{name}"] = out_path
        manifest["samples"][f"audit_{name}"] = {
            "confidence_min": band_config["confidence_min"],
            "confidence_max": band_config.get("confidence_max"),
            "requested_n": band_config["n"],
            "pool_n": len(pool),
            "sample_n": len(sample),
            "path": str(out_path),
        }

    repair_pool, repair_sample = sample_band(
        deduped,
        config["repair"],
        seed=seed,
        salt="repair:borderline_reject",
    )
    repair_out = samples_dir / "repair_borderline_reject.jsonl"
    write_jsonl(repair_sample, repair_out)
    outputs["repair_borderline_reject"] = repair_out
    manifest["samples"]["repair_borderline_reject"] = {
        "confidence_min": config["repair"]["confidence_min"],
        "confidence_max": config["repair"].get("confidence_max"),
        "requested_n": config["repair"]["n"],
        "pool_n": len(repair_pool),
        "sample_n": len(repair_sample),
        "path": str(repair_out),
    }

    manifest_path = manifests_dir / "sample_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    outputs["manifest"] = manifest_path
    return outputs, manifest


def main(argv=None):
    parser = build_arg_parser("Freeze model-quality experiment samples.")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    _, manifest = build_samples(config, data_dir=args.data_dir)
    print(
        "Frozen samples: "
        + ", ".join(
            f"{name}={details['sample_n']}/{details['pool_n']}"
            for name, details in manifest["samples"].items()
        )
    )


if __name__ == "__main__":
    main()
