import json
import random
from pathlib import Path

from scripts.chemy import Chemy
from scripts.ops.raw_reactions.presets import get_preset
from scripts.ops.raw_reactions.services import RawReactionTaskBuilder

from .helpers import build_arg_parser, load_config, resolve_path, write_jsonl


def _complexity_key(chem):
    return (
        chem.get("bertz_complexity") or 0,
        chem.get("heavy_count") or 0,
        chem.get("cmpdname") or "",
        chem.get("cid") or 0,
    )


def split_complexity_bands(chems):
    ordered = sorted(chems, key=_complexity_key)
    n = len(ordered)
    bands = {"simple": [], "medium": [], "hard": []}
    if n == 0:
        return bands
    for idx, chem in enumerate(ordered):
        frac = idx / n
        if frac < 1 / 3:
            band = "simple"
        elif frac < 2 / 3:
            band = "medium"
        else:
            band = "hard"
        bands[band].append(chem)
    return bands


def _sample(entries, n, seed, salt):
    entries = sorted(entries, key=lambda chem: (chem.get("cid"), chem.get("cmpdname")))
    rng = random.Random(f"{seed}:{salt}")
    if n >= len(entries):
        return list(entries)
    return rng.sample(entries, n)


def build_tasks(config, *, data_dir, root=None):
    root = Path.cwd() if root is None else Path(root)
    chemy_data_dir = resolve_path(config["chemy_data_dir"], root=root)
    chemy = Chemy(str(chemy_data_dir))
    task_builder = RawReactionTaskBuilder(chemy.reactions)
    out_dir = Path(data_dir)
    manifests_dir = out_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    manifest = {
        "seed": config["seed"],
        "chemy_data_dir": str(chemy_data_dir),
        "tasks_per_band": config["tasks_per_band"],
        "presets": {},
    }

    for preset_name in config["bakeoff_presets"]:
        preset = get_preset(preset_name)
        criteria = preset.build_criteria(chemy.compounds)
        eligible = [chem for chem in chemy.compounds.chems if criteria(chem)]
        bands = split_complexity_bands(eligible)
        manifest["presets"][preset_name] = {
            "eligible_n": len(eligible),
            "bands": {},
        }
        for band_name in config["complexity_bands"]:
            pool = bands.get(band_name, [])
            selected = _sample(
                pool,
                config["tasks_per_band"],
                config["seed"],
                f"{preset_name}:{band_name}",
            )
            manifest["presets"][preset_name]["bands"][band_name] = {
                "pool_n": len(pool),
                "task_n": len(selected),
            }
            for chem in selected:
                tasks.append(
                    task_builder.build_task(
                        chem,
                        preset,
                        complexity_band=band_name,
                        mode="cold_start",
                    ).to_dict()
                )

    tasks_path = out_dir / "tasks.jsonl"
    write_jsonl(tasks, tasks_path)
    manifest["task_n"] = len(tasks)
    manifest["tasks_path"] = str(tasks_path)
    manifest_path = manifests_dir / "task_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return tasks_path, manifest


def main(argv=None):
    parser = build_arg_parser("Freeze fresh-generation model bakeoff tasks.")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    _, manifest = build_tasks(config, data_dir=args.data_dir)
    print(
        f"Frozen {manifest['task_n']} tasks across "
        f"{len(manifest['presets'])} presets."
    )


if __name__ == "__main__":
    main()
