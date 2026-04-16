# Raw Reaction Model Bakeoff

This experiment compares raw-reaction generation models from scratch. It freezes
the same compound tasks for every candidate model, stratifies tasks by existing
presets and compound complexity, stores raw model outputs, normalizes generated
reactions into candidates, and audits candidates with a configurable strong model.

All experiment state is written under `experiments/model_quality/data/`.
Production `chemy.raw_reactions.fetch/validate/repair/parse` behavior uses the
same shared raw-reaction services, but keeps its existing production output file
shapes.

## Notebook Workflow

The easiest way to run the experiment is the control notebook:

```bash
jupyter notebook experiments/model_quality/notebooks/model_bakeoff.ipynb
```

The notebook centralizes config review, smoke-mode overrides, task freezing,
candidate generation, gold audit, scoring, output previews, and plots. Live model
calls are guarded by explicit switches in the setup cell:

```python
RUN_FREEZE_TASKS = False
RUN_GENERATION = False
RUN_GOLD_AUDIT = False
RUN_SCORE = True
SMOKE_MODE = True
SMOKE_GENERATION_RUN_LIMIT = 2
SMOKE_MAX_WORKERS = 2
```

Start with `SMOKE_MODE = True`, `RUN_FREEZE_TASKS = True`,
`RUN_GENERATION = True`, `RUN_GOLD_AUDIT = True`, and one cheap candidate model.
After the smoke run succeeds, disable smoke mode or edit `config.yaml` for the
full run.

Smoke mode intentionally narrows to one preset, one complexity band, one
candidate model, disables generation self-review, and caps new generation runs so
stage 3 stays interactive. Smoke outputs go to `experiments/model_quality/data_smoke/`;
full-run outputs go to `experiments/model_quality/data/`.

`max_workers` controls parallel model calls for generation and gold-audit batches.
Keep it modest because provider rate limits can make higher values slower or
noisier.

## CLI Workflow

The same steps are available as commands if you want a terminal run.

Freeze deterministic tasks:

```bash
python -m experiments.model_quality.sample
```

Run generation for all configured candidate models:

```bash
python -m experiments.model_quality.run_generation --workers 2
```

For a tiny CLI smoke generation run:

```bash
python -m experiments.model_quality.run_generation --model openai/gpt-oss-120b --limit 2 --workers 2
```

Run the strong-model gold audit:

```bash
python -m experiments.model_quality.run_gold_audit --workers 2
```

Score the outputs:

```bash
python -m experiments.model_quality.score --out experiments/model_quality/data/summary.json
```

## Task Sampling

The default config samples from:

- `wiki_crc_rp`
- `simple_inorganic_wiki_rp`
- `simple_organic_wiki_rp`
- `oxidizers_wiki_rp`
- `corrosive_wiki_rp`
- `top_rare_rp`
- `wiki_uncommon_p`

For each preset, eligible compounds are sorted by `bertz_complexity` and
`heavy_count`, then split into `simple`, `medium`, and `hard` bands. The seed and
`tasks_per_band` setting make the frozen task set reproducible.

## Outputs

- `tasks.jsonl`: frozen prompts and task metadata
- `generations.jsonl`: raw model responses, parsed reactions, usage, latency
- `candidates.jsonl`: normalized reaction candidates with parser/evaluator status
- `gold_audits.jsonl`: strong-model audit votes and aggregate confidence
- `summary.json`: optional scored metrics for notebooks and reports
