# Model Quality Experiment

This experiment estimates where stronger models are worth using in raw reaction
validation and repair. It keeps production pipeline behavior unchanged and writes
all experiment state under `experiments/model_quality/data/`.

## Workflow

Freeze deterministic samples:

```bash
python -m experiments.model_quality.sample
```

Run the fixed-round validation audit with the configured strong model:

```bash
python -m experiments.model_quality.run_audit
```

Run repair generation and validation routes:

```bash
python -m experiments.model_quality.run_repairs
```

Summarize frozen outputs:

```bash
python -m experiments.model_quality.analyze
```

## Samples

The default config samples:

- `high_accept`: 80 reactions with original confidence `>= 0.8`
- `borderline_accept`: 80 reactions with confidence `0.4 <= c < 0.6`
- `borderline_reject`: 80 reactions with confidence `0.2 <= c < 0.4`
- `low_reject`: 40 reactions with confidence `0.0 <= c < 0.1`
- `repair_borderline_reject`: 80 reactions with confidence `0.2 <= c < 0.4`

Rows are deduplicated by `rid` when present, otherwise by a stable hash of the
reaction object.

## Repair Routes

The first pass intentionally runs only the representative routes:

- original reaction + strong validator
- cheap fix + cheap validator
- cheap fix + strong validator
- strong fix + strong validator

Repairs are generated once and frozen before validation, so verifier comparisons
consume the exact same fixed reaction objects.

## Notes

Audit validation runs a fixed number of rounds instead of using the production
early stopping behavior. This makes confidence values more comparable across
models and confidence bands.
