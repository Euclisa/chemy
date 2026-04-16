import os


class RawReactionsLayout:
    def __init__(self, data_dir: str):
        self._base = os.path.join(data_dir, 'raw_reactions')

    def _preset_dir(self, preset_name: str) -> str:
        return os.path.join(self._base, preset_name)

    def _list_runs(self, preset_name: str) -> list:
        dir_ = self._preset_dir(preset_name)
        if not os.path.isdir(dir_):
            return []
        runs = []
        for fn in os.listdir(dir_):
            if fn.startswith('raw_') and fn.endswith('.jsonl'):
                try:
                    runs.append(int(fn[4:-6]))
                except ValueError:
                    pass
        return sorted(runs)

    def raw(self, preset_name: str, run=None) -> str:
        """Resolve raw file path for the given preset and run.

        run=None  — continue latest run (raw_1.jsonl if no runs exist yet).
        run='new' — create a new run (next sequential number).
        run=<int> — use that specific existing run; raises ValueError if absent.
        """
        dir_ = self._preset_dir(preset_name)
        os.makedirs(dir_, exist_ok=True)

        existing = self._list_runs(preset_name)

        if run is None:
            n = max(existing) if existing else 1
            return os.path.join(dir_, f'raw_{n}.jsonl')
        elif run == 'new':
            n = max(existing, default=0) + 1
            return os.path.join(dir_, f'raw_{n}.jsonl')
        else:
            n = int(run)
            if n not in existing:
                raise ValueError(
                    f"Run {n} does not exist for preset '{preset_name}'. "
                    f"Available runs: {existing}"
                )
            return os.path.join(dir_, f'raw_{n}.jsonl')

    def raw_all(self, preset_name: str) -> list:
        """Return paths of all raw run files for a preset, in run order."""
        dir_ = self._preset_dir(preset_name)
        return [
            os.path.join(dir_, f'raw_{n}.jsonl')
            for n in self._list_runs(preset_name)
        ]

    def verdict(self) -> str:
        """Path to the global verdict vault directory."""
        path = os.path.join(self._base, 'verdict')
        os.makedirs(path, exist_ok=True)
        return path

    def repairs(self) -> str:
        os.makedirs(self._base, exist_ok=True)
        return os.path.join(self._base, 'repairs.jsonl')
