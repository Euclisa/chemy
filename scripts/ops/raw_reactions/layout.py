import os


class RawReactionsLayout:
    def __init__(self, data_dir: str):
        self._base = os.path.join(data_dir, 'raw_reactions')

    def raw(self, preset_name: str) -> str:
        dir_ = os.path.join(self._base, preset_name)
        os.makedirs(dir_, exist_ok=True)
        return os.path.join(dir_, 'raw.jsonl')

    def verdict(self, preset_name: str) -> str:
        dir_ = os.path.join(self._base, preset_name)
        os.makedirs(dir_, exist_ok=True)
        return os.path.join(dir_, 'verdict.jsonl')
