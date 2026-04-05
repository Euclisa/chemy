import os
import json
import shutil
import glob


class JsonlStore:

    def __init__(self):
        self._file_sorting_prefs = dict()
        self._dir_vault_prefs = dict()
        self._dir_vault_entries_per_file = 5000

    def register_sorting(self, path, prefs):
        self._file_sorting_prefs[path] = prefs

    def register_vault(self, path, prefix):
        self._dir_vault_prefs[path] = prefix

    def load_jsonl(self, filename):
        if not os.path.exists(filename):
            return []

        if os.path.isdir(filename):
            if filename not in self._dir_vault_prefs:
                raise ValueError(f"Directory '{filename}' has no preferences set")

            prefix = self._dir_vault_prefs[filename]

            files = sorted(
                (
                    f for f in os.listdir(filename)
                    if f.endswith('.jsonl') and f.startswith(prefix)
                ),
                key=lambda x: int(x[len(prefix):-len('.jsonl')])
            )

            all_entries = []
            for fn in files:
                full_path = os.path.join(filename, fn)
                all_entries.extend(self.load_jsonl(full_path))

            return all_entries

        else:
            with open(filename) as f:
                content = f.read().strip()
                if not content:
                    return []

                return [json.loads(x) for x in content.split('\n')]

    def entries_to_map(self, entries, key, value):
        return {entry[key]: entry[value] for entry in entries}

    def load_jsonl_map(self, filename, key, value):
        entries = self.load_jsonl(filename)
        return self.entries_to_map(entries, key, value)

    def _apply_sorting_prefs(self, entries, sorting_prefs):
        if sorting_prefs is not None:
            if isinstance(sorting_prefs, tuple) and len(sorting_prefs) == 2:
                sorting_field = sorting_prefs[0]
                sorting_reverse = sorting_prefs[1]
            elif isinstance(sorting_prefs, str):
                sorting_field = sorting_prefs
                sorting_reverse = False
            else:
                raise Exception(f"Invalid format of sorting preferences: {str(sorting_prefs)}")

            entries = sorted(entries, key=lambda x: x[sorting_field], reverse=sorting_reverse)

        return entries

    def __prepare_dirs(self, dir_path):
        for fn in glob.glob(os.path.join(dir_path, "*.jsonl")):
            try:
                os.remove(fn)
            except Exception:
                pass

    def __backup_path(self, path):
        backup = f"{path}.backup"
        if os.path.exists(backup):
            if os.path.isdir(backup):
                shutil.rmtree(backup)
            else:
                os.remove(backup)

        if os.path.isdir(path):
            shutil.copytree(path, backup)
        else:
            shutil.copy(path, backup)

    def write_jsonl(self, entries, filename, backup=True, no_sort=False):

        if filename in self._file_sorting_prefs:
            entries = self._apply_sorting_prefs(entries, self._file_sorting_prefs[filename])
        elif not no_sort:
            pass  # caller is responsible for logging warnings

        if os.path.isdir(filename) and filename not in self._dir_vault_prefs:
            raise ValueError(f"Directory '{filename}' has no preferences set")

        if os.path.exists(filename) and backup:
            self.__backup_path(filename)

        if not os.path.exists(filename) and filename in self._dir_vault_prefs:
            os.makedirs(filename, exist_ok=True)

        if os.path.isdir(filename):
            prefix = self._dir_vault_prefs[filename]

            self.__prepare_dirs(filename)

            for batch_idx in range(0, len(entries), self._dir_vault_entries_per_file):
                batch = entries[batch_idx:batch_idx + self._dir_vault_entries_per_file]
                file_num = batch_idx // self._dir_vault_entries_per_file + 1
                fn = os.path.join(filename, f"{prefix}{file_num}.jsonl")
                self.write_jsonl(batch, fn, backup=False, no_sort=True)

        else:
            with open(filename, 'w') as f:
                for entry in entries:
                    f.write(json.dumps(entry) + '\n')

    def map_to_entries(self, in_map, key, value):
        return [{key: k, value: v} for k, v in in_map.items()]

    def write_jsonl_map(self, in_map, key, value, filename, backup=True, no_sort=False):
        entries = self.map_to_entries(in_map, key, value)
        self.write_jsonl(entries, filename, backup=backup, no_sort=no_sort)
