import threading
import os
import json
import shutil
import logging
import glob

from rich.logging import RichHandler
from rich.console import Console
from rich import traceback

from contextlib import contextmanager


class ChemsDB:

    def __init__(self, data_dir):
        self.print_lock = threading.Lock()

        self.data_dir = data_dir
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
        
        self.structures_dir = os.path.join(self.data_dir, 'structures')
        if not os.path.exists(self.structures_dir):
            os.makedirs(self.structures_dir)

        self._file_sorting_prefs = dict()
        self._dir_vault_prefs = dict()

        self._dir_vault_entries_per_file = 5000

        self._console = Console()
        traceback.install(console=self._console)

        logging.basicConfig(
            level=logging.ERROR,
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(console=self._console)]
        )
        self.__logger = logging.getLogger("ChemsDB")
        self.__logger.setLevel(logging.INFO)

        self.__no_warnings = False
    

    @contextmanager
    def no_warnings(self):
        self.__no_warnings = True
        try:
            yield
        finally:
            self.__no_warnings = False


    def _load_jsonl(self, filename):
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
                key=lambda x: int(x[len(prefix):-len('.jsonl')])  # strip prefix and '.jsonl'
            )

            all_entries = []
            for fn in files:
                full_path = os.path.join(filename, fn)
                all_entries.extend(self._load_jsonl(full_path))

            return all_entries

        else:
            with open(filename) as f:
                content = f.read().strip()
                if not content:
                    return []

                return [json.loads(x) for x in content.split('\n')]
    

    def _apply_sorting_prefs(self, entries, sorting_prefs):
        sorting_prefs = sorting_prefs
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
        # 1. Remove all existing .backup files
        for fn in glob.glob(os.path.join(dir_path, "*.backup")):
            try:
                os.remove(fn)
            except Exception as e:
                self.log_warn(f"Failed to remove old backup '{fn}': {e}")

        # 2. Create fresh backups of all .jsonl files
        for fn in glob.glob(os.path.join(dir_path, "*.jsonl")):
            backup_fn = fn + ".backup"
            try:
                shutil.copy(fn, backup_fn)
            except Exception as e:
                self.log_warn(f"Failed to create backup for '{fn}': {e}")

        # 3. Remove all original .jsonl files
        for fn in glob.glob(os.path.join(dir_path, "*.jsonl")):
            try:
                os.remove(fn)
            except Exception as e:
                self.log_warn(f"Failed to remove old file '{fn}': {e}")

    
    def _write_jsonl(self, entries, filename, backup=True, no_sort=False):

        if filename in self._file_sorting_prefs:
            entries = self._apply_sorting_prefs(entries, self._file_sorting_prefs[filename])
        elif not no_sort:
            self.log_warn(f"Writing to '{filename}' without sorting")

        if os.path.isdir(filename):
            if filename not in self._dir_vault_prefs:
                raise ValueError(f"Directory '{filename}' has no preferences set")
            
            prefix = self._dir_vault_prefs[filename]

            self.__prepare_dirs(filename)

            for batch_idx in range(0, len(entries), self._dir_vault_entries_per_file):
                batch = entries[batch_idx:batch_idx + self._dir_vault_entries_per_file]
                file_num = batch_idx // self._dir_vault_entries_per_file + 1
                fn = os.path.join(filename, f"{prefix}{file_num}.jsonl")
                self._write_jsonl(batch, fn, backup=False, no_sort=True)

        else:
            if os.path.exists(filename) and backup:
                shutil.copy(filename, f"{filename}.backup")

            with open(filename, 'w') as f:
                for entry in entries:
                    f.write(json.dumps(entry) + '\n')
    

    def print(self, message=""):
        with self.print_lock:
            self._console.print(message)


    def log(self, message=""):
        with self.print_lock:
            self.__logger.info(message)
    
    def log_warn(self, message):
        if self.__no_warnings:
            return

        with self.print_lock:
            self.__logger.warning(message)
    

    def log_err(self, message):
        with self.print_lock:
            self.__logger.error(message)

    
