import threading
import os
import json
import shutil
import logging

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

        self._console = Console()
        traceback.install(console=self._console)

        logging.basicConfig(
            level="INFO",
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(console=self._console)]
        )

        self.__logger = logging.getLogger("ChemsDB")

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

        with open(filename) as f:
            content = f.read().strip()
            if not content:
                return []

            return [json.loads(x) for x in content.split('\n')]
    
    def _write_jsonl(self, entries, filename, backup=True, suppress_warnings=False):
        staged_entries = entries

        if filename in self._file_sorting_prefs:
            sorting_prefs = self._file_sorting_prefs[filename]
            if sorting_prefs is not None:
                if isinstance(sorting_prefs, tuple) and len(sorting_prefs) == 2:
                    sorting_field = sorting_prefs[0]
                    sorting_reverse = sorting_prefs[1]
                elif isinstance(sorting_prefs, str):
                    sorting_field = sorting_prefs
                    sorting_reverse = False
                else:
                    raise Exception(f"Invalid format of sorting preferences for file '{filename}': {str(sorting_prefs)}")

                staged_entries = sorted(entries, key=lambda x: x[sorting_field], reverse=sorting_reverse)

        elif not suppress_warnings:
            self.log_warn(f"Writing to '{filename}' without sorting")


        if os.path.exists(filename) and backup:
            shutil.copy(filename, f"{filename}.backup")

        with open(filename, 'w') as f:
            for entry in staged_entries:
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

    
