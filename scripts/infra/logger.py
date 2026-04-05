import threading
import logging

from rich.logging import RichHandler
from rich.console import Console
from rich.progress import Progress
from rich import traceback

from contextlib import contextmanager


class Logger:

    def __init__(self):
        self.print_lock = threading.Lock()

        self._console = Console()
        traceback.install(console=self._console)

        logging.basicConfig(
            level=logging.ERROR,
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(console=self._console)]
        )
        self.__logger = logging.getLogger("Chemy")
        self.__logger.setLevel(logging.INFO)

        self.__no_warnings = False

    @contextmanager
    def no_warnings(self):
        self.__no_warnings = True
        try:
            yield
        finally:
            self.__no_warnings = False

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

    def track(self, iterator, description, total=None, transient=False, auto_refresh=True):
        if total is None:
            try:
                total = len(iterator)
            except TypeError:
                pass

        with self.progress(transient=transient, auto_refresh=auto_refresh) as prog:
            task = prog.add_task(description=description, total=total)
            for item in iterator:
                yield item
                prog.update(task, advance=1)
                if not auto_refresh:
                    prog.refresh()

    def progress(self, transient=False, auto_refresh=True):
        return Progress(console=self._console, auto_refresh=auto_refresh, transient=transient, refresh_per_second=1)
