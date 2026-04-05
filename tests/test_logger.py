from rich.progress import Progress

from scripts.infra import Logger


def test_no_warnings_suppresses_only_inside_context(monkeypatch):
    logger = Logger()
    captured = []

    monkeypatch.setattr(logger._Logger__logger, "warning", captured.append)

    with logger.no_warnings():
        logger.log_warn("hidden")
    logger.log_warn("visible")

    assert captured == ["visible"]


def test_track_yields_items_for_iterables_without_len():
    logger = Logger()

    items = (value for value in [1, 2, 3])

    assert list(logger.track(items, "demo")) == [1, 2, 3]


def test_progress_returns_rich_progress_instance():
    logger = Logger()

    progress = logger.progress()

    assert isinstance(progress, Progress)
