from concurrent.futures import ThreadPoolExecutor, as_completed


def _make_jobs(entries, batch_size):
    if batch_size is None:
        return list(entries)
    return [
        entries[i:i + batch_size]
        for i in range(0, len(entries), batch_size)
    ]


def iter_parallel(
    llm_client,
    entries,
    routine,
    logger,
    *,
    max_workers=1,
    routine_args=None,
    batch_size=None,
    description="Processing",
):
    """Yield routine results from a bounded worker pool.

    This is intentionally output-agnostic. Callers decide how to persist each
    result, while shared LLM exception handling and progress tracking stay in
    one place.
    """
    if routine_args is None:
        routine_args = []

    jobs = _make_jobs(entries, batch_size)
    logger.log(f"{description}: submitted {len(jobs)} entries")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(routine, job, *routine_args)
            for job in jobs
        ]
        for future in logger.track(as_completed(futures), description, total=len(futures)):
            result = llm_client.get_future_result(future, executor)
            if result is not None:
                yield result
