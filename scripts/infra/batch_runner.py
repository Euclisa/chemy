import json
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_batch(
    llm_client,
    entries,
    routine,
    out_fn,
    logger,
    max_workers=1,
    routine_args=None,
    batch_size=None,
    description="Processing",
    store=None,
):
    if routine_args is None:
        routine_args = []

    logger.log(f"{description}: submitted {len(entries)} entries")

    def make_futures(executor):
        if batch_size is None:
            return [executor.submit(routine, entry, *routine_args) for entry in entries]
        return [
            executor.submit(routine, entries[i:i + batch_size], *routine_args)
            for i in range(0, len(entries), batch_size)
        ]

    if store is not None:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = make_futures(executor)
            for future in logger.track(as_completed(futures), description, total=len(futures)):
                res = llm_client.get_future_result(future, executor)
                if res is None:
                    continue
                items = res if isinstance(res, list) else [res]
                for item in items:
                    store.stream_append(item, out_fn)
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor, open(out_fn, 'a') as f_out:
            futures = make_futures(executor)
            for future in logger.track(as_completed(futures), description, total=len(futures)):
                res = llm_client.get_future_result(future, executor)
                if res is None:
                    continue
                if isinstance(res, list):
                    for item in res:
                        f_out.write(json.dumps(item) + '\n')
                else:
                    f_out.write(json.dumps(res) + '\n')
                f_out.flush()
