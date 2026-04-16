import json

from .parallel_runner import iter_parallel


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

    if store is not None:
        for res in iter_parallel(
            llm_client,
            entries,
            routine,
            logger,
            max_workers=max_workers,
            routine_args=routine_args,
            batch_size=batch_size,
            description=description,
        ):
            items = res if isinstance(res, list) else [res]
            for item in items:
                store.stream_append(item, out_fn)
    else:
        with open(out_fn, 'a') as f_out:
            for res in iter_parallel(
                llm_client,
                entries,
                routine,
                logger,
                max_workers=max_workers,
                routine_args=routine_args,
                batch_size=batch_size,
                description=description,
            ):
                if isinstance(res, list):
                    for item in res:
                        f_out.write(json.dumps(item) + '\n')
                else:
                    f_out.write(json.dumps(res) + '\n')
                f_out.flush()
