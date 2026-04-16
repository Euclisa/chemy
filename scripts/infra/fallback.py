def with_fallback(fn, models, logger=None):
    """Try fn(model) for each model, returning the first non-None result."""
    for model in models:
        result = fn(model)
        if result is not None:
            return result
        if logger is not None:
            logger.log_warn(f"Falling back from '{model}'")
    return None
