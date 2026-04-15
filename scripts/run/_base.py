"""Shared helpers for run/ entry points."""
import argparse


def base_parser(**kwargs) -> argparse.ArgumentParser:
    """Return a parser pre-loaded with --data-dir and --api-key."""
    p = argparse.ArgumentParser(**kwargs)
    p.add_argument("--data-dir", default="data", metavar="PATH",
                   help="Root data directory (default: data)")
    p.add_argument("--api-key", default=None, metavar="KEY",
                   help="OpenRouter API key (falls back to OPENROUTER_API_KEY env var)")
    p.add_argument("--workers", type=int, default=1, metavar="N",
                   help="Number of parallel LLM workers (default: 1)")
    return p


def make_chemy(args):
    from scripts.chemy import Chemy
    return Chemy(args.data_dir, api_key=args.api_key)
