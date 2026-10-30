"""Logging and errors shared by every module."""
import time

T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


class HistoalignError(RuntimeError):
    """A run cannot continue (bad input, failed fit). The message says what to change."""
