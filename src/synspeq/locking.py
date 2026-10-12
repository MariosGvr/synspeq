"""Inter-process exclusive lock (flock; an O_EXCL lock file where flock is unavailable)."""
from __future__ import annotations
import contextlib
import os
import time


@contextlib.contextmanager
def file_lock(lock_path: str, timeout: float = 120.0):
    os.makedirs(os.path.dirname(os.path.abspath(lock_path)), exist_ok=True)
    try:
        import fcntl
    except ImportError:
        fcntl = None
    if fcntl is not None:
        fd = os.open(lock_path, os.O_CREAT | os.O_RDWR)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX)
            yield
        finally:
            try:
                fcntl.flock(fd, fcntl.LOCK_UN)
            finally:
                os.close(fd)
        return
    t0 = time.time()
    while True:
        try:
            fd = os.open(lock_path + ".excl", os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            break
        except FileExistsError:
            if time.time() - t0 > timeout:
                raise TimeoutError(f"could not lock {lock_path}")
            time.sleep(0.05)
    try:
        yield
    finally:
        os.close(fd)
        os.remove(lock_path + ".excl")
