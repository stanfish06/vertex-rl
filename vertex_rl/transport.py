"""Persistent JSON-lines workers with one monotonic deadline per exchange."""
import itertools
import json
import math
import os
from pathlib import Path
import selectors
import subprocess
import time


class WorkerExitError(RuntimeError):
    def __init__(self, worker_id, status, stderr):
        self.worker_id, self.status, self.stderr = worker_id, status, stderr
        super().__init__(f"{worker_id} exited before a complete reply (status {status}): {stderr}")


class StepTimeout(TimeoutError):
    pass


class ReadinessTimeout(TimeoutError):
    def __init__(self, pending_ids):
        self.pending_ids = list(pending_ids)
        super().__init__(f"readiness deadline expired; missing ready replies from {self.pending_ids}")


class WorkerProtocolError(RuntimeError):
    pass


def positive_timeout(value):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"timeout must be positive and finite, received {value!r}")
    return float(value)


class Worker:
    ids = itertools.count(1)

    def __init__(self, command, directory, worker_id=None):
        self.command = command
        self.directory = Path(directory)
        self.id = worker_id if worker_id is not None else f"worker-{next(self.ids)}"
        self.process = self.stderr = None
        self.ready = False

    def launch(self):
        self.stderr = (self.directory / "worker.stderr").open("wb")
        try:
            self.process = subprocess.Popen(self.command, cwd=self.directory, stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE, stderr=self.stderr, bufsize=0)
        except OSError as error:
            self.close()
            raise WorkerExitError(self.id, None, str(error)) from error
        os.set_blocking(self.process.stdin.fileno(), False)
        os.set_blocking(self.process.stdout.fileno(), False)

    def exited(self):
        process = self.process
        status = process.poll()
        if status is None:
            # stdout EOF can precede a child's exit by the time required for teardown.
            try:
                status = process.wait(timeout=0.1)
            except subprocess.TimeoutExpired:
                pass
        with (self.directory / "worker.stderr").open("rb") as f:
            f.seek(max(0, f.seek(0, 2) - 8192))
            stderr = f.read().decode(errors="replace")
        return WorkerExitError(self.id, status, stderr)

    def close(self):
        self.ready = False
        if self.process is not None:
            if self.process.poll() is None:
                self.process.kill()
            self.process.wait()
            self.process.stdin.close()
            self.process.stdout.close()
            self.process = None
        if self.stderr is not None:
            self.stderr.close()
            self.stderr = None

    def request(self, payload, timeout=None):
        deadline = None if timeout is None else time.monotonic() + positive_timeout(timeout)
        try:
            return _exchange({self: payload}, deadline, startup=False)[self]
        except BaseException:
            self.close()
            raise


def _exchange(requests, deadline, startup):
    """Interleave nonblocking writes and reads, retaining incomplete lines until newline."""
    replies = {}
    incoming = {w: bytearray() for w in requests}
    outgoing = {w: memoryview((json.dumps(q, allow_nan=False, separators=(",", ":")) + "\n").encode())
                for w, q in requests.items()}
    with selectors.DefaultSelector() as selector:
        for w in requests:
            selector.register(w.process.stdin, selectors.EVENT_WRITE, w)
            selector.register(w.process.stdout, selectors.EVENT_READ, w)
        while len(replies) < len(requests):
            remaining = None if deadline is None else deadline - time.monotonic()
            if remaining is not None and remaining <= 0:
                pending = [w.id for w in requests if w not in replies]
                if startup:
                    raise ReadinessTimeout(pending)
                raise StepTimeout(f"{pending[0]} step deadline expired before a complete reply")
            events = selector.select(0.1 if remaining is None else min(remaining, 0.1))
            for key, event in events:
                w = key.data
                if event == selectors.EVENT_WRITE:
                    try:
                        count = os.write(key.fd, outgoing[w])
                    except BlockingIOError:
                        continue
                    except BrokenPipeError as error:
                        raise w.exited() from error
                    outgoing[w] = outgoing[w][count:]
                    if not outgoing[w]:
                        selector.unregister(key.fileobj)
                else:
                    try:
                        part = os.read(key.fd, 65536)
                    except BlockingIOError:
                        continue
                    if not part:
                        raise w.exited()
                    incoming[w].extend(part)
                    if b"\n" in part:
                        line, _, trailing = incoming[w].partition(b"\n")
                        try:
                            if trailing:
                                raise ValueError("multiple or partial extra replies")
                            reply = json.loads(line)
                            if not isinstance(reply, dict) or not isinstance(reply.get("status"), str):
                                raise ValueError("reply must be an object containing status")
                            if startup and reply["status"] != "ready":
                                raise ValueError(f"expected ready, received {reply}")
                        except (ValueError, UnicodeError) as error:
                            raise WorkerProtocolError(f"{w.id}: {error}") from error
                        replies[w] = reply
                        selector.unregister(key.fileobj)
            # A descendant retaining an inherited pipe must not hide worker death.
            if not events:
                for w in requests:
                    if w not in replies and w.process.poll() is not None:
                        raise w.exited()
    return replies


def start_batch(envs, timeout=30):
    """Start every private worker and await all pings under one shared wall-clock limit."""
    envs = list(envs)
    for env in envs:
        if env._closed:
            raise RuntimeError("cannot start a closed environment")
        if env._state is not None and (not env._worker.ready or env._worker.process.poll() is not None):
            raise RuntimeError("explicit reset required to replace an episode's worker")
    _start_workers([env._worker for env in envs], timeout)


def _start_workers(workers, timeout):
    timeout = positive_timeout(timeout)
    if len({w.id for w in workers}) != len(workers):
        raise ValueError("batch workers must have distinct IDs")
    if any(w.process is not None and w.process.poll() is None and not w.ready for w in workers):
        raise RuntimeError("worker startup already in progress")
    deadline = time.monotonic() + timeout
    fresh = [w for w in workers if not w.ready or w.process.poll() is not None]
    try:
        for w in fresh:
            if time.monotonic() >= deadline:
                raise ReadinessTimeout([v.id for v in fresh])
            w.close()
            w.launch()
        _exchange({w: {"op": "ping"} for w in fresh}, deadline, startup=True)
        for w in fresh:
            w.ready = True
    except BaseException:
        for w in workers:
            w.close()
        raise
