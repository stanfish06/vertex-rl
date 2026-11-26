"""Serve one 2D solver over newline-delimited JSON on stdin/stdout."""
import argparse
import ctypes
import json
import math
import os
from pathlib import Path
import sys


def number(value, name):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return value


def integer(value, name, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


class Worker:
    def __init__(self, path, threads=1):
        self.threads = integer(threads, "threads", 1, 2**31 - 1)
        self.lib = ctypes.CDLL(str(path.resolve()))
        names = ctypes.POINTER(ctypes.c_char_p)
        values = ctypes.POINTER(ctypes.c_double)
        self.lib.vertex2d_reset.argtypes = (ctypes.c_uint, ctypes.c_int, ctypes.c_int, names, values)
        self.lib.vertex2d_advance.argtypes = (ctypes.c_double,)
        self.lib.vertex2d_export.argtypes = (ctypes.c_int,)
        self.lib.vertex2d_error.restype = self.lib.vertex2d_measure.restype = ctypes.c_char_p
        self.status = "uninitialized"
        self.limit = 0
        self.time = 0

    def state(self):
        return json.loads(self.lib.vertex2d_measure())

    def error(self):
        return self.lib.vertex2d_error().decode()

    def reset(self, request):
        seed = integer(request["seed"], "seed", 0, 2**32 - 1)
        limit = number(request["time_limit"], "time_limit")
        if not 0 < limit <= 3.4028234663852886e38:
            raise ValueError("time_limit must be positive and fit the solver clock")
        parameters = request.get("parameters", {})
        if not isinstance(parameters, dict):
            raise ValueError("parameters must be an object")
        if any("\0" in key for key in parameters):
            raise ValueError("invalid model parameter name")
        n = len(parameters)
        names = (ctypes.c_char_p * n)(*(key.encode() for key in parameters))
        values = (ctypes.c_double * n)(*(number(v, k) for k, v in parameters.items()))
        result = self.lib.vertex2d_reset(seed, self.threads, n, names, values)
        if result == -2:
            raise ValueError(self.error())
        if result:
            self.status = "uninitialized"
            raise ValueError(self.error())
        self.limit = limit
        return self.reply_state()

    def reply_state(self, failure=False):
        state = self.state()
        self.time = state["time"]
        if failure:
            self.status = "solver_error"
            state["message"] = self.error()
        elif not state["measurements"]["cell_count"]:
            self.status = "terminated"
        elif self.time == self.limit:
            self.status = "truncated"
        else:
            self.status = "ok"
        return dict(status=self.status, **state)

    def handle(self, request):
        if not isinstance(request, dict):
            raise ValueError("request must be an object")
        op = request["op"]
        allowed = {"ping": (), "close": (), "measure": (), "export": ("frame_id",),
                   "reset": ("seed", "time_limit", "parameters"), "step": ("deadline",)}
        if not isinstance(op, str) or op not in allowed:
            raise ValueError("unknown operation")
        extra = request.keys() - {"op", *allowed[op]}
        if extra:
            raise ValueError(f"unexpected request fields {sorted(extra)}")
        if op == "ping":
            return {"status": "ready"}
        if op == "close":
            return {"status": "closed"}
        if op == "reset":
            return self.reset(request)
        if self.status == "uninitialized":
            raise ValueError("reset required")
        if op == "measure":
            return dict(status=self.status, **self.state())
        if op == "export":
            frame = integer(request["frame_id"], "frame_id", 1, 2**31 - 1)
            Path("output").mkdir(exist_ok=True)
            if self.lib.vertex2d_export(frame):
                raise ValueError(self.error())
            return {"status": "ok", "frame_id": frame, "time": self.time}
        if self.status != "ok":
            raise ValueError("reset required after episode end or solver failure")
        deadline = number(request["deadline"], "deadline")
        if deadline <= self.time:
            raise ValueError("deadline must be later than the current time")
        failed = self.lib.vertex2d_advance(min(self.limit, deadline))
        return self.reply_state(failure=bool(failed))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", type=Path)
    parser.add_argument("--threads", type=int, default=1, help="solver thread count, retained on reset")
    args = parser.parse_args()
    if not 1 <= args.threads <= 2**31 - 1:
        parser.error("--threads must be a positive integer in the native int range")
    # Legacy native routines print diagnostics to fd 1. Reserve the original
    # stdout pipe exclusively for complete protocol replies.
    with os.fdopen(os.dup(sys.stdout.fileno()), "w", encoding="utf-8") as replies:
        os.dup2(sys.stderr.fileno(), sys.stdout.fileno())
        worker = Worker(args.library, args.threads)
        for line in sys.stdin:
            try:
                reply = worker.handle(json.loads(line))
            except (ValueError, KeyError, OSError, OverflowError) as error:
                reply = {"status": "error", "message": str(error)}
            print(json.dumps(reply, allow_nan=False, separators=(",", ":")), file=replies, flush=True)
            if reply["status"] == "closed":
                break


if __name__ == "__main__":
    main()
