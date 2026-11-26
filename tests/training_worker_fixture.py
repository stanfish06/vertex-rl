"""One-shot pipe faults in a deterministic JSON-lines worker for environment tests."""
import json
from pathlib import Path
import sys
import time


def reply(obj):
    print(json.dumps(obj),flush=True)


def state(t,status):
    return {"time":t,"status":status,"cells":[{"id":1,"zone":0,"age":t,"growing":True,
            "target_area":1.,"vertices":[[0,0],[1,0],[1,1],[0,1]]}],"measurements":{"cell_count":1}}


for line in sys.stdin:
    q=json.loads(line)
    with Path("requests.jsonl").open("a") as log:
        print(json.dumps(q),file=log)
    flag=Path(q["op"]+".fault")
    fault=flag.read_text() if flag.exists() else ""
    if fault:
        flag.unlink()
    if fault in ("partial","exit"):
        sys.stdout.write('{"status":')
        sys.stdout.flush()
    if fault in ("partial","stall"):
        time.sleep(15)
    if fault=="exit":
        print("exit diagnostic",file=sys.stderr,flush=True)
        sys.exit(7)
    if fault=="delay":
        time.sleep(.15)
    if q["op"]=="ping":
        reply({"status":"ready"})
    elif q["op"]=="reset":
        limit=q["time_limit"]
        reply(state(0.,"ok"))
    elif q["op"]=="step":
        if fault=="solver_error":
            reply({**state(q["deadline"]/2,"solver_error"),"message":"rejected trial"})
        else:
            reply(state(q["deadline"],"truncated" if q["deadline"]==limit else "ok"))
