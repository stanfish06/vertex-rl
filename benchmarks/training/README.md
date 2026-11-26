# Running the training benchmarks

Build the 3D solver and install the Python package in a virtual environment. From the repository root, run each benchmark with a new destination directory.

```sh
make -C model vertex3d
python -m pip install -e .
python -m benchmarks.training.run learn /tmp/vertex-learning --binary model/vertex3d
python -m benchmarks.training.run throughput /tmp/vertex-throughput --binary model/vertex3d
```

The learning run writes the training targets, all candidate scores and the frozen policy and random-search parameters before evaluating the held-out seeds. Its `evaluation.json` records each interval’s action, time, reward and episode end reason. The throughput run writes the predetermined actions and each batch’s reset and step observations, durations, wall times, replay comparison and ratio. Both commands write a machine record and return a failing exit status when their acceptance condition fails.

The search budgets, seeds, geometry and acceptance criteria are specified in the [frozen protocol](protocol.md).
