# Ordered 2D fate draws

The seeded single-thread reference and the optimized solver at 1, 2 and 4 threads emitted identical traces on the six-unit rectangle. Each trace includes the simulation time, recipient cell slot, event type, random draw and outcome. The 122,962 draws consisted of 92,259 division decisions, 30,681 death decisions and 22 division angles. There were 22 successful division decisions and two successful death decisions. The SHA-256 of each complete CSV, including its header, was `3a7e3f2f51044fc36b76150c9d0d02d38315ff47f4c0b84fffcc3fc7b5130f92`.

The candidate writes the trace when `VERTEX_FATE_TRACE` names a file. Leave this environment variable unset for timing. To reproduce the reference instrumentation from the repository root, apply the patch to an isolated copy of the baseline.

```sh
reference_src="$(mktemp -d)"
git archive c91a7b1 model | tar -x -C "$reference_src"
patch -d "$reference_src" -p1 < benchmarks/serial2d/dev-fate-trace.patch
g++ -std=c++17 -O3 -fopenmp -DVERTEX_ARRAY_CAPACITY=1024 \
    "$reference_src/model/main.c" -lreadline -o "$reference_src/dev"
VERTEX_FATE_TRACE="$reference_src/dev.csv" python3 benchmarks/serial2d/run.py \
    "$reference_src/dev" "$reference_src/serial" --seed 646
sha256sum "$reference_src/dev.csv"
```

For each candidate thread count, start the candidate CLI in a fresh run directory with the same five input files and an empty `output/`. Set `VERTEX_FATE_TRACE` to a separate CSV path, launch with `--seed 646 --reference`, then send these CLI commands on stdin, replacing N with the desired thread count.

```text
set number_of_threads N
set outFreq 0.5
set tMAX 6
run
```

Compare the whole CSV with `cmp` and run `compare.py` on the output directory. The trace preserves the float outcome of each division-angle draw before folding the angle at π/2.

## Force accumulation

One OpenMP team computes area terms, edge tension, vertex–edge contacts and shoulder contacts. A polygon owns its vertices in the 2D model. Area terms preserve the order of the two incident polygon sides. Tension terms follow ascending edge ID. Threads evaluate contiguous ranges of contact rows and record each vertex's individual contributions in increasing row and partner order. After the shoulder pass, one thread owns the final sum for each vertex. It reads rows in thread order and rounds each addition into the float accumulator. The maximum uses an OpenMP max reduction. Energy terms and contact-induced target-size changes retain serial evaluation order.

Each block of 32 potential partners caches its earliest eligible interaction time using the solver's existing per-pair schedule. A block can be skipped while every pair is ineligible. Topology changes invalidate these block times and rebuild ordered edge incidence, including when an ID is reused. Fate decisions, topology changes and the second wall projection with dependent cell-center reads retain the serial ordering.
