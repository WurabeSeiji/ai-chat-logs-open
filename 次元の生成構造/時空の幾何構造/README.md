# n=64 round12 brute-force experiment

This bundle contains the exact current experiment and a separate plotting script.

## Files

- `n64_round12_simulation.py`
  - Runs the brute-force simulation.
  - Default: `n=64`, `max_step=n*2=128`.
  - Writes `n64_round12_results.csv`.
  - Uses the current WRITE deduplication key:
    - `round(value.real, 12)`
    - `round(value.imag, 12)`

- `n64_round12_plot.py`
  - Reads the CSV.
  - Produces PNG and SVG figures.
  - It can also plot a partial CSV if the simulation is interrupted.

## Run

```bash
python3 n64_round12_simulation.py
python3 n64_round12_plot.py
```

## Current dynamics

```text
Wave(harmonic, value)

rotate:
    value <- omega**harmonic * value

pair_read:
    harmonic <- a.harmonic + b.harmonic
    value    <- a.value * b.value

step:
    current snapshot only
    -> rotate each wave
    -> rotate every unordered pair product, including self-pairs
    -> keep only the first candidate for each rounded complex-value key
```

No Fourier aggregation or modulo-state compression is used.

## Important computational note

The pair loop is O(M^2) in the number of retained states. In the previous run,
step 53 had 16,530 retained states, making the next all-pairs step roughly
136 million candidate pairs. This is intentionally not optimized away here,
because the goal is to hand over the current brute-force experiment unchanged.
