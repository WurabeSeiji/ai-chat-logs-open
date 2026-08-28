# Historical annealing arguments: NOT RECOVERED

`anneal_subsets.cpp` itself was recovered from the published package. Its CLI is:

```text
anneal_subsets N k input_csv steps seed
```

However, the commands used on 2026-08-25/26 were not stored in the package, and searches of the available Drive/package/chat-derived records did not recover the historical `steps` and `seed` values for the 21 stochastic cases (N=14 k=7..12; N=15 k=7..13; N=16 k=7..14).

Therefore **no historical parameter value is invented here**.

For future reproducible reruns, `reproduction_anneal_args.csv` defines a new explicit matrix. Those values are deliberately labeled `NEW_REPRODUCTION_NOT_HISTORICAL`; they reproduce the algorithm deterministically for a new run but are not claimed to reproduce the exact historical best subsets.
