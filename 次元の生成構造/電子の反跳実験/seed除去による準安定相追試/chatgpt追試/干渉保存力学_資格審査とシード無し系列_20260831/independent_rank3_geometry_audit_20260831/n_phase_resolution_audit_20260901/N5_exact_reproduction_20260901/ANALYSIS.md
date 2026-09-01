# N=5 exact reproduction audit

Authoritative pass2_run.py and original_engine.py were copied byte-for-byte from Drive. Authoritative parent_v.npz was copied byte-for-byte. pass2_run.py was executed unchanged as `python program/pass2_run.py hm_N5`, therefore N=5, L=124, STEPS=40000.

- parent == original state[0]: True
- parent == current state[0]: True
- original state[0] == current state[0]: True
- first state mismatch step: 1

|   step |   original_Hperp_frac |   current_Hperp_frac |   max_abs_state_diff |
|-------:|----------------------:|---------------------:|---------------------:|
|      0 |           4.75516e-69 |          4.75516e-69 |          0           |
|      1 |           4.92381e-31 |          2.11626e-31 |          4.37171e-16 |
|     25 |           6.17753e-30 |          7.13882e-30 |          1.44336e-15 |
|    100 |           6.43597e-29 |          3.80335e-29 |          5.30662e-15 |
|    500 |           6.6247e-27  |          2.79298e-27 |          6.00042e-14 |
|   5000 |           2.97765e-25 |          3.39661e-25 |          6.58853e-13 |
|  10000 |           3.07611e-24 |          1.31088e-24 |          6.46711e-12 |
|  20000 |           1.14625e-23 |          2.96601e-23 |          9.6837e-09  |
|  30000 |           4.93804e-22 |          5.39894e-22 |          1.93976e-05 |
|  35000 |           1.55043e-21 |          1.75161e-21 |          0.000964524 |
|  40000 |           1.98011e-23 |          4.02538e-21 |          0.241047    |

Original final Hperp/H = 1.98011136452219414e-23; current = 4.02538104860422272e-21. The same exact program and exact initial binary64 array diverge from step 1, establishing environment-dependent floating-point/linear-algebra trajectory selection.
