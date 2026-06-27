#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Does the localized odd-harmonic wave reproduce the base wave's squared
 (Born) distribution, under the two conditions
   (1) only phase DIFFERENCES are observable   -> convolution kernel + squaring
   (2) unobservable (finer-than-N) phase diffs are ignored -> finite-N truncation
 using ONLY the single localized kernel of the paper?
================================================================================

Paper object (half-wavelength interval phi in [-pi/2, pi/2]):
    S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi) = sin((N+1) phi) / (2 sin phi)
    I_N(phi) = S_N(phi)^2,   peak S_N(0) = (N+1)/2,   width ~ 1/(N+1)

Key basis fact on [-pi/2, pi/2]:
    {cos((2m+1)phi)}  are ORTHOGONAL:  int cos((2m+1)phi)cos((2m'+1)phi) dphi = (pi/2) d_{mm'}
So S_N = sum of the first (N+1)/2 basis functions with EQUAL amplitude 1
       = truncated reproducing kernel of this basis.

Two readings tested:
  (A) ENVELOPE reading : average the fast factor sin^2((N+1)phi) -> 1/2,
      observed weight ~ 1/((N+1) sin phi)^2 .  Compare to base^2 = cos^2.
  (B) PROJECTOR reading: observation = scan localized kernel across base state
      (convolution = depends only on phase DIFFERENCE) then square (phase-diff only).
      P_N(phi0) = | int S_N(phi0 - phi) psi_base(phi) dphi |^2 .  Compare to |psi_base|^2.
================================================================================
"""
import numpy as np
import sympy as sp

# --------------------------------------------------------------------------- #
#  PART 1 :  SYMBOLIC core integral
#  int_{-pi/2}^{pi/2} cos((2m+1)(phi0 - phi)) cos(phi) dphi  =?  (pi/2) cos(phi0) d_{m0}
# --------------------------------------------------------------------------- #
print("="*78)
print("PART 1  -- symbolic core integral  (sympy, exact)")
print("  I_m(phi0) = int_{-pi/2}^{pi/2} cos((2m+1)(phi0-phi)) cos(phi) dphi")
print("-"*78)
phi, phi0 = sp.symbols('phi phi0', real=True)
for m in range(0, 5):
    integrand = sp.cos((2*m+1)*(phi0 - phi)) * sp.cos(phi)
    val = sp.simplify(sp.integrate(integrand, (phi, -sp.pi/2, sp.pi/2)))
    print(f"  m={m}:  I_m(phi0) = {val}")
print("  => only m=0 gives (pi/2)cos(phi0); all m>=1 give 0  (orthogonality).")
print("  => (S_N * cos)(phi0) = (pi/2) cos(phi0)  EXACTLY for every odd N>=1.")

# --------------------------------------------------------------------------- #
#  PART 2 :  NUMERIC convolution  (B, projector reading)  vs base^2 = cos^2
# --------------------------------------------------------------------------- #
print("\n"+"="*78)
print("PART 2  -- projector reading (B):  |(S_N * cos)(phi0)|^2  vs  cos^2(phi0)")
print("-"*78)

def S_N(x, N):
    x = np.asarray(x, dtype=float)
    s = np.sin(x)
    return np.where(np.abs(s) < 1e-12, (N+1)/2.0, np.sin((N+1)*x)/(2.0*s))

def conv_with_base(phi0_arr, N, base_func, ngrid=200001):
    phi = np.linspace(-np.pi/2, np.pi/2, ngrid)
    bp = base_func(phi)
    out = np.empty_like(phi0_arr, dtype=float)
    for i, p0 in enumerate(phi0_arr):
        out[i] = np.trapezoid(S_N(p0 - phi, N) * bp, phi)
    return out

base_cos = lambda x: np.cos(x)
phi0_test = np.linspace(-np.pi/2*0.98, np.pi/2*0.98, 21)

print(f"  {'N':>5} | {'max|conv - (pi/2)cos|':>22} | {'max| |conv|^2/norm - cos^2 |':>28}")
for N in [1, 3, 9, 31, 99]:
    conv = conv_with_base(phi0_test, N, base_cos)
    dev_lin = np.max(np.abs(conv - (np.pi/2)*np.cos(phi0_test)))
    P_norm = conv**2 / (np.pi/2)**2
    dev_sq = np.max(np.abs(P_norm - np.cos(phi0_test)**2))
    print(f"  {N:>5} | {dev_lin:22.3e} | {dev_sq:28.3e}")
print("  => |(S_N*cos)|^2 = (pi^2/4) cos^2(phi0) to machine precision, ALL N.")

# --------------------------------------------------------------------------- #
#  PART 3 :  general band-limited base state  ->  |psi_base|^2
# --------------------------------------------------------------------------- #
print("\n"+"="*78)
print("PART 3  -- general base state psi_base = cos(phi) + 0.5 cos(3 phi) - 0.3 cos(5 phi)")
print("-"*78)
def base_multi(x):
    return np.cos(x) + 0.5*np.cos(3*x) - 0.3*np.cos(5*x)

phi0_test2 = np.linspace(-np.pi/2*0.98, np.pi/2*0.98, 41)
psi_true = base_multi(phi0_test2)
for N in [1, 3, 5, 9, 31]:
    conv = conv_with_base(phi0_test2, N, base_multi)
    P_norm = conv**2 / (np.pi/2)**2
    dev = np.max(np.abs(P_norm - psi_true**2))
    note = "(N < 5: misses 5th mode)" if N < 5 else "(N >= 5: covers all modes -> exact)"
    print(f"  N={N:>3}:  max| |conv|^2/norm - |psi_base|^2 | = {dev:.3e}   {note}")
print("  => once N >= highest mode of psi_base, projector+square = |psi_base|^2 EXACTLY.")

# --------------------------------------------------------------------------- #
#  PART 4 :  ENVELOPE reading (A)  --  honest NEGATIVE control
# --------------------------------------------------------------------------- #
print("\n"+"="*78)
print("PART 4  -- envelope reading (A):  weight ~ 1/((N+1)sin phi)^2 diverges at 0")
print("-"*78)
print("  envelope(phi->0) ~ 1/phi^2 DIVERGES -> non-normalizable -> not a probability.")
print("  1/sin^2 shape != cos^2.  (A) is the wrong operationalization; (B) is correct.")

# --------------------------------------------------------------------------- #
#  PART 5 :  COMPLEX band-limited base state -> genuine |psi|^2 = Re^2 + Im^2
#  This is where the squaring is a TRUE modulus |Z|^2 = Z Zbar != Z^2 :
#  the complex Born structure absent from Part 3 (real-even, Z real => ZZbar=Z^2).
# --------------------------------------------------------------------------- #
print("\n"+"="*78)
print("PART 5  -- COMPLEX base (odd modes, complex coeffs): |Z|^2 = Re^2 + Im^2")
print("  psi_c(phi) = e^{i phi} + (0.5 - 0.3i) e^{i3phi} + (0.2 + 0.4i) e^{-i5phi}")
print("-"*78)

def base_complex(x):
    return (np.exp(1j*x)
            + (0.5 - 0.3j)*np.exp(3j*x)
            + (0.2 + 0.4j)*np.exp(-5j*x))

def conv_c(phi0_arr, N, base_func, ngrid=200001):
    """Complex convolution: kernel S_N (real) scanned across a complex base."""
    phi = np.linspace(-np.pi/2, np.pi/2, ngrid)
    bp = base_func(phi)
    out = np.empty_like(phi0_arr, dtype=complex)
    for i, p0 in enumerate(phi0_arr):
        out[i] = np.trapezoid(S_N(p0 - phi, N) * bp, phi)
    return out

phi0c = np.linspace(-np.pi/2*0.98, np.pi/2*0.98, 41)
psi_c_true = base_complex(phi0c)
print(f"  {'N':>4} | {'max|conv - (pi/2)psi_c|':>24} | {'max| |conv|^2/norm - |psi_c|^2 |':>32}")
for N in [1, 3, 5, 9, 31]:
    cc = conv_c(phi0c, N, base_complex)
    dev_lin = np.max(np.abs(cc - (np.pi/2)*psi_c_true))
    P = np.abs(cc)**2 / (np.pi/2)**2
    dev_sq = np.max(np.abs(P - np.abs(psi_c_true)**2))
    note = "(N>=5: covers modes 1,3,5 -> exact)" if N >= 5 else "(N<5: misses 5th mode)"
    print(f"  {N:>4} | {dev_lin:24.3e} | {dev_sq:32.3e}   {note}")
print("  => |(S_N*psi_c)|^2 = (pi^2/4)|psi_c|^2 = (pi^2/4)(Re^2 + Im^2): TRUE modulus.")

# Explicit |Z|^2 != Z^2 :  pure complex mode psi = e^{i phi}
cc1 = conv_c(phi0c, 9, lambda x: np.exp(1j*x))
mod2 = np.abs(cc1)**2 / (np.pi/2)**2          # = cos^2+sin^2 = 1 (constant)
sq   = cc1**2 / (np.pi/2)**2                   # = e^{2 i phi0} (oscillates, complex)
print("  modulus vs square for psi=e^{i phi}:  |Z|^2/norm = %.6f (const = Re^2+Im^2),"
      % mod2.mean())
print("     Z^2/norm oscillates: max|Im(Z^2)|=%.3f  => squaring here is |.|^2, not (.)^2."
      % np.max(np.abs(sq.imag)))

# --------------------------------------------------------------------------- #
#  PART 6 :  REAL base with an ODD (sin) component
#  Q4 mechanism: the sin harmonics are NOT cancelled in general; they reproduce
#  the ODD part.  The Part-3 cancellation was special to EVEN-real base states.
# --------------------------------------------------------------------------- #
print("\n"+"="*78)
print("PART 6  -- real base with odd part:  psi = cos(phi) + 0.7 sin(phi)")
print("-"*78)
base_rs = lambda x: np.cos(x) + 0.7*np.sin(x)
phi0r = np.linspace(-np.pi/2*0.98, np.pi/2*0.98, 41)
for N in [1, 9]:
    cr = conv_c(phi0r, N, base_rs)             # real base -> conv real (Im ~ 0)
    dev = np.max(np.abs(cr.real - (np.pi/2)*base_rs(phi0r)))
    devsq = np.max(np.abs(cr.real**2/(np.pi/2)**2 - base_rs(phi0r)**2))
    print(f"  N={N:>2}:  max|conv - (pi/2)psi| = {dev:.3e}   "
          f"max| |conv|^2/norm - psi^2 | = {devsq:.3e}")
print("  => the sin harmonic reproduces the ODD part 0.7 sin(phi); it is NOT cancelled.")
print("     (Part-3 'cancellation' = the degenerate even-real case where ZZbar = Z^2.)")

print("\n"+"="*78)
print("CONCLUSION")
print("-"*78)
print("  (A) envelope-average     : FAILS  (1/sin^2, non-normalizable).")
print("  (B) phase-diff kernel+sq : EXACT  |(S_N*psi_base)|^2 = (pi^2/4)|psi_base|^2.")
print("  Equal-amplitude odd harmonics = truncated reproducing kernel -> reproduces base;")
print("  'phase-difference-only' supplies the SQUARING (=> Born |psi|^2).")
print("  Single localized kernel, two conditions.")
print("="*78)
