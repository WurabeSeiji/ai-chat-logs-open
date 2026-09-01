
import numpy as np
import matplotlib.pyplot as plt

def parent_plane(parent_npz):
    pz = np.load(parent_npz)
    v = np.asarray(pz["v"], dtype=np.complex128)
    p = v.real.astype(np.float64, copy=True)
    p /= np.linalg.norm(p)
    q = v.imag.astype(np.float64, copy=True)
    q -= np.dot(q, p) * p
    q /= np.linalg.norm(q)
    return p, q

def hperp_fraction(states, p, q):
    vals = np.empty(states.shape[0], dtype=np.float64)
    for i, z in enumerate(states):
        z = np.asarray(z, dtype=np.complex128)
        h = np.vdot(z, z).real
        zp = z - p*np.dot(p, z) - q*np.dot(q, z)
        hp = np.vdot(zp, zp).real
        vals[i] = hp/h
    return vals

def draw_curve(vals, title, out_png):
    steps = np.arange(vals.size)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.semilogy(steps, vals, linewidth=1.2)
    ax.set_xlim(0, 500)
    ax.set_ylim(1e-70, 1e-26)
    ax.set_xlabel("step")
    ax.set_ylabel("Hperp/H")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

A_STATE = r"/mnt/data/self_crosscheck/exact_reproduction/data/hm_N5/states_treatment.npz"
B_STATE = r"/mnt/data/denominator_controls_64bit_with124_20260901/hm_N5_den_124_states_500.npz"
PARENT_NPZ = r"/mnt/data/self_crosscheck/exact_reproduction/data/hm_N5/parent_v.npz"
OUT = r"/mnt/data/ab_plot_audit_20260901/N5_L124_500_A_vs_B_overlay.png"

ZA = np.load(A_STATE)["Z"]
ZB = np.load(B_STATE)["Z"]
assert np.array_equal(ZA, ZB), "A/B state arrays are not bitwise identical"
p, q = parent_plane(PARENT_NPZ)
a = hperp_fraction(ZA, p, q)
b = hperp_fraction(ZB, p, q)

steps = np.arange(501)
fig, ax = plt.subplots(figsize=(9, 6))
ax.semilogy(steps, a, linewidth=1.6, label="A canonical-copy")
ax.semilogy(steps, b, linestyle="--", linewidth=1.0, label="B denominator-control")
ax.set_xlim(0, 500)
ax.set_ylim(1e-70, 1e-26)
ax.set_xlabel("step")
ax.set_ylabel("Hperp/H")
ax.set_title("A vs B: N=5, L=124, 500 steps")
ax.grid(True, which="both", alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig(OUT, dpi=200)
plt.close(fig)

print("np.array_equal(Z_A,Z_B) =", np.array_equal(ZA,ZB))
print("max_abs_state_difference =", float(np.max(np.abs(ZA-ZB))))
print("max_abs_Hperp_fraction_difference =", float(np.max(np.abs(a-b))))
print("saved", OUT)
