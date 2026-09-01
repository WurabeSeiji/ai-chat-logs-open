
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

STATE_NPZ = r"/mnt/data/denominator_controls_64bit_with124_20260901/hm_N5_den_124_states_500.npz"
PARENT_NPZ = r"/mnt/data/self_crosscheck/exact_reproduction/data/hm_N5/parent_v.npz"
OUT = r"/mnt/data/ab_plot_audit_20260901/N5_L124_500_B_new_program.png"

Z = np.load(STATE_NPZ)["Z"]
assert Z.shape == (501, 10)
assert Z.dtype == np.complex128
p, q = parent_plane(PARENT_NPZ)
vals = hperp_fraction(Z, p, q)
draw_curve(vals, "B: ChatGPT denominator-control, N=5, L=124, 500 steps, complex128", OUT)
print("step0", repr(float(vals[0])))
print("step1", repr(float(vals[1])))
print("step500", repr(float(vals[-1])))
print("saved", OUT)
