#!/usr/bin/env python3
"""
de Broglie p-lambda sweep + EXACT two-slit fringes of a localized odd-harmonic
wave.  FULLY PARAMETERIZED (L, W, N, Delta, V, D, --align ... via argparse).

Physics (unchanged from Paper 1 / fig_oddharm_interference):
  source = equal-amplitude odd harmonics n = 1,3,...,N,  lambda_n = lambda/n
  exact two-slit intensity  I(s;y) = | sum_n ( e^{iPhi1^n} + e^{iPhi2^n} ) |^2 ,
    Phi_k^n = 2 pi n (r_k - y_slit,k s)/lambda ,  r_k = sqrt(L^2 + (y - y_slit,k)^2)
  via the algebraically-IDENTICAL, precision-safe factored form
    e^{iPhi1^n}+e^{iPhi2^n} = e^{i 2 pi n r_avg/lambda} * 2 cos( pi n [(r1-r2) - W s]/lambda ),
    r1 - r2 = -2 W y/(r1+r2)  (difference-of-squares; no cancellation),  r_avg=(r1+r2)/2.

de Broglie:  p = sqrt(2 m_e e V),  lambda0 = h/p  (approximate/rough wavelength).

  --align OFF (default) : the pattern is evaluated at the rough lambda0.  For N=1 the
      global alignment phase drops out (machine precision).  For N>1 the inter-harmonic
      alignment phases (r_avg/lambda ~ 1e9) do NOT line up at an arbitrary V -> the
      pattern SCRAMBLES (reported as-is; no interference).
  --align ON : lambda0 is snapped to the DETERMINISTIC aligned wavelength
      lambda* = 2 r_k/round(2 r_k/lambda0)  (module debroglie_align_lambda), so that
      r_k/lambda* = m*/2 and every odd harmonic shares a common phase -> a SHARP
      localized fringe forms even for N>1.  The adjustment |lambda*-lambda0|/lambda0
      ~1e-9 is tiny yet essential: "a hair off lambda and it will not interfere".

Both figures annotate BOTH the rough lambda0 and the aligned lambda'(=lambda*).
Filenames carry every parameter incl. the align flag (..._align / ..._raw).

Usage:
  python3 debroglie_plambda_sweep.py                 # N=1, rough lambda0
  python3 debroglie_plambda_sweep.py --N 9 --align   # N=9, aligned lambda*
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from debroglie_align_lambda import align_lambda            # deterministic solver

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

h  = 6.62607015e-34
me = 9.1093837015e-31
e  = 1.602176634e-19


def parse_list(s):
    return [float(x) for x in str(s).split(",") if x.strip() != ""]


def sci(x, digits=20):
    """High-precision scientific string as mathtext: 1.7344...e-10 -> '1.7344...\\times10^{-10}'.
    Shows `digits` mantissa decimals (float64 carries ~16 significant figures; the
    rough/aligned difference appears near the 10th-11th figure, well within that)."""
    m, ex = ("%.*e" % (digits, x)).split("e")
    return r"%s\times10^{%d}" % (m, int(ex))


def main():
    ap = argparse.ArgumentParser(description="de Broglie p-lambda sweep + exact fringes")
    # geometry
    ap.add_argument("--L", type=float, default=0.05, help="source-slit distance [m] (5 cm)")
    ap.add_argument("--W", type=float, default=5e-6, help="slit separation [m] (5 um)")
    # wave / model
    ap.add_argument("--N", type=int, default=1, help="highest odd harmonic order (1,3,...,N)")
    ap.add_argument("--Delta", type=str, default="lambda",
                    help="source-fluctuation FULL width: 'lambda' -> Delta=lambda(V), or a number [m]")
    # voltage sweep
    ap.add_argument("--V", type=str, default="50,100,150,200,300,400,500",
                    help="accelerating-voltage list [V], comma-separated")
    # alignment
    ap.add_argument("--align", action="store_true",
                    help="snap rough lambda0 to the deterministic aligned lambda* (2 r_k/m*)")
    # display (fringe figures)
    ap.add_argument("--D", type=float, default=1.0, help="display screen distance [m]")
    ap.add_argument("--nfringes", type=int, default=5, help="half-window in fringes (widest-fringe V)")
    ap.add_argument("--Vshow", type=str, default="50,100,150,300,500",
                    help="voltages to draw fringe figures for [V], comma-separated")
    # sampling
    ap.add_argument("--M", type=int, default=200, help="source samples per V")
    ap.add_argument("--Nhist", type=int, default=200000, help="source samples for push-forward histogram")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed")
    # numerics
    ap.add_argument("--ngrid", type=int, default=20001, help="s-grid points")
    ap.add_argument("--nfr", type=int, default=8, help="fringes used for spacing measurement")
    # output selection
    ap.add_argument("--out", type=str, default="both", choices=["sweep", "fringe", "both"])
    args = ap.parse_args()

    L, W, N, D = args.L, args.W, args.N, args.D
    y1, y2 = +W / 2.0, -W / 2.0
    n_list = np.arange(1, N + 1, 2)
    Vlist = parse_list(args.V)
    Vshow = parse_list(args.Vshow)
    rng = np.random.default_rng(args.seed)

    lam_of = lambda V: h / np.sqrt(2 * me * e * V)
    p_of   = lambda V: np.sqrt(2 * me * e * V)
    Delta_of = lambda lam: lam if args.Delta == "lambda" else float(args.Delta)

    def lam_used_of(V):
        """rough lambda0 and (used) lambda: lambda* if --align else lambda0."""
        lam0 = lam_of(V)
        ad = align_lambda(L, W, lam0)
        lamu = ad["lam_star"] if args.align else lam0
        return lam0, lamu, ad

    def intensity(s, y, lam):
        r1 = np.sqrt(L**2 + (y - y1)**2)
        r2 = np.sqrt(L**2 + (y - y2)**2)
        dr = -2.0 * W * y / (r1 + r2)                 # r1 - r2  (no cancellation)
        ravg = 0.5 * (r1 + r2)
        psi = np.zeros(np.shape(s), dtype=complex)
        for n in n_list:
            psi += np.exp(1j * 2*np.pi * n * ravg / lam) \
                   * 2.0 * np.cos(np.pi * n * (dr - W * s) / lam)
        return (psi * np.conj(psi)).real

    def sample_source(nn, Delta):
        yc  = rng.uniform(-Delta/2.0, Delta/2.0, size=3*nn)
        acc = rng.uniform(0.0, 1.0, size=3*nn) < np.cos(np.pi*yc/Delta)**2
        return yc[acc][:nn]

    def fringe_spacing(y, lam):
        dse = lam / W
        s = np.linspace(-args.nfr*dse, args.nfr*dse, args.ngrid)
        I = intensity(s, y, lam)
        up = (I[1:-1] >= I[:-2]) & (I[1:-1] > I[2:])
        sp = s[1:-1][up]
        return np.mean(np.diff(sp)) if len(sp) >= 2 else np.nan

    # ---- filename tag (parameters self-documented, incl. align flag) -------
    alg = "align" if args.align else "raw"
    fmtL = "L%gcm" % (L*100)
    fmtW = "W%gum" % (W*1e6)
    fmtD = "Dlam" if args.Delta == "lambda" else "D%gnm" % (float(args.Delta)*1e9)
    fmtVr = "V%g-%g" % (min(Vlist), max(Vlist))
    tag = "%s_%s_N%d_%s_%s_%s" % (fmtL, fmtW, N, fmtVr, fmtD, alg)
    outdir = os.path.dirname(os.path.abspath(__file__))

    # ---- sweep : p vs 1/lambda + push-forward -----------------------------
    if args.out in ("sweep", "both"):
        rows = []
        for V in Vlist:
            lam0, lamu, ad = lam_used_of(V); p = p_of(V)
            if args.align:
                lam_meas = ad["lam_star"]                       # deterministic aligned lambda
            else:
                ys = sample_source(args.M, Delta_of(lam0))
                ds = np.array([fringe_spacing(y, lam0) for y in ys])
                lam_meas = W * np.nanmean(ds)
            rows.append((V, p, lam0, ad["lam_star"], lam_meas, ad["m"], ad["rel"]))
        rows = np.array(rows)
        Vv, pp, l0, lst, lm, mm, rel = rows.T
        invlam = 1.0/lm
        A = np.polyfit(invlam, pp, 1)
        # push-forward of the source SHIFT at a representative V
        Vp = 150.0 if 150.0 in Vlist else Vlist[len(Vlist)//2]
        _, lamu_p, _ = lam_used_of(Vp)
        ysh = sample_source(args.Nhist, Delta_of(lam_of(Vp)))
        r1 = np.sqrt(L**2+(ysh-y1)**2); r2 = np.sqrt(L**2+(ysh-y2)**2)
        u_deg = 360.0 * (-2*W*ysh/(r1+r2)) / lamu_p

        hV = pp * lm                       # per-V recovered Planck const  h_V = p*lam_meas
        dev = hV - h
        maxdev = float(np.max(np.abs(dev)))
        half = maxdev * 1.4 if maxdev > 0 else h * 1e-15   # zoom around CODATA h
        fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(17.6, 5.7))
        # --- left: p vs 1/lambda (de Broglie relation) ---
        a1.plot(invlam, pp, "o", color="#1f5fbf", ms=8, label="simulation (per V)")
        xx = np.linspace(invlam.min(), invlam.max(), 100)
        a1.plot(xx, A[0]*xx + A[1], "-", color="#d11f2d", lw=1.6,
                label=r"fit slope $=%.4e$" % A[0])
        a1.set_xlabel(r"$1/\lambda_{\rm meas}$  (m$^{-1}$)")
        a1.set_ylabel(r"$p=\sqrt{2m_eeV}$  (kg m/s)")
        a1.set_title(r"$p$ vs $1/\lambda$  ($N=%d$, %s): slope$/h=%.4f$" % (N, alg, A[0]/h))
        a1.legend(loc="upper left", fontsize=9); a1.grid(alpha=0.3)
        for i, V in enumerate(Vv):
            a1.annotate("%dV" % V, (invlam[i], pp[i]), fontsize=7,
                        textcoords="offset points", xytext=(4, -9))
        # --- middle: recovered h_V vs 1/lambda; CODATA h dashed at centre, zoomed ---
        a2.axhline(h, ls="--", color="k", lw=1.2, zorder=2,
                   label=r"theory $h=%.8e$ (CODATA)" % h)
        a2.plot(invlam, hV, "-o", color="#1f5fbf", ms=7, lw=1.3, zorder=3,
                label=r"experiment $h_V=p\,\lambda_{\rm meas}$")
        a2.set_ylim(h - half, h + half)
        a2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        a2.set_xlabel(r"$1/\lambda_{\rm meas}$  (m$^{-1}$)")
        a2.set_ylabel(r"$h_V=p\,\lambda_{\rm meas}$  (J s)")
        a2.set_title(r"recovered $h_V$ vs standard $h$   (max dev $=%.2e$ $=$ %.3g ppb)"
                     % (maxdev, maxdev/h*1e9))
        a2.legend(loc="upper right", fontsize=8); a2.grid(alpha=0.3)
        for i, V in enumerate(Vv):
            a2.annotate("%dV" % V, (invlam[i], hV[i]), fontsize=7,
                        textcoords="offset points", xytext=(4, 5))
        # --- right: push-forward of the source shift ---
        a3.hist(u_deg, bins=120, density=True, color="#9ecae1", edgecolor="#4a90c2", lw=0.4)
        a3.set_xlabel(r"central-fringe shift $u$ (deg) at $V=%g$ V" % Vp)
        a3.set_ylabel("density"); a3.set_title(r"push-forward ($\Delta=$%s)" % args.Delta)
        a3.grid(alpha=0.3)
        # BOTH lambdas: rough lambda0 vs aligned lambda' table
        hdr = "  V[V]   lambda0 (rough)[nm]   lambda' (aligned)[nm]   d_lambda=l'-l0 [m]    l'/l0-1     m*"
        tbl = [hdr]
        for i in range(len(Vv)):
            tbl.append("  %5.0f      %.8f           %.8f         %+.3e      %+.2e   %d"
                       % (Vv[i], l0[i]*1e9, lst[i]*1e9, (lst[i]-l0[i]), rel[i], int(mm[i])))
        fig.text(0.5, 0.005, ("mode = %s   (lambda used = %s)\n" %
                              (alg, "lambda' aligned" if args.align else "lambda0 rough"))
                 + "\n".join(tbl), ha="center", va="bottom",
                 family="monospace", fontsize=7.6)
        fig.suptitle("de Broglie sweep  %s" % tag, fontsize=10.5, y=0.995)
        fig.subplots_adjust(left=0.055, right=0.99, top=0.88, bottom=0.34, wspace=0.28)
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(outdir, "debroglie_sweep_%s.%s" % (tag, ext)), dpi=170)
        plt.close(fig)
        print("sweep(N=%d,%s): slope=%.6e (h=%.6e, ratio=%.8f)"
              % (N, alg, A[0], h, A[0]/h))

    # ---- fringes : combined ('matome') + individual -----------------------
    if args.out in ("fringe", "both"):
        _, lamu_min, _ = lam_used_of(min(Vshow))
        Xwin = args.nfringes * (lamu_min / W) * D                # window in metres
        nP = len(Vshow)
        figc, axs = plt.subplots(nP, 1, figsize=(10.8, 2.05*nP), sharex=True)
        axs = np.atleast_1d(axs)
        for ax, V in zip(axs, Vshow):
            lam0, lamu, ad = lam_used_of(V)
            X = np.linspace(-Xwin, Xwin, 6000)
            Ii = intensity(X/D, 0.0, lamu); Ii = Ii/Ii.max()
            ax.plot(X*1e3, Ii, color="#1f5fbf", lw=0.9)
            ax.set_ylim(0, 1.12); ax.set_ylabel(r"$I/I_{\max}$", fontsize=9); ax.grid(alpha=0.3)
            if args.align:
                ulab, uval, olab, oval = r"\lambda'", ad["lam_star"], r"\lambda_0", lam0
            else:
                ulab, uval, olab, oval = r"\lambda_0", lam0, r"\lambda'", ad["lam_star"]
            dlam_abs = abs(lam0 - ad["lam_star"])
            ax.text(0.995, 0.955,
                    (r"$V=%g$ V   (USED $\lambda=%s$, %s)" "\n"
                     r"USED $%s=%s$ m      $%s=%s$ m" "\n"
                     r"$|\Delta\lambda|=|\lambda_0-\lambda'|=%s$ m   ($m^*=%d$)")
                    % (V, ulab, ("aligned" if args.align else "rough"),
                       ulab, sci(uval, 20), olab, sci(oval, 20),
                       sci(dlam_abs, 6), ad["m"]),
                    transform=ax.transAxes, ha="right", va="top", fontsize=6.9,
                    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.92))
        axs[-1].set_xlabel(r"screen position $X$ (mm)   ($D=%g$ m)" % D)
        figc.suptitle(r"Exact two-slit fringes ($N=%d$, %s: $\lambda$ used $=%s$), combined   %s"
                      % (N, alg, r"\lambda'" if args.align else r"\lambda_0", tag), fontsize=10.3)
        figc.tight_layout(rect=[0, 0, 1, 0.975])
        for ext in ("png", "svg"):
            figc.savefig(os.path.join(outdir, "debroglie_fringes_combined_%s.%s" % (tag, ext)), dpi=170)
        plt.close(figc)
        for V in Vshow:
            lam0, lamu, ad = lam_used_of(V)
            X = np.linspace(-Xwin, Xwin, 6000)
            Ii = intensity(X/D, 0.0, lamu); Ii = Ii/Ii.max()
            f, a = plt.subplots(figsize=(9.6, 3.6)); a.plot(X*1e3, Ii, color="#1f5fbf", lw=1.0)
            a.set_xlim(-Xwin*1e3, Xwin*1e3); a.set_ylim(0, 1.12)
            a.set_xlabel(r"screen position $X$ (mm)   ($D=%g$ m)" % D)
            a.set_ylabel(r"$I/I_{\max}$")
            if args.align:
                ulab, uval, olab, oval = r"\lambda'", ad["lam_star"], r"\lambda_0", lam0
            else:
                ulab, uval, olab, oval = r"\lambda_0", lam0, r"\lambda'", ad["lam_star"]
            dlam_abs = abs(lam0 - ad["lam_star"])
            a.set_title((r"Exact fringes  $N=%d$ (%s),  $V=%g$ V   —   USED $\lambda=%s$ (%s)" "\n"
                         r"USED $%s=%s$ m      $%s=%s$ m" "\n"
                         r"$|\Delta\lambda|=|\lambda_0-\lambda'|=%s$ m   ($m^*=%d$)")
                        % (N, alg, V, ulab, ("aligned" if args.align else "rough"),
                           ulab, sci(uval, 20), olab, sci(oval, 20),
                           sci(dlam_abs, 6), ad["m"]),
                        fontsize=8.5)
            a.grid(alpha=0.3); f.tight_layout()
            fn = "debroglie_fringe_%s_%s_N%d_V%03d_%s_%s" % (fmtL, fmtW, N, int(V), fmtD, alg)
            for ext in ("png", "svg"):
                f.savefig(os.path.join(outdir, "%s.%s" % (fn, ext)), dpi=170)
            plt.close(f)
        print("fringes(N=%d,%s): combined + %d individual saved  [%s]" % (N, alg, len(Vshow), tag))


if __name__ == "__main__":
    main()
