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

from debroglie_align_lambda import align_lambda_search     # honest bounded search

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
    ap.add_argument("--dlam", type=float, default=0.01,
                    help="per-electron wavelength fluctuation half-width fraction: "
                         "lambda0' = lambda0 + eta, eta ~ cos^2 over +/- dlam*lambda0 "
                         "(independent of the position fluctuation).  0 disables it.")
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

    _ens_cache = {}

    def lam_used_of(V):
        """rough lambda0 and (used) lambda.  In align mode the alignment/NG is judged
        FAITHFULLY per source position (ensemble over y ~ cos^2, +/- lambda/2); the
        reported lambda' is the mean over OK electrons and 'ng_frac' the NG fraction.
        Cached per V so the sweep and fringe use the same ensemble."""
        lam0 = lam_of(V)
        if args.align:
            if V not in _ens_cache:
                _ens_cache[V] = align_ensemble(lam0)
            e = _ens_cache[V]
            ad = {"lam_star": (e["lam_mean"] if e["n_ok"] else None),
                  "m": e["m"], "rel": ((e["lam_mean"]/lam0 - 1.0) if e["n_ok"] else None),
                  "n_eval": e["n_eval"], "status": ("OK" if e["n_ok"] else "NG"),
                  "ng_frac": e["ng_frac"], "n_ok": e["n_ok"], "n": e["n"],
                  "lam_std": e["lam_std"], "lam_se": e["lam_se"],
                  "ok_y": e["ok_y"], "ok_lam": e["ok_lam"]}
            lamu = ad["lam_star"]
        else:
            ad = dict(align_lambda_search(L, W, lam0, N, y=0.0))   # reference at y=0
            ad["ng_frac"] = 0.0; ad["n_ok"] = args.M
            ad["lam_std"] = np.nan; ad["lam_se"] = np.nan          # (computed in sweep for raw)
            lamu = lam0
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

    def sample_wavelength(nn, center):
        """lambda0' = center + eta,  eta ~ cos^2(pi eta/(2 delta)) over [-delta, delta],
        delta = dlam*center.  INDEPENDENT MC draw (own rng calls) -> uncorrelated with
        the position sampling (no false correlation)."""
        delta = args.dlam * center
        if delta <= 0.0:
            return np.full(nn, center)
        ec  = rng.uniform(-delta, delta, size=3*nn)
        acc = rng.uniform(0.0, 1.0, size=3*nn) < np.cos(np.pi*ec/(2.0*delta))**2
        return center + ec[acc][:nn]

    def fringe_spacing(y, lam):
        dse = lam / W
        s = np.linspace(-args.nfr*dse, args.nfr*dse, args.ngrid)
        I = intensity(s, y, lam)
        up = (I[1:-1] >= I[:-2]) & (I[1:-1] > I[2:])
        sp = s[1:-1][up]
        return np.mean(np.diff(sp)) if len(sp) >= 2 else np.nan

    def accum_intensity(Xarr, lam):
        """Born push-forward, computed exactly (NOT a single centred source):
        each electron gets an independent position y ~ cos^2(+/-lam/2) AND initial
        wavelength lam0' ~ cos^2(+/-dlam*lam); RECOMPUTE the exact two-slit
        interference at (y, lam0') and accumulate -> the |psi|^2 build-up."""
        ys  = sample_source(args.M, Delta_of(lam))
        l0s = sample_wavelength(args.M, lam)
        n = min(len(ys), len(l0s))
        Itot = np.zeros_like(np.asarray(Xarr, dtype=float))
        for y, l0 in zip(ys[:n], l0s[:n]):
            Itot += intensity(Xarr / D, y, l0)       # source at y, its own wavelength lam0'
        return Itot

    def align_ensemble(lam0):
        """FAITHFUL per-electron alignment.  Each electron gets, INDEPENDENTLY,
        a source position y ~ cos^2 (+/- lambda0/2) AND an initial wavelength
        lambda0' ~ cos^2 (+/- dlam*lambda0); the alignment/NG is judged at THAT
        (y, lambda0') geometry -- no y=0 / fixed-lambda0 shortcut.  Returns the
        per-electron OK results (position y_j, aligned lambda'_j) and ensemble
        stats (mean lambda', NG fraction)."""
        ys  = sample_source(args.M, Delta_of(lam0))           # position, cos^2, +/- lam0/2
        l0s = sample_wavelength(args.M, lam0)                 # wavelength lam0', cos^2, +/- dlam
        n = min(len(ys), len(l0s))
        ok_y, ok_lam = [], []
        n_ng, nev, m_repr = 0, 0, None
        for y, l0 in zip(ys[:n], l0s[:n]):
            adj = align_lambda_search(L, W, l0, N, y=y)       # judged at (position y, wavelength lam0')
            nev += adj["n_eval"]
            if adj["status"] == "OK":
                ok_y.append(y); ok_lam.append(adj["lam_star"])
                if m_repr is None:
                    m_repr = adj["m"]
            else:
                n_ng += 1
        ok_lam = np.array(ok_lam); ok_y = np.array(ok_y)
        n_ok = len(ok_lam)
        lam_mean = (ok_lam.mean() if n_ok else np.nan)
        lam_std = (ok_lam.std(ddof=1) if n_ok > 1 else 0.0)   # sample std of lambda'
        lam_se = (lam_std / np.sqrt(n_ok) if n_ok > 0 else np.nan)  # std error of the mean
        return {"lam0": lam0, "ok_y": ok_y, "ok_lam": ok_lam, "n": n,
                "n_ok": n_ok, "n_ng": n_ng, "ng_frac": n_ng / max(n, 1),
                "lam_mean": lam_mean, "lam_std": lam_std, "lam_se": lam_se,
                "m": m_repr, "n_eval": nev}

    def accum_ensemble(Xarr, ens):
        """Accumulate the exact pattern over the ensemble's OK electrons, each at
        its OWN position y_j and its OWN aligned wavelength lambda'_j."""
        Itot = np.zeros_like(np.asarray(Xarr, dtype=float))
        for y, lam in zip(ens["ok_y"], ens["ok_lam"]):
            Itot += intensity(Xarr / D, y, lam)
        return Itot

    # ---- filename tag (parameters self-documented, incl. align flag) -------
    alg = "align" if args.align else "raw"
    fmtL = "L%gcm" % (L*100)
    fmtW = "W%gum" % (W*1e6)
    fmtD = "Dlam" if args.Delta == "lambda" else "D%gnm" % (float(args.Delta)*1e9)
    fmtVr = "V%g-%g" % (min(Vlist), max(Vlist))
    fmtDL = "dl%gpct" % (args.dlam * 100)                     # wavelength fluctuation
    tag = "%s_%s_N%d_%s_%s_%s_%s" % (fmtL, fmtW, N, fmtVr, fmtD, alg, fmtDL)
    outdir = os.path.dirname(os.path.abspath(__file__))

    def _cap(V, lam0, ad, align):
        """Per-panel caption.  align: per-position ensemble (mean lambda', NG frac)."""
        lstar = ad["lam_star"]
        if align and lstar is None:                          # all electrons NG
            return ((r"$V=%g$ V   (per-position align: ALL NG, %d/%d OK)" "\n"
                     r"$\lambda'=$NG,   $\lambda_0=%s$ m" "\n"
                     r"A-evals$=%d$;  pattern at $\lambda_0$ (does NOT interfere)")
                    % (V, ad["n_ok"], ad["n"], sci(lam0, 20), ad["n_eval"]))
        if align:
            spct = 100.0 * ad["lam_std"] / lstar; sepct = 100.0 * ad["lam_se"] / lstar
            return ((r"$V=%g$ V  (per-electron align: OK %d/%d, NG$_{\rm frac}=%.3f$)" "\n"
                     r"$\langle\lambda'\rangle=%s$ m  ($\sigma=%.4f$%%, SE$=%.4f$%%),   $\lambda_0=%s$ m" "\n"
                     r"$|\lambda_0-\langle\lambda'\rangle|=%s$ m  ($m^*=%d$, A-evals$=%d$)")
                    % (V, ad["n_ok"], ad["n"], ad["ng_frac"], sci(lstar, 20), spct, sepct,
                       sci(lam0, 20), sci(abs(lam0 - lstar), 6), ad["m"], ad["n_eval"]))
        if lstar is None:                                    # raw, y=0 reference NG
            return ((r"$V=%g$ V   (USED $\lambda_0$, rough)" "\n"
                     r"USED $\lambda_0=%s$ m      ref $\lambda'=$NG" "\n"
                     r"(reference at $y{=}0$; $\geq 2$ resonances)") % (V, sci(lam0, 20)))
        return ((r"$V=%g$ V   (USED $\lambda_0$, rough)" "\n"
                 r"USED $\lambda_0=%s$ m      ref $\lambda'=%s$ m" "\n"
                 r"($m^*=%d$)") % (V, sci(lam0, 20), sci(lstar, 20), ad["m"]))

    # ---- sweep : p vs 1/lambda + push-forward -----------------------------
    if args.out in ("sweep", "both"):
        Vv, pp, l0, lst, lm, mm, rel, nev, ngfr, nok, lstd, lse = ([] for _ in range(12))
        for V in Vlist:
            lam0, lamu, ad = lam_used_of(V); p = p_of(V)
            isng = args.align and ad["status"] == "NG"        # all electrons NG
            if args.align:
                lam_meas = np.nan if isng else ad["lam_star"]  # mean over OK; NaN if all NG
                sd, se = ad["lam_std"], ad["lam_se"]
            else:
                ys  = sample_source(args.M, Delta_of(lam0))
                l0s = sample_wavelength(args.M, lam0)
                nn = min(len(ys), len(l0s))
                lam_arr = W * np.array([fringe_spacing(y, l0) for y, l0 in zip(ys[:nn], l0s[:nn])])
                good = ~np.isnan(lam_arr); ng = int(good.sum())
                lam_meas = float(np.nanmean(lam_arr)) if ng else np.nan
                sd = float(np.nanstd(lam_arr, ddof=1)) if ng > 1 else 0.0
                se = sd / np.sqrt(ng) if ng > 0 else np.nan
            Vv.append(V); pp.append(p); l0.append(lam0); lst.append(ad["lam_star"])
            lm.append(lam_meas); mm.append(ad["m"]); rel.append(ad["rel"]); nev.append(ad["n_eval"])
            ngfr.append(ad["ng_frac"]); nok.append(ad["n_ok"]); lstd.append(sd); lse.append(se)
        Vv = np.array(Vv, float); pp = np.array(pp, float); l0 = np.array(l0, float)
        lm = np.array(lm, float); nev = np.array(nev, int)
        lstd = np.array(lstd, float); lse = np.array(lse, float)
        valid = ~np.isnan(lm)
        n_ng = int((~valid).sum())
        invlam = 1.0 / lm                                     # NaN where NG (align)
        A = np.polyfit(invlam[valid], pp[valid], 1)
        # push-forward of the source SHIFT at a representative (non-NG) V
        Vok = [float(v) for v, ok in zip(Vv, valid) if ok]
        Vp = 150.0 if 150.0 in Vok else (Vok[len(Vok)//2] if Vok else float(Vv[0]))
        _, lamu_p, _ = lam_used_of(Vp)
        if lamu_p is None:
            lamu_p = lam_of(Vp)
        ysh = sample_source(args.Nhist, Delta_of(lam_of(Vp)))
        r1 = np.sqrt(L**2+(ysh-y1)**2); r2 = np.sqrt(L**2+(ysh-y2)**2)
        u_deg = 360.0 * (-2*W*ysh/(r1+r2)) / lamu_p

        hV = pp * lm                       # per-V recovered h_V = p*<lambda'> (NaN where NG)
        yerr = pp * lse                    # SE of h_V  = p * SE(lambda')
        dev = hV - h
        span = np.abs(dev) + np.nan_to_num(yerr)             # reach of point + error bar
        maxdev = float(np.nanmax(span[valid])) if valid.any() else 0.0
        half = maxdev * 1.4 if maxdev > 0 else h * 1e-15   # zoom around CODATA h
        fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(17.6, 5.7))
        # --- left: p vs 1/lambda (de Broglie relation) ---
        a1.plot(invlam[valid], pp[valid], "o", color="#1f5fbf", ms=8, label="simulation (per V)")
        xx = np.linspace(invlam[valid].min(), invlam[valid].max(), 100)
        a1.plot(xx, A[0]*xx + A[1], "-", color="#d11f2d", lw=1.6,
                label=r"fit slope $=%.4e$" % A[0])
        a1.set_xlabel(r"$1/\lambda_{\rm meas}$  (m$^{-1}$)")
        a1.set_ylabel(r"$p=\sqrt{2m_eeV}$  (kg m/s)")
        a1.set_title(r"$p$ vs $1/\lambda$  ($N=%d$, %s): slope$/h=%.4f$   [NG excl.: %d]" % (N, alg, A[0]/h, n_ng))
        a1.legend(loc="upper left", fontsize=9); a1.grid(alpha=0.3)
        for i in range(len(Vv)):
            if valid[i]:
                a1.annotate("%dV" % int(Vv[i]), (invlam[i], pp[i]), fontsize=7,
                            textcoords="offset points", xytext=(4, -9))
        # --- middle: recovered h_V vs 1/lambda; CODATA h dashed at centre, zoomed ---
        a2.axhline(h, ls="--", color="k", lw=1.2, zorder=2,
                   label=r"theory $h=%.8e$ (CODATA)" % h)
        a2.errorbar(invlam[valid], hV[valid], yerr=yerr[valid], fmt="-o", color="#1f5fbf",
                    ms=6, lw=1.2, capsize=3, zorder=3,
                    label=r"experiment $h_V=p\langle\lambda'\rangle \pm p\cdot$SE")
        a2.set_ylim(h - half, h + half)
        a2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        a2.set_xlabel(r"$1/\lambda_{\rm meas}$  (m$^{-1}$)")
        a2.set_ylabel(r"$h_V=p\,\lambda_{\rm meas}$  (J s)")
        a2.set_title(r"recovered $h_V$ vs standard $h$   (max dev $=%.2e$ $=$ %.3g ppb)"
                     % (maxdev, maxdev/h*1e9))
        a2.legend(loc="upper right", fontsize=8); a2.grid(alpha=0.3)
        for i in range(len(Vv)):
            if valid[i]:
                a2.annotate("%dV" % int(Vv[i]), (invlam[i], hV[i]), fontsize=7,
                            textcoords="offset points", xytext=(4, 5))
        # --- right: push-forward of the source shift ---
        a3.hist(u_deg, bins=120, density=True, color="#9ecae1", edgecolor="#4a90c2", lw=0.4)
        a3.set_xlabel(r"central-fringe shift $u$ (deg) at $V=%g$ V" % Vp)
        a3.set_ylabel("density"); a3.set_title(r"push-forward ($\Delta=$%s)" % args.Delta)
        a3.grid(alpha=0.3)
        # BOTH lambdas: rough lambda0 vs aligned lambda' table
        hdr = ("  V[V]  lambda0[nm]  <lambda'>[nm]   <l'>/l0-1   sigma%%   SE%%     m*    "
               "NG_frac  n_OK/M  A-evals")
        tbl = [hdr]
        for i in range(len(Vv)):
            if lst[i] is None:
                tbl.append("  %5.0f  %.8f  %-13s %-10s %-7s %-7s %-9s %.3f  %d/%d  %d"
                           % (Vv[i], l0[i]*1e9, "NG(all)", "NG", "NG", "NG", "NG",
                              ngfr[i], nok[i], args.M, int(nev[i])))
            else:
                spct = 100.0*lstd[i]/lst[i]; sepct = 100.0*lse[i]/lst[i]
                tbl.append("  %5.0f  %.8f  %.8f  %+.2e  %.4f  %.4f  %-9d %.3f  %d/%d  %d"
                           % (Vv[i], l0[i]*1e9, lst[i]*1e9, rel[i], spct, sepct,
                              int(mm[i]), ngfr[i], nok[i], args.M, int(nev[i])))
        fig.text(0.5, 0.005, ("mode = %s   (per electron: position y~cos^2(+/-lam/2) AND wavelength "
                              "lam0'~cos^2(+/-%.1f%%), independent; alignment/NG judged FAITHFULLY at each "
                              "(y,lam0'); <lambda'>=mean over OK; NG_frac=NG/M; h only in lambda0 centre)\n"
                              % (alg, args.dlam*100))
                 + "\n".join(tbl), ha="center", va="bottom",
                 family="monospace", fontsize=7.4)
        fig.suptitle("de Broglie sweep  %s" % tag, fontsize=10.5, y=0.995)
        fig.subplots_adjust(left=0.055, right=0.99, top=0.88, bottom=0.34, wspace=0.28)
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(outdir, "debroglie_sweep_%s.%s" % (tag, ext)), dpi=170)
        plt.close(fig)
        print("sweep(N=%d,%s): slope=%.6e ratio=%.8f  (%d NG excluded)"
              % (N, alg, A[0], A[0]/h, n_ng))

    # ---- fringes : combined ('matome') + individual -----------------------
    if args.out in ("fringe", "both"):
        Xwin = args.nfringes * (lam_of(min(Vshow)) / W) * D      # window in metres (lambda0)
        nP = len(Vshow)
        figc, axs = plt.subplots(nP, 1, figsize=(10.8, 2.05*nP), sharex=True)
        axs = np.atleast_1d(axs)
        for ax, V in zip(axs, Vshow):
            lam0, lamu, ad = lam_used_of(V)
            X = np.linspace(-Xwin, Xwin, 6000)
            if args.align and ad["n_ok"]:
                Ii = accum_ensemble(X, ad)               # per-position, per-lambda' OK electrons
            else:
                Ii = accum_intensity(X, lam0)            # raw, or align all-NG -> lambda0
            Ii = Ii/Ii.max()
            ax.plot(X*1e3, Ii, color="#1f5fbf", lw=0.9)
            ax.set_ylim(0, 1.12); ax.set_ylabel(r"$I/I_{\max}$", fontsize=9); ax.grid(alpha=0.3)
            ax.text(0.995, 0.955, _cap(V, lam0, ad, args.align),
                    transform=ax.transAxes, ha="right", va="top", fontsize=6.7,
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
            if args.align and ad["n_ok"]:
                Ii = accum_ensemble(X, ad)               # per-position, per-lambda' OK electrons
            else:
                Ii = accum_intensity(X, lam0)            # raw, or align all-NG -> lambda0
            Ii = Ii/Ii.max()
            f, a = plt.subplots(figsize=(9.6, 3.9)); a.plot(X*1e3, Ii, color="#1f5fbf", lw=1.0)
            a.set_xlim(-Xwin*1e3, Xwin*1e3); a.set_ylim(0, 1.12)
            a.set_xlabel(r"screen position $X$ (mm)   ($D=%g$ m)" % D)
            a.set_ylabel(r"$I/I_{\max}$")
            a.set_title((r"Exact fringes  $N=%d$ (%s)" "\n" "%s")
                        % (N, alg, _cap(V, lam0, ad, args.align)), fontsize=8.5)
            a.grid(alpha=0.3); f.tight_layout()
            fn = "debroglie_fringe_%s_%s_N%d_V%03d_%s_%s" % (fmtL, fmtW, N, int(V), fmtD, alg)
            for ext in ("png", "svg"):
                f.savefig(os.path.join(outdir, "%s.%s" % (fn, ext)), dpi=170)
            plt.close(f)
        print("fringes(N=%d,%s): combined + %d individual saved  [%s]" % (N, alg, len(Vshow), tag))


if __name__ == "__main__":
    main()
