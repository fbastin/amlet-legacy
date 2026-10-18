#!/usr/bin/env python3
"""Recompute AMLET's log-likelihoods and standard errors independently.

Runs AMLET on models of the examples, then recomputes with numpy and scipy,
from the data and the uniform draws that AMLET saves (drawsname):

- australia/driver_mmnl (mixed logit, one normal coefficient, 1000 draws):
  the simulated log-likelihood, the robust standard errors (H^-1 B H^-1,
  with the exact Hessian) and, with robust = false, the classical ones
  (H^-1);
- simulated/driver_panel with 200 draws and a fixed sample size (panel
  mixed logit, four normal coefficients and a lognormal one): the same,
  and, with bias = 1, the simulated log-likelihood without the correction
  of the simulation bias, that the results also give;
- australia/driver_mnl with the first coefficient fixed at 4: the
  estimates, the log-likelihood and the robust standard errors.

    python3 verify_likelihood.py                 # uses ./local/bin/amlet
    AMLET=$HOME/amlet/bin/amlet python3 verify_likelihood.py
    RUN=/tmp/amlet-verify python3 verify_likelihood.py
"""
import os, re, shutil, subprocess, sys

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

ROOT = os.path.dirname(os.path.abspath(__file__))
AMLET = os.path.abspath(os.environ.get("AMLET", os.path.join(ROOT, "local", "bin", "amlet")))
RUN = os.path.abspath(os.environ.get("RUN", os.path.join(ROOT, "run", "verify")))

failed = 0


def report(model, what, amlet, python, tol):
    """Print one comparison; tol is the largest relative difference."""
    global failed
    amlet, python = np.atleast_1d(amlet).astype(float), np.atleast_1d(python).astype(float)
    if amlet.shape != python.shape:
        diff = np.inf
    else:
        diff = np.max(np.abs(amlet - python) / np.maximum(np.abs(python), 1e-300))
    ok = bool(diff <= tol)          # a missing value (NaN) fails
    failed += not ok
    print("%-16s %-26s %15.10g %15.10g %9.1e  %s"
          % (model, what, amlet[0], python[0], diff, "OK" if ok else "FAILED"))


def run_amlet(example, driver, name, edits):
    """Copy an example, edit its driver, run AMLET, return (dir, results)."""
    d = os.path.join(RUN, name)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(os.path.join(ROOT, "examples", example), d)
    path = os.path.join(d, driver)
    text = open(path, encoding="latin-1").read()
    for old, new in edits:
        if old not in text:
            sys.exit("ERROR: %r is not in %s" % (old, driver))
        text = text.replace(old, new, 1)
    text = re.sub(r'^output\s+"[^"]*"', 'output "results"', text, flags=re.M)
    open(path, "w", encoding="latin-1").write(text)
    p = subprocess.run([AMLET, driver], cwd=d, stdin=subprocess.DEVNULL,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    open(os.path.join(d, "amlet.out"), "wb").write(p.stdout)
    if p.returncode:
        sys.exit("ERROR: amlet %s failed (exit code %d), see %s" % (driver, p.returncode, d))
    return d, open(os.path.join(d, "results"), encoding="latin-1").read()


def estimates(res):
    """Estimates and standard errors of the estimated parameters."""
    rows = re.findall(r"^[ \t]+\d+\.[ \t]+(\S+)(?:[ \t]+(\S+))?", res, re.M)
    theta = np.array([float(e) for e, s in rows if s])
    se = np.array([float(s) for e, s in rows if s])
    return theta, se


def total_ll(res):
    return float(re.search(r"^Optimal\s+(\S+)", res, re.M).group(1))


def sandwich(scores, hessian):
    """Robust and classical standard errors."""
    hinv = np.linalg.inv(-hessian)
    return (np.sqrt(np.diag(hinv @ scores.T @ scores @ hinv)),
            np.sqrt(np.diag(hinv)))


def fd_hessian(grad, theta):
    """Hessian by central differences of an analytic gradient."""
    n = len(theta)
    h = np.empty((n, n))
    for k in range(n):
        e = np.zeros(n)
        e[k] = 1e-5 * max(1.0, abs(theta[k]))
        h[:, k] = (grad(theta + e) - grad(theta - e)) / (2 * e[k])
    return (h + h.T) / 2


# Mixed logit, Australia: 210 individuals, 4 alternatives, 7 rows per
# record (6 variables, then the choice), normal 5th coefficient.

def australia_data(d):
    raw = np.loadtxt(os.path.join(d, "model_australia")).reshape(210, 7, 4)
    return raw[:, :6, :].transpose(0, 2, 1), raw[:, 6, :].argmax(axis=1)


def mmnl_scores(theta, X, y, Z):
    """Log of the simulated probabilities and their gradients."""
    n, R = Z.shape
    B = np.empty((n, R, 6))
    B[:, :, :4] = theta[:4]
    B[:, :, 4] = theta[4] + theta[5] * Z
    B[:, :, 5] = theta[6]
    V = np.einsum("iak,irk->ira", X, B)
    P = np.exp(V - V.max(axis=2, keepdims=True))
    P /= P.sum(axis=2, keepdims=True)
    p = P[np.arange(n), :, y]                                   # (n, R)
    sp = p.mean(axis=1)
    dv = X[np.arange(n), y][:, None, :] - np.einsum("ira,iak->irk", P, X)
    w = (p / (R * sp[:, None]))[:, :, None] * dv                # d log sp / d beta, per draw
    s = np.empty((n, 7))
    s[:, :4] = w[:, :, :4].sum(axis=1)
    s[:, 4] = w[:, :, 4].sum(axis=1)
    s[:, 5] = (w[:, :, 4] * Z).sum(axis=1)
    s[:, 6] = w[:, :, 5].sum(axis=1)
    return np.log(sp), s


def verify_australia_mmnl():
    common = [('drawsname = "draws";', 'drawsname = "draws";')]   # saves the draws
    d, res = run_amlet("australia", "driver_mmnl", "mmnl_robust", common)
    X, y = australia_data(d)
    Z = norm.ppf(np.loadtxt(os.path.join(d, "draws")).reshape(210, -1))
    theta, se = estimates(res)
    lsp, s = mmnl_scores(theta, X, y, Z)
    hess = fd_hessian(lambda t: mmnl_scores(t, X, y, Z)[1].sum(axis=0), theta)
    robust, classical = sandwich(s, hess)
    report("australia mmnl", "log-likelihood", total_ll(res), lsp.sum(), 1e-8)
    report("australia mmnl", "robust standard errors", se, robust, 1e-3)
    d2, res2 = run_amlet("australia", "driver_mmnl", "mmnl_classical",
                         common + [("robust = true;", "robust = false;")])
    report("australia mmnl", "classical standard errors", estimates(res2)[1], classical, 1e-3)


# Panel mixed logit: simulated data written by AMLET (Train99), 2500
# individuals, 5000 observations, 5 alternatives, 5 random coefficients
# (normal, the last one lognormal: exp(m + s z)).

NOBS, NALT, NVAR, NP = 5000, 5, 5, 2500


def panel_data(d, R):
    x = np.array(open(os.path.join(d, "panel_data_x")).read().split(), float)
    x = x.reshape(NOBS, NVAR + 1, NALT)                        # 6th row: availability
    if not (x[:, NVAR, :] == 1).all():
        sys.exit("ERROR: the panel data have unavailable alternatives")
    y = np.loadtxt(os.path.join(d, "panel_data_y"), dtype=int) - 1
    times = np.loadtxt(os.path.join(d, "panel_data_times"), dtype=int)
    ind = np.repeat(np.arange(NP), times)
    u = np.array(open(os.path.join(d, "pdraws")).read().split(), float)
    return x[:, :NVAR, :].transpose(0, 2, 1), y, ind, norm.ppf(u.reshape(NP, R, NVAR))


def panel_scores(theta, X, y, ind, Z):
    R = Z.shape[1]
    m, s = theta[0::2], theta[1::2]
    B = m + s * Z                                              # (NP, R, NVAR)
    B[:, :, 4] = np.exp(B[:, :, 4])
    Bo = B[ind]                                                # (NOBS, R, NVAR)
    V = np.einsum("oak,ork->ora", X, Bo)
    V -= V.max(axis=2, keepdims=True)
    E = np.exp(V)
    P = E / E.sum(axis=2, keepdims=True)
    lp = np.log(P[np.arange(NOBS), :, y])                      # (NOBS, R)
    lpn = np.zeros((NP, R))
    np.add.at(lpn, ind, lp)                                    # product over the panel
    mx = lpn.max(axis=1, keepdims=True)
    L = np.exp(lpn - mx)
    lsp = mx[:, 0] + np.log(L.mean(axis=1))
    w = L / L.sum(axis=1, keepdims=True)                       # (NP, R)
    dv = X[np.arange(NOBS), y][:, None, :] - np.einsum("ora,oak->ork", P, X)
    D = np.zeros((NP, R, NVAR))
    np.add.at(D, ind, dv)                                      # d log prod P / d beta
    D[:, :, 4] *= B[:, :, 4]                                   # lognormal: d beta = beta d(m + s z)
    sc = np.empty((NP, 2 * NVAR))
    sc[:, 0::2] = (w[:, :, None] * D).sum(axis=1)
    sc[:, 1::2] = (w[:, :, None] * D * Z).sum(axis=1)
    return lsp, sc


def verify_panel():
    R = 200
    edits = [("nrep = 1000;", 'nrep = %d; drawsname = "pdraws";' % R),
             ("METHOD = MCBTRDA;", "METHOD = MCBTR;"),
             ("statistics = true;", "statistics = false;")]
    d, res = run_amlet("simulated", "driver_panel", "panel",
                       edits + [("bias=1;", "bias=0;")])
    X, y, ind, Z = panel_data(d, R)
    theta, se = estimates(res)
    lsp, s = panel_scores(theta, X, y, ind, Z)
    hess = fd_hessian(lambda t: panel_scores(t, X, y, ind, Z)[1].sum(axis=0), theta)
    robust, classical = sandwich(s, hess)
    report("panel mmnl", "log-likelihood", total_ll(res), lsp.sum(), 1e-8)
    report("panel mmnl", "robust standard errors", se, robust, 1e-3)
    # bias = 1: the results also give the log-likelihood without the correction
    d, res = run_amlet("simulated", "driver_panel", "panel_bias", edits)
    X, y, ind, Z = panel_data(d, R)
    m = re.search(r"the simulated log-likelihood is (\S+)", res)
    report("panel mmnl", "uncorrected LL (bias = 1)", float(m.group(1)) if m else np.nan,
           panel_scores(estimates(res)[0], X, y, ind, Z)[0].sum(), 1e-8)


# Multinomial logit, Australia, with the first coefficient fixed at 4.

def verify_fixed_mnl():
    d, res = run_amlet("australia", "driver_mnl", "mnl_fixed",
                       [("btype = { CONSTANT,", "btype = { FIXED,"),
                        ("b = { 0.0,", "b = { 4.0,")])
    X, y = australia_data(d)
    Xc = X[np.arange(210), y]

    def parts(free):
        V = X @ np.concatenate([[4.0], free])
        P = np.exp(V - V.max(axis=1, keepdims=True))
        P /= P.sum(axis=1, keepdims=True)
        s = (Xc - np.einsum("ia,iak->ik", P, X))[:, 1:]
        return np.log(P[np.arange(210), y]).sum(), s, P

    r = minimize(lambda f: (-parts(f)[0], -parts(f)[1].sum(axis=0)), np.zeros(5),
                 jac=True, method="BFGS", options={"gtol": 1e-10})
    ll, s, P = parts(r.x)
    Xf = X[:, :, 1:] - np.einsum("ia,iak->ik", P, X[:, :, 1:])[:, None, :]
    hess = -np.einsum("ia,iak,ial->kl", P, Xf, Xf)
    theta, se = estimates(res)
    report("mnl, 1 fixed", "estimates", theta, r.x, 1e-5)
    report("mnl, 1 fixed", "log-likelihood", total_ll(res), ll, 1e-8)
    report("mnl, 1 fixed", "robust standard errors", se, sandwich(s, hess)[0], 1e-3)


def main():
    if not os.access(AMLET, os.X_OK):
        sys.exit("ERROR: AMLET executable not found at %s (run build.sh, or set AMLET=...)" % AMLET)
    os.makedirs(RUN, exist_ok=True)
    print("%-16s %-26s %15s %15s %9s  %s" % ("model", "quantity", "AMLET", "Python",
                                            "rel. diff", "status"))
    print("(standard errors: first parameter shown, largest difference over all)")
    verify_australia_mmnl()
    verify_panel()
    verify_fixed_mnl()
    print("\n%s. Outputs in %s" % ("All agree" if not failed else "%d FAILED" % failed, RUN))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
