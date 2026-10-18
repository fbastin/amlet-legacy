# Changelog

Changes in this repository, the 2026 rebuild of AMLET 0.13.1 (with ORATIO
0.7.5 and OPHELIA 0.7.4), newest first. The notes of each version are those
of its [GitHub release](https://github.com/fbastin/amlet-legacy/releases).
The releases are archived on Zenodo from v0.13.1-4:
[doi:10.5281/zenodo.23454086](https://doi.org/10.5281/zenodo.23454086)
designates all of them.

## v0.13.1-5 — 2026-10-18

Archived on Zenodo: [doi:10.5281/zenodo.23458742](https://doi.org/10.5281/zenodo.23458742).

AMLET 0.13.1 (2026 rebuild), sixth release. It fixes `draws = old;` in batch jobs, computes the log-likelihood at constants, gives the p-value of the likelihood ratio test, and adds a script that recomputes the log-likelihoods and standard errors with numpy and scipy. **The estimates, log-likelihoods and standard errors of the 19 examples do not change.**

### Changes

- **`draws = old;`** When the file of the draws exists but does not hold the draws of the model (fewer or more values), AMLET now stops with a message. It used to ask on the standard input whether to overwrite the file, and overwrote it when there was no answer, as in a batch job. Without `drawsname`, `draws = old;` is an error; when the file does not exist, new draws are generated and saved in it, with a warning, as before.
- **Log-likelihood at constants.** It was printed as 0 (it was never computed). When every alternative is available to every observation, it is now computed: the probabilities of the model with alternative specific constants only are the shares of the alternatives in the choices (−283.7588 for the Australia data). Otherwise, the results say "not computed".
- **Likelihood ratio test.** The results gave a "level of rejection", 1 − p-value, which was 1 for any well-fitted model. They now give the p-value (from the upper incomplete gamma function, checked against scipy) and the degrees of freedom.

### Checks

- **`verify_likelihood.py`** runs AMLET, then recomputes with numpy and scipy, from the data and the draws that AMLET saves, the simulated log-likelihoods and the standard errors of three models: the Australia mixed logit (robust and classical standard errors), the panel model of `driver_panel` with 200 draws (five random coefficients, one of them lognormal; also the log-likelihood without the bias correction that the results give with `bias = 1`), and a multinomial logit with a fixed coefficient. Everything agrees, to 10⁻⁸ for the log-likelihoods and 3·10⁻⁶ for the standard errors. The CI runs it in the four build jobs.
- **`run_checks.py`** has 33 checks (log-likelihood at constants, p-value, availability rows in the Train96 format, `draws = old;`); v0.13.1-4 fails 7 of them. The availability of the alternatives (`censor = 1`), which no example used, was checked against an independent estimation.
- The standard errors of the Australia mixed logit examples are checked to within 3% (10% for `driver_mmnl_spline`): the point where the optimization stops varies slightly with the BLAS and the compiler.

### Other changes

- `CHANGELOG.md` gathers the notes of all the releases.
- `CITATION.cff` gives the Zenodo DOI and the ORCID of the author; the manual documents `draws` and `drawsname`.

The releases are archived on Zenodo: [doi:10.5281/zenodo.23454086](https://doi.org/10.5281/zenodo.23454086) designates all of them. See the [README](https://github.com/fbastin/amlet-legacy/blob/v0.13.1-5/README.md) for details and known limitations.

## v0.13.1-4 — 2026-10-18

Archived on Zenodo: [doi:10.5281/zenodo.23454087](https://doi.org/10.5281/zenodo.23454087).

AMLET 0.13.1 (2026 rebuild), fifth release. The simulated log-likelihood and the robust standard errors were recomputed independently, in Python, from the data and the draws that AMLET saves: they are right. This release fixes the results that were not, and the options that did nothing. **The standard errors of the mixed logit examples change**: they now use the exact Hessian (see the first item), and `simulated/driver_panel_halton` also changes its log-likelihood (see "Correction of the simulation bias").

### Results that change

- **Standard errors with `hessian_approx = true`.** This option computes the standard errors with the BFGS approximation of the Hessian that the optimization builds. The approximation depends on the optimization path, so the standard errors changed between builds, and it can be far from the Hessian: on `simulated/driver_panel`, the robust standard errors were 10 to 40% too small (0.0253 instead of 0.0425 for the first parameter; 0.0550 instead of 0.0661 on `swissmetro/driver_mmnl`). The default, `hessian_approx = false`, uses the exact Hessian, whose standard errors match an independent computation. AMLET now warns when `hessian_approx = true`, and the examples use the exact Hessian, except the two `_legacy` examples, which reproduce AMLET 0.13.1.
- **Fixed coefficients** (`FIXED` in `btype`). In a mixed logit, every t-statistic after a fixed coefficient was wrong: the t-statistic of the j-th estimated parameter divided the j-th value of the vector of all parameters, fixed ones included (276.4 instead of 3.54 on the Australia model). In a multinomial logit, fixed coefficients crashed AMLET; they are now supported, and the log-likelihood, estimates and robust standard errors match an independent computation. Without `b`, the fixed values were undefined; this is now an error.
- **Correction of the simulation bias** (`bias = 1`, the default). The optimization uses the simulated log-likelihood corrected for the simulation bias, and the results gave the corrected value without saying so. They now say it, and give the simulated log-likelihood without the correction (−5796.274 for `driver_panel`, whose results give −5791.972). The correction estimates the error of the simulated probabilities as for independent draws: with quasi-Monte Carlo draws, it is now turned off, with a warning. `driver_panel_halton` (200 Halton points) gives LL = −5802.493 instead of −5781.358, and the standard error of the first parameter is 0.0252 instead of 0.0388 (0.0253 for `driver_panel`, with 1000 Monte Carlo draws).
- **`level_accuracy`**, documented in the manual and used by the methods with a variable sample size, was not known to the parser, and was ignored. It is now read.
- **`draws = old;`** with `qmcrep` > 1 read only the draws of the first replication.

### Options that did nothing

`rescale`, `rescalidx`, `rescalval`, `ignore`, `synthetic` and `varrep` are read, but not implemented in AMLET 0.13.1. They now give a warning.

### Other fixes

- `alternatives_sampling` (experimental and not documented: it writes the data with sampled alternatives to `intermediaire`, then stops) swapped the arguments of the integer generator, which read past a variable on the stack and drew integers in [0, 1). With 3 alternatives or more it crashed or looped forever, and it freed memory twice. It now draws from `seed1`.
- Results file: a list of the estimates, without a title and left from debugging, preceded the table of estimates, and the header "Standard deviation" ran into the next one. The t-statistics are absolute values, and the header now says `|t-statistic|`.
- A warning of `validation = true` lacked its argument.

### Checks

The simulated log-likelihood of `australia/driver_mmnl` is identical, to all printed digits, to a computation in Python (numpy and scipy) from the data and the saved draws, and its robust standard errors (with the exact Hessian) agree to 5 digits; with `robust = false`, AMLET gives the inverse of the Hessian. The same holds for the panel model of `driver_panel` (five random coefficients, one of them lognormal), with 200 draws. On the Australia model, `mcbfgsda` reaches another local maximum of the same simulated likelihood (LL = −177.978, against −178.996 for `mcbtrda`): the sign of the standard deviation is not identified.

`run_checks.py` now has 26 checks; the previous release fails 7 of them. The other 18 examples give the same log-likelihoods and estimates as v0.13.1-3; the standard errors of the mixed logit examples change with the Hessian. The CI runs both on Ubuntu 24.04 (GCC 13) and 26.04 (GCC 15), with reference BLAS and OpenBLAS, and under AddressSanitizer and UBSan.

See the [README](https://github.com/fbastin/amlet-legacy/blob/v0.13.1-4/README.md) for details and known limitations.

## v0.13.1-3 — 2026-10-18

AMLET 0.13.1 (2026 rebuild), fourth release. The driver parser and the data readers did not check what they read: some mistakes in a driver or a data file gave wrong results without any message, and others crashed. This release reports them. **The results of correct drivers do not change**: the 19 examples give the same results as v0.13.1-2, except `driver_errc_import` (see below).

### Results that could be silently wrong

- **Values of the wrong type.** A value was stored as read, in a union, whatever the type of the parameter: `niter = 100.0;` gave 0 iterations, `nrep = 1e3;` was read as 0, and `eps = 1;` left `eps` almost unchanged. A real with an integer value is now accepted for an integer parameter, and an integer for a real one; any other mismatch is an error.
- **Unknown names**, such as `nitre = 100;`, were ignored silently. AMLET now prints a warning.
- **Train96 choices.** A choice row without a 1 left the choice undefined: the Australia data with one such record gave LL = −200.42 instead of −199.13, without any message. A row with several 1 kept the last one. Both are now errors.
- **Synthetic data.** `alttypes[i]` and `altvalues[i]` shorter than the model needs were read past their end, and gave other data. This is now an error.
- **Precision of the files.** Exported synthetic data had 6 significant digits, so estimating from them gave other results; they now have 17, and `driver_errc_import` gives exactly the results of `driver_errc`. The draws saved with `drawsname` had 16 digits: with `draws = old;`, the Australia mixed logit ended after 66 iterations instead of 70, with LL = −178.9961262 instead of −178.9961242. They now have 17, and the estimation is reproduced exactly.
- **Bias of simulation.** With the methods with a fixed sample size and without `statistics`, the results could print a "Bias of simulation" with an undefined value.

### Crashes, now errors

- `niter = {100};` (a list for a single value), `b = 0;` (a single value for a list), `b[1] = 0.5;`.
- `generate` before `nalt`, or `nalt` changed after `generate`.
- No `nalt`, `nalt = 1`, `nobs = 0`, `censor = 1` without a valid `idcensor`, and a negative `seed1` (which gave an invalid MRG32k3a state).
- In Train99, a choice outside 1 to `nalt`.
- `mcrtrda` when `sizes.txt` could not be created in the current directory.

`nests = 1;`, the name of `nnests` in AMLET ≤ 0.10, crashed the parser; it is now read as `nnests`, with a warning.

### Other fixes

- OPHELIA: in the active sets of the monotonicity constraints (`mcbtrc`), a misplaced parenthesis, `fabs(p <= DBL_EPSILON)`, accepted any negative direction as zero. The spline example gives the same results.
- `input "file";` opened the file and did nothing; it is now ignored with a warning.

### Tests

`run_checks.py` (19 checks, a few seconds) changes the drivers or the data of the examples, and checks that each error is reported with a message and a non-zero exit code, and that `draws = old;` reproduces an estimation. The previous release fails 18 of them. The CI runs them in every job, including under AddressSanitizer and UBSan.

See the [README](https://github.com/fbastin/amlet-legacy/blob/v0.13.1-3/README.md) for details and known limitations.

## v0.13.1-2 — 2026-10-18

AMLET 0.13.1 (2026 rebuild), third release. Every option of the driver file was run under AddressSanitizer and UndefinedBehaviorSanitizer; this release fixes the options that did not work. **Several of them change results**: see the first section if you used quasi-Monte Carlo draws, `robust = false`, `validation`, `mlhs`, `lattice` or `randomization`.

### Results that change

- **Quasi-Monte Carlo draws with a variable sample size.** `mcbtrda` (the default method) and `mcrtrda` estimate the simulation error as for independent draws. With Halton, Sobol, lattice or MLHS draws, they stalled (and, since `niter` was ignored, iterated forever) or stopped with as few as 32 draws. AMLET now uses `mcbtr` or `mcrtr` instead, with a warning.
- **`mlhs`** gave Monte Carlo draws; it now uses ORATIO's modified Latin hypercube sampling (Hess et al.).
- **`lattice`** used the generating vector of the largest tabulated lattice not greater than `nrep` (887 points for 1000) with `nrep` points. `nrep` is now set to the number of points of that lattice.
- **`randomization`** was ignored. All individuals use the same Halton points, each with its own random shift; `scrambling` now adds the Braaten-Weller permutations, and `shuffling` shuffles the order of the points of each individual.
- **`robust = false`** printed the robust covariance matrix, labelled "Standard estimation". It now prints the inverse of the Hessian.
- **`validation = true`** was always disabled: the driver checks ran before AMLET knew the number of random coefficients, so every model was checked as a multinomial logit.
- **Truncated normal coefficients** (`truncated_normal`, `truncated_normal_b` in `btype`) were estimated wrongly (no derivatives; uniform draws for `truncated_normal_b`). They are now an error in `btype`, and remain available in `bgentype` to generate synthetic data.

With `legacy_rng = true`, the MRG32k3a draws and the random shifts of AMLET 0.13.1 are still reproduced, and `randomization` is still ignored; the other fixes above apply too (fixed sample size with quasi-Monte Carlo draws, MLHS, number of lattice points).

### Errors that were not reported

- `nnests` other than 1 (nested logit models are not supported), `mcbtrc` without spline coefficients (now `mcbtr`), and a `btype`, `bgentype` or `bgen` of the wrong length wrote past their arrays or crashed.
- In the Train96 format, a missing `iddep` made AMLET read the choice row as a variable; truncated Train96 and Train99 files gave undefined values.
- A `b` of the wrong length ended with "normal termination" and exit code 0.

AMLET now stops with a message and a non-zero exit code.

### Other fixes

- `niter` was ignored by `mcbtrda` and `mcrtrda`.
- The results printed `mcrtrda` for `mcbtrda`, and vice versa.
- `faure`, not implemented, gave Monte Carlo draws; it is now an error.
- The `ecuyer` generator (L'Ecuyer 1988 with a shuffle table) crashed; it now gives the draws of Numerical Recipes' `ran2`.
- ORATIO: an integer overflow in the padding of the random shifts of quasi-Monte Carlo points; scrambled Halton sequences beyond dimension 16 read past their table; three unused functions called themselves forever or ignored their increment.
- Messages longer than 1023 characters were truncated.

### Tests

`run_examples.sh` now runs 19 examples (about 1.5 minutes): the quasi-Monte Carlo methods, a uniform coefficient with the classical covariance, a spline coefficient with monotonicity constraints, the `ecuyer` generator, two other optimization methods, scrambled Halton draws on a panel, and the estimation from the data that a synthetic example exports. The CI runs them on Ubuntu 24.04 (GCC 13) and 26.04 (GCC 15), with reference BLAS and OpenBLAS, and under AddressSanitizer and UBSan. The examples of v0.13.1-1 give the same results; only the name of the method changes in their results files (`mcbtrda`, previously printed as `mcrtrda`).

See the [README](https://github.com/fbastin/amlet-legacy/blob/v0.13.1-2/README.md) for details and known limitations.

## v0.13.1-1 — 2026-10-17

AMLET 0.13.1 (2026 rebuild), second release. It builds with GCC 15 (Ubuntu 26.04) and fixes two bugs that affect results: please read the first two points if you used v0.13.1.

### Results that change

- **Standard errors could be wrong.** ORATIO's `f2c.h` declared Fortran integers as `long` (8 bytes), while gfortran's `INTEGER` has 4. After the Cholesky factorisation of the Hessian, half of LAPACK's `info` was left from the stack: when it was not zero, a successful factorisation was taken for a failure, and only the diagonal of the Hessian was inverted. The standard errors were then wrong (up to 5 times too small in the mixed logit examples), while the estimates and the log-likelihood were right. This depends on the compiler and its options: the v0.13.1 builds with GCC 13 and the default options were not affected, but GCC 15 on Ubuntu 26.04 was. Fixed.
- **The random number generator is now L'Ecuyer's MRG32k3a.** ORATIO's MRG32k3a used x(n−1) instead of x(n−2) in its first component, and its jumps to the next stream or substream were not reduced modulo m. Both are fixed: 10^6 draws are identical to those of RandomDataStreams.jl, and streams and substreams agree with R's `L'Ecuyer-CMRG`. **Results that depend on random draws (mixed logit, synthetic data) therefore differ from v0.13.1.** To reproduce them exactly, add `legacy_rng = true;` to the driver file.

### Other changes

- Builds with GCC 14 and 15 (C23): missing prototypes (`yylex()`, `nlp_set_derivatives_type()`), a `config.guess` from about 2000 replaced in ORATIO and OPHELIA, and their autotools scripts made executable (needed on Ubuntu 26.04, whose `mkdir` comes from uutils).
- Memory errors found with AddressSanitizer: two file-name buffers one byte too short, and a read past the array of draws (unused value, no effect on results).
- `build.sh` always compiles with `-ffp-contract=off`, so that `legacy_rng = true` gives the same draws with or without FMA instructions.
- The driver accepts `rangen = mrg32k3a;`, which was missing.
- `run_examples.sh` runs 8 examples (two of them with `legacy_rng = true`) and checks the log-likelihood and a standard error.
- CI on Ubuntu 24.04 (GCC 13) and 26.04 (GCC 15), with reference BLAS and OpenBLAS, plus a job under AddressSanitizer and UBSan.
- The patches against the original tarballs cover AMLET, OPHELIA and ORATIO. The AMLET manual documents `legacy_rng`.

### Quick start

```sh
sudo apt install build-essential libblas-dev liblapack-dev python3 curl
bash build.sh          # installs into ./local
bash run_examples.sh   # 8 examples, about 45 seconds
```

See the [README](https://github.com/fbastin/amlet-legacy/blob/v0.13.1-1/README.md) for details and known limitations.

## v0.13.1 — 2026-10-17

The last release of AMLET (*Another Mixed Logit Estimation Tool*), with the ORATIO 0.7.5 and OPHELIA 0.7.4 libraries, made to build and run on current Linux systems (tested with GCC 13.3 on Ubuntu 24.04). The code is no longer actively developed.

### Changes with respect to amlet-0.13.1.tar.gz (2015)

- **Builds with current compilers.** The 2015 tarball did not compile. Fixed: missing `src/exec.h`, undefined `_max` / `_min`, global variables defined in headers (GCC ≥ 10), C99 `inline` semantics, missing prototypes.
- **Bias correction fixed.** `simulation_bias()` was implicitly declared as `int`, so in builds of the 0.13.x tarballs the bias correction (`bias = 1`) used a garbage value. Versions ≤ 0.12 are not affected.
- **Panel data** (`np < nobs`): the "Total" log-likelihood, ρ² and the likelihood-ratio statistic now use the correct count. Parameter estimates were never affected.
- **`build.sh`** builds the three packages into a local prefix, without root access. Bison/flex, autotools and LaTeX are not needed.
- **`run_examples.sh`** runs 6 examples (Australia, simulated data, Swissmetro) and checks them against reference values. The Swissmetro logit reproduces Biogeme's `01logit` example.
- **License**: MIT for the files whose only copyright holder is Fabian Bastin. Files co-authored with others keep the Open Software License v2.1 (see the README).

### Quick start

```sh
sudo apt install build-essential libblas-dev liblapack-dev python3 curl
bash build.sh          # installs into ./local
bash run_examples.sh   # 6 examples, about 40 seconds
```

See the [README](https://github.com/fbastin/amlet-legacy/blob/v0.13.1/README.md) for details and known limitations, and `doc/amlet_manual.pdf` for the driver file format.
