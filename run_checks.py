#!/usr/bin/env python3
"""Check that AMLET reports the errors of driver and data files.

Each check copies an example, changes its driver or its data, runs AMLET,
and checks the exit code and a message of the output. Any report of
AddressSanitizer or UndefinedBehaviorSanitizer is a failure too.

    python3 run_checks.py                      # uses ./local/bin/amlet
    AMLET=$HOME/amlet/bin/amlet python3 run_checks.py
    RUN=/tmp/amlet-checks python3 run_checks.py  # where outputs are written
"""
import os, re, shutil, subprocess, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
AMLET = os.path.abspath(os.environ.get("AMLET", os.path.join(ROOT, "local", "bin", "amlet")))
RUN = os.path.abspath(os.environ.get("RUN", os.path.join(ROOT, "run", "checks")))
MNL_LL = "-199.1283687"   # examples/reference_values.tsv, australia/driver_mnl


def sub(old, new):
    """Replace the text old, which must be in the driver, by new."""
    def edit(text):
        if old not in text:
            raise ValueError("%r is not in the driver" % old)
        return text.replace(old, new, 1)
    return edit


def statement(name, new=""):
    """Replace the whole statement 'name = {...};' by new."""
    pattern = re.compile(re.escape(name) + r"\s*=\s*\{[^}]*\};")
    def edit(text):
        if not pattern.search(text):
            raise ValueError("no statement %s in the driver" % name)
        return pattern.sub(new, text, 1)
    return edit


def both(*edits):
    def edit(text):
        for e in edits:
            text = e(text)
        return text
    return edit


def no_choice(d):
    """Train96: record 3 of the Australia data without a chosen alternative
    (7 rows per record; the choice is the 7th)."""
    path = os.path.join(d, "model_australia")
    lines = open(path).read().split("\n")
    lines[3*7 - 1] = "0 0 0 0"
    open(path, "w").write("\n".join(lines))


def tiny_train99(d):
    """Train99 data, 2 individuals and 2 alternatives, with a choice of 3."""
    open(os.path.join(d, "tiny_x"), "w").write("1 2\n2 1\n")
    open(os.path.join(d, "tiny_y"), "w").write("1\n3\n")
    open(os.path.join(d, "tiny_times"), "w").write("1\n1\n")
    open(os.path.join(d, "driver_tiny"), "w").write(
        'output "results_tiny";\nformat = Train99;\nnobs = 2;\nnp = 2;\n'
        'nalt = 2;\nnvar = 1;\nbtype = { CONSTANT };\nb = { 0 };\n'
        'method = MCBTR;\nload "tiny";\n')


# name, example, driver, edit of the driver, edit of the data,
# expected exit ("ok" or "error"), message expected in the output (or, after
# "!", a message that must not be there), log-likelihood expected in the
# results file (or None)
CHECKS = [
    ("integral_real", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100.0;"), None,
     "ok", None, MNL_LL),
    ("real_for_integer", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100.5;"), None,
     "error", "niter expects an integer", None),
    ("string_for_integer", "australia", "driver_mnl",
     sub("method = MCBTR;", 'method = "mcbtr";'), None,
     "error", "method expects an integer", None),
    ("list_for_scalar", "australia", "driver_mnl",
     sub("niter = 100;", "niter = {100};"), None,
     "error", "niter expects a single value", None),
    ("scalar_for_list", "australia", "driver_mnl",
     sub("b = { 0.0, 0., 0., 0., 0., 0. };", "b = 0;"), None,
     "error", "b expects a list of values", None),
    ("indexed_parameter", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100; b[1] = 0.5;"), None,
     "error", "b cannot be indexed", None),
    ("unknown_name", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100; nitre = 3;"), None,
     "ok", "unknown parameter nitre is ignored", MNL_LL),
    ("old_nests", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100; nests = 1;"), None,
     "ok", "nests is the name of nnests", MNL_LL),
    ("no_nalt", "australia", "driver_mnl",
     sub("nalt = 4;", ""), None,
     "error", "nalt (the number of alternatives) must be given", None),
    ("no_observation", "australia", "driver_mnl",
     both(sub("nobs = 210;", "nobs = 0;"), sub("np = 210;", "np = 0;")), None,
     "error", "nobs (the number of observations) must be at least 1", None),
    ("censor_without_idcensor", "australia", "driver_mnl",
     sub("censor = 0;", "censor = 1;"), None,
     "error", "idcensor must give the row of the availabilities", None),
    ("negative_seed", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100; seed1 = -5;"), None,
     "error", "seed1 must be positive", None),
    ("train96_no_choice", "australia", "driver_mnl",
     None, no_choice,
     "error", "record 3: the choice row (row 7) must contain one 1", None),
    ("train99_bad_choice", "australia", "driver_tiny",
     None, tiny_train99,
     "error", "observation 2: the choice must be an alternative, from 1 to 2", None),
    ("alttypes_too_short", "simulated", "driver_errc",
     statement("alttypes[2]", "alttypes[2] = { NORMAL, NORMAL, NORMAL, NORMAL };"), None,
     "error", "alttypes[2] has 4 values; the model has 5 variables", None),
    ("altvalues_too_short", "simulated", "driver_errc",
     statement("altvalues[3]", "altvalues[3] = { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };"), None,
     "error", "altvalues[3] has 8 values; alttypes[3] needs 10", None),
    ("missing_alternative", "simulated", "driver_errc",
     both(statement("alttypes[4]"), statement("altvalues[4]")), None,
     "error", "missing definition for alternative 4", None),
    ("fixed_mnl", "australia", "driver_mnl",
     both(sub("btype = { CONSTANT,", "btype = { FIXED,"),
          sub("b = { 0.0,", "b = { 4.0,")), None,
     "ok", None, "-200.4075665"),
    ("fixed_without_b", "australia", "driver_mnl",
     both(sub("btype = { CONSTANT,", "btype = { FIXED,"),
          sub("b = { 0.0, 0., 0., 0., 0., 0. };", "")), None,
     "error", "b must be given, with the values of the fixed coefficients", None),
    ("level_accuracy", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100; level_accuracy = 0.05;"), None,
     "ok", "!unknown parameter", MNL_LL),
    ("not_implemented", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100; ignore = {2};"), None,
     "ok", "ignore is not implemented", MNL_LL),
    ("qmc_bias", "australia", "driver_mmnl_halton",
     both(sub("bias=0;", "bias=1;"), sub("nrep = 1000;", "nrep = 100;"),
          sub("niter = 500;", "niter = 3;")), None,
     "ok", "the correction of the simulation bias requires Monte Carlo draws", None),
    ("alternatives_sampling", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100; alternatives_sampling = 3;"), None,
     "ok", "Alternatives sampling completed", None),
    ("too_many_alternatives", "australia", "driver_mnl",
     sub("niter = 100;", "niter = 100; alternatives_sampling = 5;"), None,
     "error", "from 2 to nalt = 4 alternatives can be sampled", None),
    ("generate_before_nalt", "simulated", "driver_errc",
     both(sub("nalt = 5;", ""),
          sub('generate "errc_data";', 'generate "errc_data"; nalt = 5;')), None,
     "error", "nalt must be given before generate", None),
]

SANITIZERS = re.compile(r"AddressSanitizer|LeakSanitizer|runtime error:")


def run(d, driver):
    p = subprocess.run([AMLET, driver], cwd=d, stdin=subprocess.DEVNULL,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       timeout=1200)
    out = p.stdout.decode("latin-1")
    open(os.path.join(d, driver + ".log"), "w", encoding="latin-1").write(out)
    return p.returncode, out


def results(d, driver):
    """The 'Optimal' line and the number of iterations of the results file."""
    text = open(os.path.join(d, driver), encoding="latin-1").read()
    m = re.search(r'^output\s+"([^"]+)"', text, re.M)
    try:
        res = open(os.path.join(d, m.group(1)), encoding="latin-1").read()
    except (AttributeError, OSError):
        return None, None
    opt = re.search(r"^Optimal\s+(\S+)", res, re.M)
    it = re.search(r"^Number of iterations: (\d+)", res, re.M)
    return (opt.group(1) if opt else None), (it.group(1) if it else None)


def prepare(name, example):
    d = os.path.join(RUN, name)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(os.path.join(ROOT, "examples", example), d)
    return d


def check(name, example, driver, edit, data, expect, message, ll):
    d = prepare(name, example)
    if data:
        data(d)
    if edit:
        path = os.path.join(d, driver)
        text = open(path, encoding="latin-1").read()
        open(path, "w", encoding="latin-1").write(edit(text))
    rc, out = run(d, driver)
    problems = []
    if SANITIZERS.search(out):
        problems.append("sanitizer report")
    if expect == "ok" and rc != 0:
        problems.append("exit code %d" % rc)
    if expect == "error" and rc == 0:
        problems.append("exit code 0")
    if message and message.startswith("!"):
        if message[1:] in out:
            problems.append("unexpected message")
    elif message and message not in out:
        problems.append("message not found")
    if ll and results(d, driver)[0] != ll:
        problems.append("LL %s, expected %s" % (results(d, driver)[0], ll))
    return problems


def check_draws_old():
    """Draws saved with drawsname and read back with draws = old give the
    same estimation."""
    d = prepare("draws_old", "australia")
    problems = []
    rc1, out1 = run(d, "driver_mmnl")
    first = results(d, "driver_mmnl")
    path = os.path.join(d, "driver_mmnl")
    text = open(path, encoding="latin-1").read()
    if "#draws = old;" not in text:
        return ["driver_mmnl has no '#draws = old;' line"]
    open(path, "w", encoding="latin-1").write(text.replace("#draws = old;", "draws = old;"))
    rc2, out2 = run(d, "driver_mmnl")
    second = results(d, "driver_mmnl")
    if SANITIZERS.search(out1 + out2):
        problems.append("sanitizer report")
    if rc1 or rc2:
        problems.append("exit codes %d, %d" % (rc1, rc2))
    if first != second or None in first:
        problems.append("LL and iterations %s, then %s with draws = old" % (first, second))
    return problems


def main():
    if not os.access(AMLET, os.X_OK):
        sys.exit("ERROR: AMLET executable not found at %s (run build.sh, or set AMLET=...)" % AMLET)
    os.makedirs(RUN, exist_ok=True)
    fmt = "%-26s %-6s %s"
    print(fmt % ("check", "exit", "status"))
    failed = 0
    for name, example, driver, edit, data, expect, message, ll in CHECKS:
        problems = check(name, example, driver, edit, data, expect, message, ll)
        print(fmt % (name, expect, "FAILED (%s)" % ", ".join(problems) if problems else "OK"))
        failed += bool(problems)
    problems = check_draws_old()
    print(fmt % ("draws_old", "ok", "FAILED (%s)" % ", ".join(problems) if problems else "OK"))
    failed += bool(problems)
    total = len(CHECKS) + 1
    print("\n%d passed, %d failed. Outputs in %s" % (total - failed, failed, RUN))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
