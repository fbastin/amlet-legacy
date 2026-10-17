#!/usr/bin/env python3
"""
Convert the public Swissmetro SP data (Bierlaire, Axhausen & Abay, 2001) to
the AMLET "Train96" format, using the specification of Biogeme's reference
example 01logit (logit with ASC_TRAIN, ASC_CAR, B_TIME, B_COST).

The data are NOT shipped here; download them from the Biogeme website:
    https://transp-or.epfl.ch/data/swissmetro.dat

Usage:
    python3 swissmetro_to_amlet.py swissmetro.dat swissmetro_train96

Sample selection (as in Biogeme): PURPOSE in {1, 3} (commuters, business)
and CHOICE != 0. Alternatives: 1 = train, 2 = Swissmetro, 3 = car.

Each observation is written as a block of 6 rows x 3 columns:
    row 1  ASC_TRAIN dummy         1        0       0
    row 2  ASC_CAR dummy           0        0       1
    row 3  travel time / 100       TRAIN_TT SM_TT   CAR_TT
    row 4  cost / 100              TRAIN_CO SM_CO   CAR_CO   (0 if GA owner, for train and SM)
    row 5  availability            TRAIN_AV*(SP!=0) SM_AV CAR_AV*(SP!=0)
    row 6  choice (dependent)      0/1 indicator
so the AMLET driver uses nvar = 6, idcensor = 5, iddep = 6.
"""

import csv
import sys


def main(src, dst):
    with open(src, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    n = 0
    with open(dst, "w") as out:
        for r in rows:
            v = {k: float(x) for k, x in r.items()}
            if v["PURPOSE"] not in (1.0, 3.0) or v["CHOICE"] == 0:
                continue
            ga = v["GA"] != 0
            sp = v["SP"] != 0
            block = [
                [1, 0, 0],
                [0, 0, 1],
                [v["TRAIN_TT"] / 100, v["SM_TT"] / 100, v["CAR_TT"] / 100],
                [0 if ga else v["TRAIN_CO"] / 100,
                 0 if ga else v["SM_CO"] / 100,
                 v["CAR_CO"] / 100],
                [int(v["TRAIN_AV"] and sp), int(v["SM_AV"]), int(v["CAR_AV"] and sp)],
                [int(v["CHOICE"] == k) for k in (1, 2, 3)],
            ]
            for line in block:
                out.write(" ".join(f"{x:g}" for x in line) + "\n")
            n += 1

    print(f"{n} observations written to {dst}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
