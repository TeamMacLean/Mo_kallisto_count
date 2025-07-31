#!/usr/env python

import sys
from collections import defaultdict
import os

files = sys.argv[1:]
samples = [os.path.basename(os.path.dirname(os.path.dirname(f))) for f in files]
targets = defaultdict(list)


def get_vals(file):
    v = []
    with open(file, "r") as f:
        for line in f.readlines():
            line=line.strip()
            if not line.startswith("target"):
                els = line.split("\t")
                v.append([els[0], els[4]])
    return v

for f in files:
    vals = get_vals(f)
    for v in vals:
        targets[v[0]].append(v[1])

print("target,{}".format(",".join(samples)))

for key, el in targets.items():
    print("{},{}".format(key, ",".join(el)))

