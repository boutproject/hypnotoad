#!/usr/bin/env python3

from boututils.run_wrapper import shell
from glob import glob
from pathlib import Path
from sys import exit

tests_dir = str(Path(__file__).parent)
test_dirs = list(glob(tests_dir + "/*/"))

# no test in grid_files directory
skip_dirs = ["grid_files", "__pycache__"]
test_dirs = [d for d in test_dirs if not any(x in d for x in skip_dirs)]

results = []

for d in test_dirs:
    name = str(Path(d).name)
    print("\nrunning", name, "...\n", flush=True)
    retcode, _ = shell(d + "/runtest.py")
    results.append((name, retcode))

for r in results:
    if r[1] == 0:
        print(r[0], "passed")
    else:
        print(r[0], "failed")

if all([r[1] == 0 for r in results]):
    exit(0)
else:
    exit(1)
