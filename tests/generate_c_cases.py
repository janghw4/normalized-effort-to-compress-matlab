"""Regenerate C parity fixtures: python generate_c_cases.py /path/to/libetc.so"""
import ctypes
import csv
import hashlib
import json
from pathlib import Path
import random
import sys

library_path = Path(sys.argv[1]).resolve()
lib = ctypes.CDLL(str(library_path))
lib.compute_etc.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.c_int]
lib.compute_etc.restype = ctypes.c_int
rng = random.Random(20261022)
rows = []
for k in [2, 3, 5, 10, 20, 50]:
    for length in [20, 100, 500, 2000]:
        for trial in range(5):
            sequence = [rng.randrange(k) for _ in range(length)]
            value = lib.compute_etc((ctypes.c_int * length)(*sequence), length)
            assert value >= 0
            rows.append([' '.join(map(str, sequence)), value])
out = Path(__file__).resolve().parent
with (out/'c_reference_cases.csv').open('w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f); writer.writerow(['sequence', 'etc']); writer.writerows(rows)
metadata = {'seed': 20261022, 'cases': len(rows),
            'library_sha256': hashlib.sha256(library_path.read_bytes()).hexdigest(),
            'fixture_sha256': hashlib.sha256((out/'c_reference_cases.csv').read_bytes()).hexdigest()}
(out/'c_reference_provenance.json').write_text(json.dumps(metadata, indent=2)+'\n', encoding='utf-8')
print(json.dumps(metadata))
