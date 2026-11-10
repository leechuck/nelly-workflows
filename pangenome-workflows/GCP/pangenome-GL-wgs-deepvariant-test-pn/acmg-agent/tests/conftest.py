from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT.parent / "benchmark"
for p in (str(ROOT), str(BENCH), str(ROOT.parent)):
    if p not in sys.path:
        sys.path.insert(0, p)
for sp in ROOT.glob(".venv/lib/python*/site-packages"):
    sp_str = str(sp)
    if sp_str not in sys.path:
        sys.path.append(sp_str)
