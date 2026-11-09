from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT.parent / "benchmark"
for p in (str(ROOT), str(BENCH), str(ROOT.parent)):
    if p not in sys.path:
        sys.path.insert(0, p)
