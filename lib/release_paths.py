"""V52-R1 result roots; paths are independent of the working directory.
V52_RESULTS_ROOT is set by run_all.py --output. It is not an input-data source.
"""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RESULTS=Path(os.environ.get('V52_RESULTS_ROOT',str(ROOT/'results/current'))).resolve()
