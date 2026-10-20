#!/usr/bin/env python3
"""
Rebuild pickle cache from HyPhy JSON results.
This is useful when pickles become incompatible with a new pandas version.
"""

import argparse
import sys
from pathlib import Path
from hyphy_results_parser import HyphyResultsManager

# Define paths
REPO_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = REPO_ROOT / "results" / "hyphy_results_cache"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--suffix', default='',
                    help="Tag inserted before .pkl, e.g. 'python3.13' -> relax_results.python3.13.pkl")
args = parser.parse_args()
SUFFIX = f".{args.suffix.lstrip('.')}" if args.suffix else ""

if not RESULTS_DIR.exists():
    print(f"Error: Results directory not found: {RESULTS_DIR}")
    sys.exit(1)

print(f"Rebuilding pickle cache from JSONs in {RESULTS_DIR}")
print("-" * 60)

manager = HyphyResultsManager()

# Load RELAX results from JSON subdirectory
print("Loading RELAX results...")
try:
    relax_json_dir = RESULTS_DIR / "relax"
    if not relax_json_dir.exists():
        raise FileNotFoundError(f"RELAX JSON directory not found: {relax_json_dir}")
    relax_result = manager.load_relax_from_json(str(relax_json_dir), name='relax')
    relax_pkl = RESULTS_DIR / f"relax_results{SUFFIX}.pkl"
    relax_result.save_to_pickle(str(relax_pkl))
    print(f"✓ Saved RELAX results to {relax_pkl}")
except Exception as e:
    print(f"✗ Failed to load/save RELAX results: {e}")
    sys.exit(1)

# Load BUSTED-PH (orb) results from JSON subdirectory
print("Loading BUSTED-PH (orb) results...")
try:
    busted_orb_json_dir = RESULTS_DIR / "busted_ph_orb"
    if not busted_orb_json_dir.exists():
        raise FileNotFoundError(f"BUSTED-PH (orb) JSON directory not found: {busted_orb_json_dir}")
    busted_ph_orb = manager.load_busted_ph_from_json(str(busted_orb_json_dir), name='busted_ph_orb')
    busted_pkl = RESULTS_DIR / f"busted_ph_orb_results{SUFFIX}.pkl"
    busted_ph_orb.save_to_pickle(str(busted_pkl))
    print(f"✓ Saved BUSTED-PH (orb) results to {busted_pkl}")
except Exception as e:
    print(f"✗ Failed to load/save BUSTED-PH (orb) results: {e}")
    sys.exit(1)

# Load BUSTED-PH (non-orb) results from JSON subdirectory
print("Loading BUSTED-PH (non-orb) results...")
try:
    busted_non_orb_json_dir = RESULTS_DIR / "busted_ph_non_orb"
    if not busted_non_orb_json_dir.exists():
        raise FileNotFoundError(f"BUSTED-PH (non-orb) JSON directory not found: {busted_non_orb_json_dir}")
    busted_ph_non_orb = manager.load_busted_ph_from_json(str(busted_non_orb_json_dir), name='busted_ph_non_orb')
    busted_pkl = RESULTS_DIR / f"busted_ph_non_orb_results{SUFFIX}.pkl"
    busted_ph_non_orb.save_to_pickle(str(busted_pkl))
    print(f"✓ Saved BUSTED-PH (non-orb) results to {busted_pkl}")
except Exception as e:
    print(f"✗ Failed to load/save BUSTED-PH (non-orb) results: {e}")
    sys.exit(1)

print("-" * 60)
print("✓ All pickles rebuilt successfully!")

