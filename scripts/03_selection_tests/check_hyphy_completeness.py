# Run using the orb-selection conda environment:
# conda run -n orb-selection python scripts/03_selection_tests/check_hyphy_completeness.py

"""Pre-flight check for HyPhy JSON output directories.

Reuses HyphyResultsManager's own per-file row-building methods
(_parse_relax_record / _parse_busted_ph_record) to check that every JSON in a
results directory has the fields the real parser needs -- without building a
DataFrame, running LOC conversion, or saving pickles. Meant to be run right
after syncing new HyPhy output down from the cluster, to catch incomplete
files (e.g. a timed-out BUSTED-PH task left holding a plain BUSTED fit) before
spending time on the full parse_hyphy_results.py pipeline.

Exits 1 if any directory has incomplete files, 0 otherwise.
"""

import json
import os
import sys
import importlib.util
from pathlib import Path

repo_root = Path(__file__).parent.parent.parent
src_dir = repo_root / "src"
stage_dir = Path(__file__).parent
for path in (src_dir, stage_dir):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def _load_module(module_name: str, module_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(module_path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module {module_name} from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


hyphy_results_parser = _load_module(
    "hyphy_results_parser", stage_dir / "hyphy_results_parser.py"
)
HyphyResultsManager = hyphy_results_parser.HyphyResultsManager


def check_directory(json_directory: Path, parse_record):
    """Try to parse every JSON file in json_directory; return (n_ok, failures)."""
    json_files = sorted(f for f in os.listdir(json_directory) if f.endswith('.json'))
    failures = []
    n_ok = 0
    for json_file in json_files:
        with open(json_directory / json_file, 'r') as f:
            try:
                data = json.load(f)
            except json.JSONDecodeError as e:
                failures.append((json_file, f"invalid JSON: {e}"))
                continue
        try:
            parse_record(json_file, data)
            n_ok += 1
        except KeyError as e:
            failures.append((json_file, str(e)))
    return n_ok, failures


def main():
    results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else repo_root / "results" / "hyphy_results_cache"
    manager = HyphyResultsManager()

    checks = [
        ('relax', results_dir / 'relax', manager._parse_relax_record),
        ('busted_ph_orb', results_dir / 'busted_ph_orb', manager._parse_busted_ph_record),
        ('busted_ph_non_orb', results_dir / 'busted_ph_non_orb', manager._parse_busted_ph_record),
    ]

    print(f"=== HyPhy Results Completeness Check ({results_dir}) ===\n")

    any_failures = False
    for label, directory, parse_record in checks:
        if not directory.exists():
            print(f"✗ {label}: {directory} not found, skipping\n")
            continue

        n_ok, failures = check_directory(directory, parse_record)
        total = n_ok + len(failures)
        if total == 0:
            print(f"? {label}: {directory} has no JSON files yet, skipping")
        else:
            status = "✓" if not failures else "⚠"
            print(f"{status} {label}: {n_ok}/{total} complete")
            if failures:
                any_failures = True
                for json_file, reason in failures:
                    print(f"    - {json_file}: {reason}")
        print()

    print("=== Check Complete ===")
    sys.exit(1 if any_failures else 0)


if __name__ == "__main__":
    main()
