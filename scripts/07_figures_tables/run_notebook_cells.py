"""Run chosen cells of a notebook headlessly in the project's Jupyter kernel, so the notebook does not have to be reopened and rerun by hand.

Only the cells given in --cells are executed (in the order given, in one kernel, so later cells can use the variables of earlier ones). The notebook file is
never modified; the executed copy, with its outputs, goes to --out. Figures are saved by the cells themselves, so --pres redirects them to another folder
(it overrides the notebook's `pres` right after the setup cell, so nothing in the notebook's usual output folder is overwritten).

Example (the final odds ratio figures of 'Odds Ratio Test Plots.ipynb'):
  python scripts/07_figures_tables/run_notebook_cells.py --cells 3,5,17,21,30,32 --noop-save 17 \\
      --pres "/Users/calvin/.../talks/ort_final_nov17" --out /tmp/executed.ipynb
"""
import argparse
import os
import time

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--notebook", default=os.path.join(HERE, "Odds Ratio Test Plots.ipynb"))
    ap.add_argument("--cells", required=True, help="comma-separated cell indices (0-based, as in the notebook), run in this order")
    ap.add_argument("--pres", default=None, help="folder for the saved figures (sets `pres` after the setup cell that defines it)")
    ap.add_argument("--noop-save", default="", help="comma-separated cell indices whose figures are NOT saved (cell is run only for the variables it defines)")
    ap.add_argument("--out", required=True, help="where to write the executed copy of the selected cells")
    ap.add_argument("--kernel", default="orb-selection")
    ap.add_argument("--timeout", type=int, default=3600)
    args = ap.parse_args()

    nb = nbformat.read(args.notebook, as_version=4)
    wanted = [int(c) for c in args.cells.split(",")]
    noop = {int(c) for c in args.noop_save.split(",") if c.strip()}
    cells, labels = [], []
    for idx in wanted:
        src = nb.cells[idx]
        if src.cell_type != "code":
            raise SystemExit(f"cell {idx} is not a code cell")
        if idx in noop:
            cells.append(nbformat.v4.new_code_cell("_real_save = save\nsave = lambda *a, **k: None")); labels.append(f"(mute save before {idx})")
        cells.append(nbformat.v4.new_code_cell(src.source)); labels.append(f"cell {idx}")
        if idx in noop:
            cells.append(nbformat.v4.new_code_cell("save = _real_save")); labels.append(f"(restore save after {idx})")
        if args.pres and "pres = " in src.source:
            cells.append(nbformat.v4.new_code_cell(f"pres = {args.pres!r}\nos.makedirs(pres, exist_ok=True)\nprint('figures go to', pres)")); labels.append("(set output folder)")
    run = nbformat.v4.new_notebook(cells=cells, metadata=nb.metadata)

    client = NotebookClient(run, kernel_name=args.kernel, timeout=args.timeout, resources={"metadata": {"path": os.path.dirname(os.path.abspath(args.notebook))}}, allow_errors=False)
    t0 = time.time()
    try:
        client.execute()
        status = "ok"
    except CellExecutionError as e:
        status = f"FAILED: {str(e)[-1500:]}"
    nbformat.write(run, args.out)
    for lab, c in zip(labels, run.cells):
        for o in c.get("outputs", []):
            if o.get("output_type") == "stream":
                print(f"[{lab}]", o["text"].rstrip()[-400:])
    print(f"{status} ({time.time() - t0:.0f} s); executed copy: {args.out}")


if __name__ == "__main__":
    main()
