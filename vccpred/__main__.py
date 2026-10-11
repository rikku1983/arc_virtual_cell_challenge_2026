"""Predict with a registered model and package the submission.

  .venv/bin/python -m vccpred --model control_resample
  .venv/bin/python -m vccpred --model control_resample --param residual=0.3 --tag r30
  .venv/bin/python -m vccpred --model control_resample --limit-perts 3   # quick smoke test, no packaging
"""
import argparse
import ast
import gc
import json
import time
from pathlib import Path

from vccpred.data import Bundle
from vccpred.models import MODELS, get_model
from vccpred.submission import vcc_prep, write_prediction

OUT_DIR = Path(__file__).resolve().parents[1] / "predictions"


def parse_param(s: str) -> tuple[str, object]:
    key, _, value = s.partition("=")
    try:
        return key, ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return key, value  # bare string


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, choices=sorted(MODELS))
    ap.add_argument("--param", action="append", default=[], type=parse_param, metavar="KEY=VALUE",
                    help="model constructor argument (Python literal), repeatable")
    ap.add_argument("--tag", default="", help="suffix for the run name")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--limit-perts", type=int, default=None, help="predict only the first N perturbations; skips packaging")
    ap.add_argument("--no-package", action="store_true", help="write the .h5ad only")
    args = ap.parse_args()

    bundle = Bundle(args.data_dir) if args.data_dir else Bundle()
    params = dict(args.param)
    model = get_model(args.model, **params)
    perts = bundle.perts[: args.limit_perts] if args.limit_perts else bundle.perts

    run = args.model + (f"_{args.tag}" if args.tag else "") + ("_smoke" if args.limit_perts else "")
    run_dir = OUT_DIR / run
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.json").write_text(json.dumps(
        {"model": args.model, "params": params, "seed": args.seed, "n_perts": len(perts),
         "panel_id": bundle.manifest["panel_id"]}, indent=1))

    t0 = time.time()

    def blocks():
        for name in bundle.contexts:
            print(f"[{run}] predict context {name} ({time.time() - t0:.0f}s)")
            ctx = bundle.load(name)
            for pert, X in model.predict_context(ctx, perts, bundle.cells_per_pert, args.seed):
                yield name, pert, X
            del ctx
            gc.collect()

    h5ad = run_dir / "prediction.h5ad"
    write_prediction(h5ad, bundle, blocks())
    print(f"[{run}] wrote {h5ad} ({h5ad.stat().st_size / 1e9:.1f} GB, {time.time() - t0:.0f}s)")

    if args.limit_perts or args.no_package:
        return
    vcc = run_dir / "prediction.vcc"
    vcc_prep(h5ad, bundle, vcc)
    print(f"[{run}] packaged {vcc} ({vcc.stat().st_size / 1e9:.1f} GB, {time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
