#!/usr/bin/env python3
"""SwiftCov: correctness check against bedtools v2.31.1 + benchmark.

    python3 evaluate.py                    # build, random-input check, benchmark-set check + timing
    python3 evaluate.py --holdout          # the same on the four holdout sets
    python3 evaluate.py --quick            # build + random-input check only
    python3 evaluate.py --threads 1        # timing with OMP_NUM_THREADS=1

The check has two parts, and both must pass:
  * Fuzz: 300 cases, every run: 72 directed cases covering every output mode and option class, malformed and
    zero-length records and every non-BED input format, then random BED files from a fresh seed (printed). Input
    files get random names, each process runs in its own directory, and SwiftCov's thread count is drawn from
    {1, 2, 3, 4, 7}. The reference runs on the same files. Standard output and any created files must match byte
    for byte, and both programs must exit 0 or both non-zero. Stderr wording is ignored.
  * Benchmark sets: SwiftCov's output must match the recorded sha256 and byte count of the reference output.
Timing: wall clock of the whole process with output written to a fresh file, deleted before each repetition
outside the timed window; best of 3 below 5 s. The reference times come from bench/manifest.json.
"""
import argparse
import json
import math
import os
import pathlib
import random
import statistics
import subprocess
import tempfile
import time

import gen_data
from prepare import sha

ROOT = pathlib.Path(__file__).resolve().parent
REF = ROOT / "ext/bedtools/bin/bedtools"
CAND = ROOT / "solver/swiftcov"


def run(binary, flags, cwd, threads, timeout=20):
    try:
        r = subprocess.run(
            [str(binary), "coverage", *flags],
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            env={**os.environ, "OMP_NUM_THREADS": str(threads), "LC_ALL": "C"},
        )
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return "timeout", b"", b""


def compare_case(binary, case, seed, index):
    files, flags, tag = case
    rng = random.Random(seed + index)
    names = {k: f"in_{rng.randrange(2**40):010x}" for k in files}
    actual = [names.get(a[1:-1], a) if a.startswith("{") else a for a in flags]
    with tempfile.TemporaryDirectory(prefix="swiftcov_fuzz_") as d:
        root = pathlib.Path(d)
        rd = root / "r"
        cd = root / "c"
        rd.mkdir()
        cd.mkdir()
        for p in [rd, cd]:
            for k, data in files.items():
                (p / names[k]).write_bytes(data)
        threads = rng.choice([1, 2, 3, 4, 7])
        rr, ro, re = run(REF, actual, rd, 1)
        cr, co, ce = run(binary, actual, cd, threads)
        rf = {p.name: p.read_bytes() for p in rd.iterdir() if p.is_file() and p.name not in names.values()}
        cf = {p.name: p.read_bytes() for p in cd.iterdir() if p.is_file() and p.name not in names.values()}
        # A reference timeout fails the case; it is never mistaken for the reference rejecting the input.
        ok = (
            rr != "timeout"
            and cr != "timeout"
            and ro == co
            and rf == cf
            and ((rr == 0 and cr == 0) or (isinstance(rr, int) and rr != 0 and isinstance(cr, int) and cr != 0))
        )
        detail = dict(
            index=index,
            tag=tag,
            flags=actual,
            names=names,
            threads=threads,
            ref_rc=rr,
            cand_rc=cr,
            ref_stdout=ro[:1000].decode(errors="replace"),
            cand_stdout=co[:1000].decode(errors="replace"),
            ref_stderr=re[:700].decode(errors="replace"),
            cand_stderr=ce[:700].decode(errors="replace"),
        )
        if not ok:
            fail = ROOT / "failures" / f"{seed}_{index}"
            fail.mkdir(parents=True, exist_ok=True)
            for k, data in files.items():
                (fail / names[k]).write_bytes(data)
            (fail / "case.json").write_text(json.dumps(detail, indent=2))
        return ok, detail


def fuzz_gate(n, seed, binary=CAND):
    rng = random.Random(seed)
    fails = []
    passed = 0
    tags = {}
    for k in range(n):
        case = gen_data.fuzz_case(rng, k)
        tags[case[2]] = tags.get(case[2], 0) + 1
        ok, detail = compare_case(binary, case, seed, k)
        if ok:
            passed += 1
        else:
            fails.append(detail)
            print("FUZZ FAIL", json.dumps(detail), flush=True)
            if len(fails) >= 5:
                break
    return dict(seed=seed, requested=n, executed=passed + len(fails), passed=passed, failures=len(fails), tags=tags)


class Bench:
    def __init__(self, row, td):
        self.row = row
        self.td = pathlib.Path(td)
        self.k = 0
        self.flags = []
        for k in ("a", "b"):
            name = f"in_{random.SystemRandom().randrange(2**48):012x}"
            os.symlink(ROOT / row["paths"][k], self.td / name)
            self.flags += ["-" + k, name]
        self.flags += row["args"]
        self.timeout = max(60, 5 * row["ref_j1_s"])

    def once(self, binary, threads):
        self.k += 1
        output = self.td / f"stdout_{self.k}"
        output.unlink(missing_ok=True)  # outside the timer, before every repetition
        env = {**os.environ, "OMP_NUM_THREADS": str(threads), "LC_ALL": "C"}
        with output.open("wb") as f:
            start = time.perf_counter()
            try:
                r = subprocess.run(
                    [str(binary), "coverage", *self.flags],
                    cwd=self.td,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    timeout=self.timeout,
                    env=env,
                )
                rc = r.returncode
            except subprocess.TimeoutExpired:
                rc = "timeout"
            elapsed = time.perf_counter() - start
        ok = (
            rc == 0
            and output.stat().st_size == self.row["expected_bytes"]
            and sha(output) == self.row["expected_sha256"]
        )
        if not ok:
            print("BENCH FAIL", self.row["name"], str(binary), rc, flush=True)
        output.unlink(missing_ok=True)
        return elapsed, ok


def benchmark(rows, threads):
    result = []
    fails = 0
    for row in rows:
        with tempfile.TemporaryDirectory(prefix="swiftcov_bench_") as td:
            b = Bench(row, td)
            times = []
            ok = True
            for _ in range(3):
                dt, ok = b.once(CAND, threads)
                if not ok:
                    break
                times.append(dt)
                if dt >= 5:
                    break
            if not ok:
                fails += 1
                result.append(dict(name=row["name"], ok=False))
                continue
            best = min(times)
            result.append(
                dict(
                    name=row["name"],
                    real=row["real"],
                    ok=True,
                    seconds=best,
                    repetitions_s=times,
                    speedup=row["ref_j1_s"] / best,
                    speedup_vs_upstream_best=row["ref_best_s"] / best,
                )
            )
    return result, fails


def build_solver(root):
    binary = root / "solver/swiftcov"
    binary.unlink(missing_ok=True)
    try:
        r = subprocess.run(
            ["bash", str(root / "solver/build.sh")], cwd=root, capture_output=True, text=True, timeout=900
        )
    except subprocess.TimeoutExpired:
        return False
    if r.returncode or not binary.exists():
        print(r.stdout[-3000:], r.stderr[-3000:])
        return False
    return True


def gm(values):
    return math.exp(statistics.mean(math.log(x) for x in values))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--holdout", action="store_true")
    ap.add_argument("--threads", type=int)
    ap.add_argument("--fuzz", type=int, default=300)
    args = ap.parse_args()
    if args.fuzz < len(gen_data.directed_cases()):
        ap.error("--fuzz must cover all directed cases")
    if args.threads is not None and args.threads < 1:
        ap.error("--threads must be positive")
    manifest = json.loads((ROOT / "bench/manifest.json").read_text())
    rows = manifest["holdout" if args.holdout else "dev"]
    bad = [p for p, h in manifest["scripts_sha256"].items() if not (ROOT / p).exists() or sha(ROOT / p) != h]
    if sha(REF) != manifest["reference"]["sha256"]:
        bad.append("reference binary")
    for row in rows:
        for k, p in row["paths"].items():
            if sha(ROOT / p) != row["input_sha256"][k]:
                bad.append(p)
    if bad:
        print("MODIFIED since prepare.py:", bad, "Restore them, or delete bench/ and rerun prepare.py.")
        raise SystemExit(2)
    start = time.perf_counter()
    if not build_solver(ROOT):
        print("gate: BUILD_FAIL")
        raise SystemExit(1)
    threads = args.threads or manifest["threads"]
    seed = random.SystemRandom().randrange(10**9)
    fuzz = fuzz_gate(args.fuzz, seed)
    mode = "quick" if args.quick else "holdout" if args.holdout else "dev"
    summary = dict(mode=mode, threads=threads, fuzz=fuzz, gate="FAIL")
    failures = 0
    if not fuzz["failures"] and not args.quick:
        result, failures = benchmark(rows, threads)
        if not failures:
            summary.update(
                geomean_s=gm([r["seconds"] for r in result]),
                speedup_vs_ref=gm([r["speedup"] for r in result]),
                speedup_vs_upstream_best=gm([r["speedup_vs_upstream_best"] for r in result]),
            )
            real = [r for r in result if r["real"]]
            if real:
                summary.update(
                    real_speedup_vs_ref=gm([r["speedup"] for r in real]),
                    real_speedup_vs_upstream_best=gm([r["speedup_vs_upstream_best"] for r in real]),
                )
        summary["per_dataset"] = result
    summary.update(
        gate="PASS" if not fuzz["failures"] and not failures else "FAIL",
        bench_failures=failures,
        eval_seconds=time.perf_counter() - start,
    )
    (ROOT / ("last_eval_" + mode + ".json")).write_text(json.dumps(summary, indent=2))
    print("gate: " + summary["gate"])
    print(f"fuzz: {fuzz['passed']}/{fuzz['requested']} PASS seed={seed}")
    for key in [
        "geomean_s",
        "speedup_vs_ref",
        "speedup_vs_upstream_best",
        "real_speedup_vs_ref",
        "real_speedup_vs_upstream_best",
        "threads",
        "eval_seconds",
    ]:
        if key in summary:
            print(f"{key}: {summary[key]}")
    print("per_dataset: " + json.dumps(summary.get("per_dataset", [])))
    raise SystemExit(0 if summary["gate"] == "PASS" else 1)


if __name__ == "__main__":
    main()
