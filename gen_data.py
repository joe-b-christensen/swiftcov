#!/usr/bin/env python3
"""SwiftCov: dataset specifications and test-case generators for bedtools coverage.

Two jobs:
  * write_dataset(root, name, spec, out): the benchmark BED files (dev + holdout). Real sets are rebuilt from the
    compact streams in data/ and checked against data/curation.json; synthetic sets come from fixed seeds.
  * fuzz_case(rng, k): 72 directed cases, then tiny random BED files, with bedtools coverage options, for the
    random-input check in evaluate.py.
"""
import gzip
import hashlib
import json
import lzma
import os
import random
import struct
import zlib
from pathlib import Path

DEV_SPECS = {
    "real_transcript_default": dict(real=True, a="refseq_transcripts", b="dev_reads", args=[]),
    "real_exome_depth": dict(real=True, a="exome_targets", b="dev_reads", args=["-d"]),
    "sim_dense_default": dict(real=False, seed=710, n=2000, m=30000, length=40000, args=[]),
}
HOLDOUT_SPECS = {
    "hold_real_transcript_mean": dict(real=True, a="refseq_transcripts", b="holdout_reads", args=["-mean"]),
    "hold_real_transcript_hist_sorted": dict(
        real=True, a="refseq_transcripts", b="holdout_reads", args=["-hist", "-sorted"]
    ),
    "hold_real_exon_hist": dict(real=True, a="refseq_exons", b="holdout_reads", args=["-hist"]),
    "hold_sim_fraction_strand": dict(real=False, seed=991, n=1600, m=22000, length=15000, args=["-f", "0.05", "-s"]),
}


def materialize_real(root, key, out):
    """Stream exact BED reconstruction once; benchmark aliases share its inode."""
    target = out / ("source_" + key + ".bed")
    if target.exists():
        return target
    record = json.loads((root / "data/curation.json").read_text())["datasets"][key]
    digest = hashlib.sha256()
    rows = 0
    partial = target.with_suffix(".partial")
    with partial.open("wb") as dst:
        if record["format"] == "bed.xz":
            with lzma.open(root / "data" / record["file"], "rb") as src:
                while block := src.read(1 << 20):
                    dst.write(block)
                    digest.update(block)
                    rows += block.count(b"\n")
        else:
            chrom = record["chromosome"]
            position = 0
            with lzma.open(root / "data" / record["delta"], "rb") as starts, lzma.open(
                root / "data" / record["span"], "rb"
            ) as spans, lzma.open(root / "data" / record["meta"], "rb") as metas:
                while a := starts.read(4 * 65536):
                    b = spans.read(len(a))
                    c = metas.read(len(a) // 2)
                    if len(a) % 4 or len(b) != len(a) or len(c) * 2 != len(a):
                        raise ValueError("Truncated compact fields: " + key)
                    block = bytearray()
                    for (delta,), (span,), (meta,) in zip(
                        struct.iter_unpack("<I", a), struct.iter_unpack("<I", b), struct.iter_unpack("<H", c)
                    ):
                        position += delta
                        rows += 1
                        strand = "-" if meta & 256 else "+"
                        block += f"{chrom}\t{position}\t{position+span}\tHG002_{chrom}_{rows}\t{meta&255}\t{strand}\n".encode()
                    dst.write(block)
                    digest.update(block)
                if spans.read(1) or metas.read(1):
                    raise ValueError("Unequal compact columns: " + key)
    if (
        rows != record["records"]
        or digest.hexdigest() != record["bed_sha256"]
        or partial.stat().st_size != record["bed_bytes"]
    ):
        partial.unlink()
        raise ValueError("Real-data reconstruction SHA/count mismatch: " + key)
    partial.replace(target)
    return target


def write_dataset(root, name, spec, out):
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    if spec["real"]:
        for k in ("a", "b"):
            source = materialize_real(root, spec[k], out)
            p = out / (name + "_" + k + ".bed")
            if not p.exists():
                os.link(source, p)
            paths[k] = str(p.relative_to(root))
    else:
        rng = random.Random(spec["seed"])
        for k, count, length in [("a", spec["n"], spec["length"]), ("b", spec["m"], spec["length"] // 2)]:
            rows = []
            for i in range(count):
                s = rng.randrange(0, 1000000)
                e = s + rng.randrange(max(1, length // 2), length + 1)
                rows.append(("chr1", s, e, f"{k}{i}", 0, rng.choice("+-")))
            rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
            p = out / (name + "_" + k + ".bed")
            p.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
            paths[k] = str(p.relative_to(root))
    return paths


def bed(rows, nl=b"\n"):
    return nl.join("\t".join(map(str, r)).encode() for r in rows) + (nl if rows else b"")


def tiny_bam(rows, cigars=None):
    """Valid uncompressed BAM payload inside standard BGZF; tiny fixtures only."""
    refs = ["chr1", "chr2", "chr10"]
    header = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{c}\tLN:10000\n" for c in refs)
    raw = b"BAM\1" + struct.pack("<i", len(header)) + header.encode() + struct.pack("<i", len(refs))
    for c in refs:
        raw += struct.pack("<i", len(c) + 1) + c.encode() + b"\0" + struct.pack("<i", 10000)
    for index, (c, s, e, n, score, strand) in enumerate(rows):
        ops = cigars[index] if cigars is not None else [(max(1, e - s), 0)]
        length = sum(n for n, op in ops if op in (0, 1, 4, 7, 8))
        name = n.encode() + b"\0"
        flag = 16 if strand == "-" else 0
        core = struct.pack(
            "<iiIIiiii",
            refs.index(c),
            s,
            (4681 << 16) | (60 << 8) | len(name),
            (flag << 16) | len(ops),
            length,
            -1,
            -1,
            0,
        )
        data = (
            core
            + name
            + b"".join(struct.pack("<I", (n << 4) | op) for n, op in ops)
            + b"\x11" * ((length + 1) // 2)
            + b"\xff" * length
        )
        raw += struct.pack("<i", len(data)) + data

    def block(data):
        co = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = co.compress(data) + co.flush()
        return (
            b"\x1f\x8b\x08\x04"
            + b"\0" * 4
            + b"\0\xff"
            + struct.pack("<H", 6)
            + b"BC"
            + struct.pack("<HH", 2, len(payload) + 25)
            + payload
            + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
        )

    return block(raw) + block(b"")


def directed_cases():
    a = b"chr1\t0\t10\ta\t0\t+\nchr1\t10\t20\tdup\t0\t-\nchr2\t5\t8\tdup\t0\t+\n"
    b = b"chr1\t2\t7\tb\t0\t+\nchr1\t5\t12\tb\t0\t-\nchr2\t7\t9\tc\t0\t+\n"
    base = {"a": a, "b": b, "g": b"chr1\t1000\nchr2\t1000\nchr10\t1000\n"}
    cases = []

    def add(tag, flags=(), files=None):
        cases.append((files or base.copy(), ["-a", "{a}", "-b", "{b}", *flags], tag))

    for mode in ([], ["-counts"], ["-hist"], ["-d"], ["-mean"]):
        add("mode_" + ("default" if not mode else mode[0]), mode)
        add("sorted_" + ("default" if not mode else mode[0]), [*mode, "-sorted", "-g", "{g}"])
    for flags in (
        ["-s"],
        ["-S"],
        ["-f", "0.5"],
        ["-F", "0.5"],
        ["-f", "0.5", "-r"],
        ["-f", "0.7", "-F", "0.1", "-e"],
        ["-nonamecheck"],
        ["-nobuf"],
        ["-iobuf", "64K"],
        ["-bed"],
        ["-ubam"],
        ["-fbam"],
        ["-sortout"],
        ["-prec", "3"],
        ["-seed", "17"],
    ):
        add("flag_" + "_".join(flags), flags)
    for flag in ["-wa", "-wb", "-wo", "-wao", "-loj", "-u", "-v", "-c", "-C", "-filenames"]:
        add("inherited_" + flag, [flag])
    add("names", ["-names", "database"])
    add("header", ["-header"], {**base, "a": b"# first\ntrack name=x\n" + a})
    add("multi_b", ["-b", "{b2}", "-hist"], {**base, "b2": b})
    add("duplicate_intervals", ["-hist"], {**base, "b": b + b})
    add("crlf", [], {**base, "a": a.replace(b"\n", b"\r\n"), "b": b.replace(b"\n", b"\r\n")})
    add("no_final_newline", [], {**base, "a": a.rstrip(), "b": b.rstrip()})
    add("zero_length", ["-hist"], {**base, "a": b"chr1\t5\t5\tx\t0\t+\n", "b": b"chr1\t5\t5\tb\t0\t+\n"})
    add("touching", [], {**base, "a": b"chr1\t0\t1\n", "b": b"chr1\t1\t2\n"})
    add("empty_a", [], {**base, "a": b""})
    add("empty_b", ["-hist"], {**base, "b": b""})
    add("negative_error", [], {**base, "a": b"chr1\t-1\t3\n"})
    add("reversed_error", [], {**base, "b": b"chr1\t20\t3\n"})
    add("malformed_error", [], {**base, "b": b"chr1\tnonsense\t3\n"})
    add("missing_field", [], {**base, "a": b"chr1\t2\n"})
    add("mutually_exclusive", ["-counts", "-mean"])
    add("fraction_error", ["-f", "0"])
    add("strand_error", ["-s", "-S"])
    add("sorted_error", ["-sorted"], {**base, "a": b"chr1\t10\t20\nchr1\t0\t10\n"})
    add("gzip", [], {**base, "b": gzip.compress(b, mtime=0)})
    blocks = b"chr1\t0\t30\tspliced\t0\t+\t0\t30\t0\t2\t5,5,\t0,25,\n"
    add(
        "bed12_split",
        ["-split", "-hist"],
        {**base, "a": blocks, "b": b"chr1\t2\t28\tr\t0\t+\t2\t28\t0\t2\t4,4,\t0,22,\n"},
    )
    add("gff", ["-mean"], {**base, "a": b'chr1\tRefSeq\texon\t1\t10\t.\t+\t.\tgene_id "x";\n'})
    vcf = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t3\t.\tA\tT\t.\tPASS\t.\nchr1\t6\tdel\tA\t<DEL>\t.\tPASS\tEND=12\n"
    add("vcf", ["-hist"], {**base, "a": vcf})
    rows = [("chr1", 2, 7, "read", 0, "+"), ("chr1", 5, 12, "read", 0, "-")]
    add("bam_b", ["-mean"], {**base, "b": tiny_bam(rows)})
    add("bam_a", ["-counts"], {**base, "a": tiny_bam(rows)})
    spliced = tiny_bam([("chr1", 2, 32, "splice", 0, "+")], [[(5, 0), (10, 3), (2, 1), (5, 0), (5, 2), (5, 0)]])
    add("bam_split_cigar_hist", ["-split", "-hist"], {**base, "b": spliced})
    add("bam_split_cigar_depth", ["-split", "-d"], {**base, "a": spliced})
    cases.append(({**base, "a": tiny_bam(rows)}, ["-abam", "{a}", "-b", "{b}", "-mean"], "legacy_abam"))
    for flags in [
        ["-i", "{a}"],
        ["-n", "2"],
        ["-o", "sum"],
        ["-ops", "sum"],
        ["-opCols", "4"],
        ["-null", "NA"],
        ["-delim", ","],
        ["-h"],
        ["--help"],
    ]:
        add("parser_" + "_".join(flags), flags)
    add(
        "large_coordinate",
        ["-sorted"],
        {**base, "a": b"chr1\t3000000000\t3000000002\n", "b": b"chr1\t3000000001\t3000000003\n"},
    )
    return cases


def fuzz_case(rng, k=None):
    directed = directed_cases()
    if k is not None and k < len(directed):
        return directed[k]
    chroms = ["chr1", "chr2", "chr10"]
    aa = []
    bb = []
    for rows, n in [(aa, rng.randrange(1, 15)), (bb, rng.randrange(0, 35))]:
        for i in range(n):
            s = rng.randrange(0, 200)
            e = s + rng.choice([0, 1, 2, 5, 10, 30, 70])
            rows.append((rng.choice(chroms), s, e, rng.choice(["dup", f"x{i}"]), 0, rng.choice(["+", "-", "."])))
        rows.sort(key=lambda r: (chroms.index(r[0]), r[1], r[2]))
    files = {"a": bed(aa), "b": bed(bb), "g": b"chr1\t1000\nchr2\t1000\nchr10\t1000\n"}
    flags = rng.choice([[], ["-counts"], ["-mean"], ["-hist"], ["-d"]]).copy()
    if rng.random() < 0.35:
        flags += ["-sorted", "-g", "{g}"]
    if rng.random() < 0.4:
        flags += [rng.choice(["-s", "-S"])]
    if rng.random() < 0.45:
        flags += ["-f", rng.choice(["0.000001", "0.1", "0.5", "0.99", "1"])]
        flags += rng.choice([[], ["-r"], ["-F", "0.3"], ["-F", "0.8", "-e"]])
    if rng.random() < 0.12:
        files["b2"] = files["b"]
        flags += ["-b", "{b2}"]
    if rng.random() < 0.08:
        files["b"] = gzip.compress(files["b"], mtime=0)
    if rng.random() < 0.08:
        files["a"] = b"# a header\n" + files["a"]
        flags += ["-header"]
    if rng.random() < 0.08:
        files["a"] = files["a"].replace(b"\n", b"\r\n")
    return files, ["-a", "{a}", "-b", "{b}", *flags], "random"
