#!/usr/bin/env python3
"""Rebuild the compact streams and BED files in data/ from the downloaded sources, and check every hash."""

import argparse, gzip, hashlib, json, lzma, struct, shutil, time
from pathlib import Path


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        while b := f.read(8 << 20):
            h.update(b)
    return h.hexdigest()


def curate_chrom(src, out, chrom, r):
    count = seen = skip = 0
    previous = 0
    first = None
    maxend = 0
    totalspan = 0
    spansha = hashlib.sha256()
    namesha = hashlib.sha256()
    bedsha = hashlib.sha256()
    started = time.perf_counter()
    paths = [out / (chrom + "." + c + ".xz") for c in ["delta32", "span32", "meta16"]]
    compressors = [lzma.LZMACompressor(preset=6) for _ in paths]
    bufs = [bytearray(), bytearray(), bytearray()]
    with gzip.open(src / (chrom + ".bgzf"), "rb") as stream, (out / (chrom + ".bed")).open("wb") as bed:
        files = [p.open("wb") for p in paths]
        stream.read(r["first_virtual"] & 65535)
        while True:
            n = struct.unpack("<i", stream.read(4))[0]
            data = stream.read(n)
            assert len(data) == n
            ref, pos, bmq, fc, lseq = struct.unpack_from("<iiIIi", data)
            if ref != r["refid"]:
                break
            seen += 1
            flag = fc >> 16
            nops = fc & 65535
            namelen = bmq & 255
            mapq = (bmq >> 8) & 255
            if flag & 4 or pos < 0:
                skip += 1
                continue
            cigar = struct.unpack_from("<" + "I" * nops, data, 32 + namelen)
            span = sum(v >> 4 for v in cigar if (v & 15) in (0, 2, 3, 7, 8))
            if span <= 0:
                skip += 1
                continue
            assert pos >= previous, (chrom, pos, previous)
            strand = "-" if flag & 16 else "+"
            row = f"{chrom}\t{pos}\t{pos+span}\tHG002_{chrom}_{count+1}\t{mapq}\t{strand}\n".encode()
            bed.write(row)
            bedsha.update(row)
            bufs[0] += struct.pack("<I", pos - previous)
            bufs[1] += struct.pack("<I", span)
            bufs[2] += struct.pack("<H", mapq | ((flag & 16) << 4))
            spansha.update(struct.pack("<iiHB", pos, span, mapq, bool(flag & 16)))
            namesha.update(data[32 : 32 + namelen])
            count += 1
            previous = pos
            first = pos if first is None else first
            maxend = max(maxend, pos + span)
            totalspan += span
            if count % 100000 == 0:
                for fp, co, buf in zip(files, compressors, bufs):
                    fp.write(co.compress(buf))
                    buf.clear()
            if count % 2000000 == 0:
                print(chrom, count, round(time.perf_counter() - started, 1), flush=True)
        for fp, co, buf in zip(files, compressors, bufs):
            fp.write(co.compress(buf) + co.flush())
            fp.close()
    result = dict(
        chromosome=chrom,
        records=count,
        source_records_seen=seen,
        excluded_unmapped_or_empty=skip,
        first_start=first,
        max_end=maxend,
        total_reference_span=totalspan,
        span_depth_over_chromosome=totalspan / r["length"],
        bed_sha256=bedsha.hexdigest(),
        bed_bytes=(out / (chrom + ".bed")).stat().st_size,
        source_names_sha256=namesha.hexdigest(),
        interval_fields_sha256=spansha.hexdigest(),
        compact={p.name: {"sha256": sha(p), "bytes": p.stat().st_size} for p in paths},
    )
    (out / (chrom + ".json")).write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)
    return result


def annotate(src, out):
    transcripts = []
    exons = []
    accessions = set()
    chroms = {"chr" + str(i) for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    with gzip.open(src / "hg38_refGene.txt.gz", "rt") as f:
        for line in f:
            x = line.rstrip().split("\t")
            name, chrom, strand = x[1:4]
            if chrom not in chroms:
                continue
            accessions.add(name)
            transcripts.append((chrom, int(x[4]), int(x[5]), name, 0, strand))
            for i, (a, b) in enumerate(zip(x[9].rstrip(",").split(","), x[10].rstrip(",").split(",")), 1):
                exons.append((chrom, int(a), int(b), name + "_exon" + str(i), 0, strand))
    result = {}
    for name, rows in [("refseq_transcripts", transcripts), ("refseq_exons", exons)]:
        rows.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        raw = "".join("\t".join(map(str, x)) + "\n" for x in rows).encode()
        (out / (name + ".bed")).write_bytes(raw)
        (out / (name + ".bed.xz")).write_bytes(lzma.compress(raw, preset=9))
        result[name] = {
            "records": len(rows),
            "sum_bp": sum(x[2] - x[1] for x in rows),
            "bed_sha256": hashlib.sha256(raw).hexdigest(),
            "bed_bytes": len(raw),
            "compressed_bytes": (out / (name + ".bed.xz")).stat().st_size,
        }
    (out / "refseq_accessions.txt.xz").write_bytes(lzma.compress(("\n".join(sorted(accessions)) + "\n").encode()))
    if (src / "TwistCore.bed").exists():
        rows = [
            line.strip().split("\t")
            for line in (src / "TwistCore.bed").read_text().splitlines()
            if line.strip() and not line.startswith(("#", "track"))
        ]
        rows.sort(key=lambda r: (r[0], int(r[1]), int(r[2])))
        raw = ("\n".join("\t".join(r) for r in rows) + "\n").encode()
        (out / "exome_targets.bed").write_bytes(raw)
        (out / "exome_targets.bed.xz").write_bytes(lzma.compress(raw, preset=9))
        result["exome_targets"] = {
            "records": len(rows),
            "sum_bp": sum(int(x[2]) - int(x[1]) for x in rows),
            "bed_sha256": hashlib.sha256(raw).hexdigest(),
            "bed_bytes": len(raw),
            "compressed_bytes": (out / "exome_targets.bed.xz").stat().st_size,
        }
    (out / "annotations.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    data_dir = Path(__file__).resolve().parent.parent  # data/
    out = args.output_dir.resolve()
    if out == data_dir or data_dir in out.parents:
        raise ValueError("Use an empty scratch output directory outside data/")
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be empty")
    out.mkdir(parents=True, exist_ok=True)
    meta = json.loads((data_dir / "curation.json").read_text())
    for key in ["dev_reads", "holdout_reads"]:
        expected = meta["datasets"][key]
        c = expected["chromosome"]
        actual = curate_chrom(args.source_dir, out, c, meta["source"]["ranges"][c])
        for field in ["records", "bed_sha256", "bed_bytes", "interval_fields_sha256"]:
            if actual[field] != expected[field]:
                raise ValueError("Recuration mismatch " + c + " " + field)
        if actual["compact"] != expected["compact"]:
            raise ValueError("Compact source differs " + c)
    annotate(args.source_dir, out)
    for key, row in json.loads((out / "annotations.json").read_text()).items():
        if row["bed_sha256"] != meta["datasets"][key]["bed_sha256"]:
            raise ValueError("Annotation mismatch " + key)
    print("All exact source-to-BED and compact checks passed", flush=True)
