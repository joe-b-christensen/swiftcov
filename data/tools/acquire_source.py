#!/usr/bin/env python3
"""Download the exact public source ranges behind data/ (about 1.9 GB) and check their hashes."""

import argparse, hashlib, json, struct, urllib.request, zlib
from pathlib import Path


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        while b := f.read(8 << 20):
            h.update(b)
    return h.hexdigest()


def fetch(url, path, expected, headers=None):
    if path.exists():
        if sha(path) != expected:
            raise ValueError("Existing source hash mismatch: " + str(path))
        return
    partial = path.with_suffix(path.suffix + ".partial")
    with urllib.request.urlopen(
        urllib.request.Request(url, headers=headers or {}), timeout=180
    ) as response, partial.open("wb") as out:
        if headers and response.status != 206:
            raise ValueError("Range request not honored")
        while chunk := response.read(4 << 20):
            out.write(chunk)
    if sha(partial) != expected:
        raise ValueError("Downloaded source hash mismatch: " + str(path))
    partial.replace(path)


def bigbed3_to_bed(source, destination):
    """Decode the downloaded BED3 BigBed via its B+ chromosome tree and zlib blocks."""
    raw = source.read_bytes()
    header = struct.unpack_from("<IHHQQQHHQQIQ", raw)
    assert header[0] == 0x8789F2EB and header[6] == 3 and header[7] == 3
    tree, data, index = header[3:6]
    magic, block, keysize, valsize, n, _ = struct.unpack_from("<IIIIQQ", raw, tree)
    assert magic == 0x78CA8C91 and valsize == 8
    chromosomes = {}

    def node(offset):
        leaf, _, count = struct.unpack_from("<BBH", raw, offset)
        offset += 4
        for i in range(count):
            key = raw[offset : offset + keysize].rstrip(b"\0").decode()
            offset += keysize
            if leaf:
                ident, length = struct.unpack_from("<II", raw, offset)
                chromosomes[ident] = key
            else:
                node(struct.unpack_from("<Q", raw, offset)[0])
            offset += 8

    node(tree + 32)
    expected = struct.unpack_from("<Q", raw, data)[0]
    pending = raw[data + 8 : index]
    rows = []
    while pending:
        decoder = zlib.decompressobj()
        block = decoder.decompress(pending)
        assert decoder.eof
        pending = decoder.unused_data
        offset = 0
        while offset < len(block):
            chrom, start, end = struct.unpack_from("<III", block, offset)
            offset += 12
            extra_end = block.index(b"\0", offset)
            assert extra_end == offset
            offset = extra_end + 1
            rows.append((chromosomes[chrom], start, end))
    assert len(rows) == expected
    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    destination.write_text("".join(f"{c}\t{a}\t{b}\n" for c, a, b in rows), newline="\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent  # data/
    meta = json.loads((root / "curation.json").read_text())
    out = args.source_dir.resolve()
    if out == root or root in out.parents:
        raise ValueError("Use a scratch directory outside data/")
    out.mkdir(parents=True, exist_ok=True)
    source = meta["source"]
    generation = source["source_headers"]["x-goog-generation"]
    for chrom, row in source["ranges"].items():
        fetch(
            source["bam_url"] + "?generation=" + generation,
            out / (chrom + ".bgzf"),
            row["download_sha256"],
            {"Range": f"bytes={row['range_start']}-{row['range_end']}"},
        )
        print("verified " + chrom, flush=True)
    for name in ["hg38_refGene.txt.gz", "TwistCore.bb"]:
        item = source["annotations"][name]
        fetch(item["url"], out / name, item["sha256"])
    bigbed3_to_bed(out / "TwistCore.bb", out / "TwistCore.bed")
    if sha(out / "TwistCore.bed") != meta["datasets"]["exome_targets"]["bed_sha256"]:
        raise ValueError("Exome target conversion differs")
    (out / "acquisition.json").write_text(json.dumps(source, indent=2))
    print("Exact sources verified; next run curate_source.py --source-dir HERE --output-dir EMPTY_SCRATCH", flush=True)


if __name__ == "__main__":
    main()
