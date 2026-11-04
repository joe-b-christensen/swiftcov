#!/bin/bash
# Builds solver/swiftcov: the bedtools v2.31.1 tree in solver/bedtools with the SwiftCov fast path.
# bedtools objects use the same release flags as the reference (prepare.py: -O2 -g0, local zlib, no bzip2/lzma CRAM
# codecs); swiftcov_fast.cpp uses -O3 (src/swiftcov/Makefile.frag). SWIFTCOV_BUILD_JOBS=N sets make -j (default 2).
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$(pwd)"
DEPS="$ROOT/ext/deps"
if [ ! -f "$DEPS/lib/libz.a" ]; then
  python3 -c "import sys; sys.path.insert(0, '$ROOT'); from prepare import deps; deps('$ROOT/ext/deps')"
fi
SRC="$ROOT/solver/bedtools"
CONFIG="$SRC/src/utils/htslib/config.h"
WANT=$'#define _GNU_SOURCE 1\n#define HAVE_DRAND48 1\n#define HAVE_FSEEKO 1\n#define HAVE_FSYNC 1\n#define HAVE_FTELLO 1\n#define HAVE_GETPAGESIZE 1\n#define HAVE_MMAP 1'
if [ ! -f "$CONFIG" ] || [ "$(cat "$CONFIG")" != "$WANT" ]; then printf '%s\n' "$WANT" > "$CONFIG"; fi
# htslib's Makefile runs ./version.sh; the bit can be lost when a tree is copied or committed from Windows.
chmod +x "$SRC/src/utils/htslib/version.sh"
JOBS="${SWIFTCOV_BUILD_JOBS:-2}"
# Statically linked (start-up about 1 ms faster than with a dynamic libstdc++); the link is redone every build so a
# flag change cannot leave a stale binary.
rm -f "$SRC/bin/bedtools"
CPATH="$DEPS/include" LIBRARY_PATH="$DEPS/lib" make -C "$SRC" -j"$JOBS" bin/bedtools \
  CXXFLAGS="-O2 -g0" CFLAGS="-O2 -g0" LDFLAGS="-static" BT_LIBS="-lz -lm -lpthread" > "$ROOT/solver/.build.log" 2>&1 || { tail -50 "$ROOT/solver/.build.log"; exit 1; }
cp "$SRC/bin/bedtools" "$ROOT/solver/swiftcov.tmp" && mv -f "$ROOT/solver/swiftcov.tmp" "$ROOT/solver/swiftcov"
