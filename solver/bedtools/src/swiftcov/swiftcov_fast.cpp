/*
 * swiftcov_fast.cpp -- SwiftCov fast path for `bedtools coverage`.
 *
 * Part of SwiftCov, a re-implementation of the hot path of bedtools coverage
 * v2.31.1 (Quinlan & Hall).  Copyright (c) 2026 Joe Christensen.  MIT licence,
 * as for bedtools v2.31.1; see LICENSE at the top of the repository.
 *
 * Design
 * ------
 * swiftcov_fast_coverage() is called from bedtools' main() before the
 * upstream driver.  It accepts a deliberately narrow, fully validated subset
 * of invocations (plain tab-delimited BED query/database files and a fixed
 * set of flags).  Everything is parsed and validated *before* any byte is
 * written.  If any argument, file or record falls outside the subset whose
 * upstream behaviour has been characterised, the function returns -1 and the
 * unmodified upstream implementation (compiled into the same binary) handles
 * the invocation, including all error paths.
 *
 * Within the subset, the upstream semantics reproduced here are:
 *   - zero-length records are widened to [start-1, end+1) for all overlap
 *     and depth computations (Record::adjustZeroLength) but printed verbatim;
 *   - a hit is any database record on the same chromosome name with
 *     min(end) > max(start) in adjusted coordinates, filtered by -s/-S and by
 *     float overlap fractions (Record::sameChromIntersects);
 *   - the per-A depth array is the sum of clipped hits (makeDepthCount);
 *   - output formats of doCounts/doDefault/doMean/doHist/doPerBase and the
 *     final "all" histogram, with float division and "%0.7f" formatting.
 * Under -sorted the upstream chromsweep equals this ideal hit set only when
 * every input is contiguous per chromosome, has non-decreasing original and
 * adjusted starts, and visits chromosomes in strictly increasing
 * lexicographic (or -g genome-file) order; other inputs fall back.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cerrno>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <emmintrin.h>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <memory>

namespace swiftcov {

typedef int64_t pos_t;
/* stored database coordinates (records needing more fall back to upstream) */
typedef int32_t rec_t;

static bool g_debug = false;
static bool g_prof = false;
#include <time.h>
static double nowSec() { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
static double g_t0 = 0, g_tl = 0;
static void prof(const char *what) {
	if (!g_prof) return;
	double t = nowSec();
	fprintf(stderr, "swiftcov: prof %-28s %8.3f ms (total %8.3f)\n", what, (t - g_tl) * 1e3, (t - g_t0) * 1e3);
	g_tl = t;
}

static bool fallback(const char *why) {
	if (g_debug) fprintf(stderr, "swiftcov: fallback: %s\n", why);
	return false;
}

/* ------------------------------------------------------------------ */
/* Configuration                                                       */
/* ------------------------------------------------------------------ */

enum Mode { M_DEFAULT, M_COUNTS, M_PERBASE, M_HIST, M_MEAN };

struct Config {
	std::string aPath;
	std::vector<std::string> bPaths;
	std::string gPath;
	Mode mode = M_DEFAULT;
	int nModes = 0;
	bool sorted = false, hasG = false;
	bool same = false, diff = false;
	bool haveFA = false, haveFB = false, recip = false, either = false;
	float fa = 0.0f, fb = 0.0f;
	int threads = 1;
};

static bool startsWithDash(const char *s) { return s[0] == '-'; }

static bool parseArgs(int argc, char **argv, Config &c) {
	int nA = 0;
	for (int i = 2; i < argc; i++) {
		const char *a = argv[i];
		if (!strcmp(a, "-a")) {
			if (i + 1 >= argc || startsWithDash(argv[i + 1])) return fallback("-a value");
			c.aPath = argv[++i]; nA++;
		} else if (!strcmp(a, "-b")) {
			if (i + 1 >= argc || startsWithDash(argv[i + 1])) return fallback("-b value");
			do {
				c.bPaths.push_back(argv[i + 1]);
				i++;
			} while (i + 1 < argc && argv[i + 1][0] != '-');
		} else if (!strcmp(a, "-counts")) { c.mode = M_COUNTS; c.nModes++; }
		else if (!strcmp(a, "-d")) { c.mode = M_PERBASE; c.nModes++; }
		else if (!strcmp(a, "-hist")) { c.mode = M_HIST; c.nModes++; }
		else if (!strcmp(a, "-mean")) { c.mode = M_MEAN; c.nModes++; }
		else if (!strcmp(a, "-sorted")) c.sorted = true;
		else if (!strcmp(a, "-g")) {
			if (i + 1 >= argc || startsWithDash(argv[i + 1])) return fallback("-g value");
			if (c.hasG) return fallback("-g twice");
			c.gPath = argv[++i]; c.hasG = true;
		}
		else if (!strcmp(a, "-s")) c.same = true;
		else if (!strcmp(a, "-S")) c.diff = true;
		else if (!strcmp(a, "-f")) {
			if (i + 1 >= argc || startsWithDash(argv[i + 1])) return fallback("-f value");
			c.haveFA = true; c.fa = (float)atof(argv[++i]);
		} else if (!strcmp(a, "-F")) {
			if (i + 1 >= argc || startsWithDash(argv[i + 1])) return fallback("-F value");
			c.haveFB = true; c.fb = (float)atof(argv[++i]);
		}
		else if (!strcmp(a, "-r")) c.recip = true;
		else if (!strcmp(a, "-e")) c.either = true;
		else if (!strcmp(a, "-header")) { /* no header lines are accepted by the fast path: no-op */ }
		else if (!strcmp(a, "-nonamecheck") || !strcmp(a, "-nobuf") || !strcmp(a, "-bed") ||
		         !strcmp(a, "-sortout") || !strcmp(a, "-ubam") || !strcmp(a, "-fbam")) { /* no output effect */ }
		else if (!strcmp(a, "-seed")) {
			if (i + 1 >= argc || startsWithDash(argv[i + 1])) return fallback("-seed value");
			i++;
		}
		else return fallback("unhandled argument");
	}
	if (nA != 1 || c.bPaths.empty()) return fallback("need exactly one -a and some -b");
	if (c.nModes > 1) return fallback("mutually exclusive modes");
	if (c.same && c.diff) return fallback("-s and -S");
	if (c.haveFA && (c.fa <= 0.0 || c.fa > 1.0)) return fallback("-f range");
	if (c.haveFB && (c.fb <= 0.0 || c.fb > 1.0)) return fallback("-F range");
	if (c.haveFB && !c.haveFA && c.either) return fallback("-e without -f");
	if (c.haveFB && c.recip) return fallback("-r with -F");
	if (c.either && !c.haveFA) return fallback("-e without -f");
	if (c.recip && !c.haveFA) return fallback("-r without -f");
	if (c.haveFA && c.recip) c.fb = c.fa;
	/* -f X -e without -F: zero-overlap (touching) records satisfy the B
	 * fraction 0 >= 0 test upstream, and which touching records the bin
	 * tree / sweep visits is implementation specific. */
	if (c.either && !c.haveFB && !c.recip) return fallback("-e with zero B fraction");
	if (c.hasG && !c.sorted) return fallback("-g without -sorted");
	if (c.aPath == "-" || c.aPath == "stdin") return fallback("stdin");
	for (auto &b : c.bPaths) if (b == "-" || b == "stdin") return fallback("stdin");
	if (c.hasG && (c.gPath == "-" || c.gPath == "stdin")) return fallback("stdin genome");
	const char *t = getenv("OMP_NUM_THREADS");
	if (t && *t) {
		int v = atoi(t);
		c.threads = v >= 1 ? (v > 64 ? 64 : v) : 1;
	}
	return true;
}

/* ------------------------------------------------------------------ */
/* Files                                                               */
/* ------------------------------------------------------------------ */

struct MappedFile {
	const char *data = nullptr;
	size_t size = 0;
	bool open(const std::string &path) {
		int fd = ::open(path.c_str(), O_RDONLY);
		if (fd < 0) return fallback("open failed");
		struct stat st;
		if (fstat(fd, &st) != 0 || !S_ISREG(st.st_mode)) { ::close(fd); return fallback("not a regular file"); }
		if (st.st_size <= 0) { ::close(fd); return fallback("empty file"); }
		size = (size_t)st.st_size;
		/* Reserve the file length rounded up to a page plus one zero guard page,
		 * then map the file over it: the parser may read up to 4 KiB past EOF and
		 * always sees zero bytes there. */
		size_t pg = 4096, mapLen = (size + pg - 1) & ~(pg - 1);
		void *base = mmap(nullptr, mapLen + pg, PROT_READ, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if (base == MAP_FAILED) { ::close(fd); return fallback("mmap reserve failed"); }
		void *m = mmap(base, size, PROT_READ, MAP_PRIVATE | MAP_FIXED | MAP_POPULATE, fd, 0);
		::close(fd);
		if (m == MAP_FAILED || m != base) return fallback("mmap failed");
		data = (const char *)m;
		return true;
	}
};

/* Map a regular file read-only (MAP_POPULATE) followed by a zero guard page, so
 * the scanner may read up to 4 KiB past EOF. */
static const char *mapWithGuard(int fd, size_t size) {
	/* The reservation is 2 MiB-aligned so that large page-cache folios can be
	 * mapped by PMD entries: measured for the 603 MB chr21 file, MAP_POPULATE
	 * takes ~1 ms aligned vs ~10 ms at a 4 KiB-aligned address, and a scan over
	 * the mapping is ~13% faster (fewer TLB misses). */
	const size_t pg = 4096, HP = (size_t)2 << 20;
	const size_t mapLen = (size + pg - 1) & ~(pg - 1);
	const size_t resLen = mapLen + pg + HP;
	char *res = (char *)mmap(nullptr, resLen, PROT_READ, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (res == MAP_FAILED) return nullptr;
	char *base = (char *)(((uintptr_t)res + HP - 1) & ~(uintptr_t)(HP - 1));
	if (base > res) munmap(res, (size_t)(base - res));
	char *keepEnd = base + mapLen + pg; /* file pages + one zero guard page */
	if (res + resLen > keepEnd) munmap(keepEnd, (size_t)(res + resLen - keepEnd));
	void *m = mmap(base, size, PROT_READ, MAP_PRIVATE | MAP_FIXED | MAP_POPULATE, fd, 0);
	if (m == MAP_FAILED || m != (void *)base) return nullptr;
	madvise(m, size, MADV_SEQUENTIAL);
	return (const char *)m;
}

/* Large zero-initialised arrays from anonymous memory, transparent huge pages
 * where available; capacity is an upper bound and untouched pages cost nothing. */
template <class T>
struct BigArray {
	T *p = nullptr;
	size_t n = 0, cap = 0;
	bool alloc(size_t c) {
		size_t bytes = c * sizeof(T);
		if (bytes < ((size_t)16 << 20)) {
			/* small inputs: plain heap memory, no huge-page zeroing */
			p = (T *)malloc(bytes + sizeof(T));
			cap = c; n = 0;
			return p != nullptr;
		}
		bytes = (bytes + (2u << 20) - 1) & ~(size_t)((2u << 20) - 1);
		if (bytes == 0) bytes = 2u << 20;
		void *m = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
		if (m == MAP_FAILED) return false;
#ifdef MADV_HUGEPAGE
		madvise(m, bytes, MADV_HUGEPAGE);
#endif
		p = (T *)m; cap = bytes / sizeof(T); n = 0;
		return true;
	}
};

/* ------------------------------------------------------------------ */
/* Chromosome names                                                    */
/* ------------------------------------------------------------------ */

struct ChromTable {
	std::unordered_map<std::string, int> ids;
	std::vector<std::string> names;
	int get(const char *p, size_t n) {
		std::string k(p, n);
		auto it = ids.find(k);
		if (it != ids.end()) return it->second;
		int id = (int)names.size();
		ids.emplace(k, id);
		names.push_back(k);
		return id;
	}
};

/* Genome file: strict subset of NewGenomeFile parsing. */
static bool loadGenome(const std::string &path, std::unordered_map<std::string, int> &gid) {
	MappedFile f;
	if (!f.open(path)) return false;
	const char *p = f.data, *end = f.data + f.size;
	int id = 0;
	while (p < end) {
		const char *nl = (const char *)memchr(p, '\n', end - p);
		const char *le = nl ? nl : end;
		if (le == p) return fallback("genome: empty line");
		const char *t = (const char *)memchr(p, '\t', le - p);
		if (!t || t == p) return fallback("genome: no tab");
		const char *q = t + 1;
		const char *t2 = (const char *)memchr(q, '\t', le - q);
		const char *fe = t2 ? t2 : le;
		if (fe == q || fe - q > 18) return fallback("genome: size field");
		if (*q == '0') return fallback("genome: size zero/leading zero");
		for (const char *r = q; r < fe; r++) if (*r < '0' || *r > '9') return fallback("genome: size digits");
		for (const char *r = p; r < le; r++) if (*r == '\r' || *r == ' ') return fallback("genome: CR/space");
		std::string name(p, t - p);
		if (gid.count(name)) return fallback("genome: duplicate chrom");
		gid.emplace(name, id++);
		p = le + 1;
	}
	if (id == 0) return fallback("genome: no entries");
	return true;
}

/* ------------------------------------------------------------------ */
/* BED parsing                                                         */
/* ------------------------------------------------------------------ */

enum { STR_UNKNOWN = 0, STR_FWD = 1, STR_REV = 2 };

struct FileFmt {
	int nfields = 0;
	bool strandParsed = false; /* BED6 / BED6_PLUS record types */
};

static inline bool isStrandText(const char *p, size_t n) {
	return n == 1 && (p[0] == '+' || p[0] == '-' || p[0] == '.' || p[0] == '*');
}
static bool isNumericText(const char *p, size_t n) {
	bool dig = false;
	for (size_t i = 0; i < n; i++) {
		char c = p[i];
		if (!((c >= '0' && c <= '9') || c == '-' || c == '.' || c == '+' || c == 'e' || c == 'E')) return false;
		dig |= (c >= '0' && c <= '9');
	}
	return dig;
}

/* Conservative superset of the x86 isHeaderLine() test. */
static inline bool maybeHeader(const char *p, size_t len) {
	char c0 = p[0];
	if (c0 == '#' || c0 == '>' || c0 == '!') return true;
	if (len > 4) {
		uint32_t peek;
		memcpy(&peek, p, 4);
		peek |= 0x20202020u;
		switch (peek) {
		case 0x6f726863u: case 0x20726863u: case 0x09726863u:
		case 0x776f7262u: case 0x63617274u: case 0x69736976u:
			return true;
		default: break;
		}
	}
	return false;
}

/* Text sanity for the upstream type-checker's scan window (first bytes of the file). */
static bool textPrefixOk(const char *data, size_t avail) {
	const unsigned char *p = (const unsigned char *)data;
	size_t n = avail < 65536 ? avail : 65536;
	if (n == 0) return fallback("empty file");
	if (p[0] == 0x1f) return fallback("gzip/bgzf magic");
	if (n >= 4 && (!memcmp(p, "BAM\1", 4) || !memcmp(p, "CRAM", 4))) return fallback("BAM/CRAM magic");
	for (size_t i = 0; i < n; i++) {
		unsigned char c = p[i];
		if (!((c >= 0x20 && c < 0x7f) || c == '\t' || c == '\n')) return fallback("non-printable byte in scan window");
	}
	return true;
}

struct InFile {
	int fd = -1;
	size_t size = 0;
	bool open(const std::string &path) {
		fd = ::open(path.c_str(), O_RDONLY);
		if (fd < 0) return fallback("open failed");
		struct stat st;
		if (fstat(fd, &st) != 0 || !S_ISREG(st.st_mode)) return fallback("not a regular file");
		if (st.st_size <= 0) return fallback("empty file");
		size = (size_t)st.st_size;
		return true;
	}
	/* read [off, off+n) fully; returns bytes read */
	size_t readAt(char *dst, size_t n, size_t off) const {
		size_t got = 0;
		while (got < n) {
			ssize_t r = pread(fd, dst + got, n - got, (off_t)(off + got));
			if (r < 0) { if (errno == EINTR) continue; break; }
			if (r == 0) break;
			got += (size_t)r;
		}
		return got;
	}
};

/* Parse an unsigned canonical decimal field terminated by `term` or by `le`. */
static inline bool parseCoord(const char *&q, const char *le, pos_t &v) {
	const char *s = q;
	uint64_t x = 0;
	while (q < le && (unsigned)(*q - '0') < 10u) { x = x * 10 + (uint64_t)(*q - '0'); q++; }
	size_t n = q - s;
	if (n == 0 || n > 18) return false;
	if (n > 1 && *s == '0') return false;
	v = (pos_t)x;
	return true;
}

struct LineRec {
	const char *line;
	size_t len;
	const char *chrom;
	size_t chromLen;
	pos_t start, end; /* original coordinates */
	int strand;
};

/* Scan from q (just after the end coordinate) to the end of the line.
 * Returns the newline position (or EOF) and the number of tabs before it;
 * false if a NUL byte occurs inside the data (upstream's strchr-based field
 * splitting would stop there).  Reads past EOF hit the zero guard page. */
static inline bool scanRest(const char *q, const char *end, const char *&nl, int &tabs) {
	const __m128i NL = _mm_set1_epi8('\n'), TB = _mm_set1_epi8('\t'), Z = _mm_setzero_si128();
	int t = 0;
	const char *r = q;
	for (;;) {
		__m128i v = _mm_loadu_si128((const __m128i *)r);
		unsigned mn = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, NL));
		unsigned mz = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, Z));
		unsigned mt = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, TB));
		unsigned stop = mn | mz;
		if (stop) {
			unsigned k = (unsigned)__builtin_ctz(stop);
			t += __builtin_popcount(mt & ((1u << k) - 1u));
			const char *x = r + k;
			if (x >= end) { nl = end; tabs = t; return true; }
			if (*x == '\0') return false;
			nl = x; tabs = t; return true;
		}
		t += __builtin_popcount(mt);
		r += 16;
	}
}

/* Eight ASCII digits (first digit in the lowest byte) to an integer (Lemire). */
static inline uint64_t parse8(uint64_t val) {
	val = (val & 0x0F0F0F0F0F0F0F0Full) * 2561 >> 8;
	val = (val & 0x00FF00FF00FF00FFull) * 6553601 >> 16;
	return (val & 0x0000FFFF0000FFFFull) * 42949672960001ull >> 32;
}

/* Parse one canonical unsigned decimal; q is left on the first non-digit.
 * Reads 16 bytes at q (safe: guard page). */
static inline bool fastCoord(const char *&q, pos_t &v) {
	__m128i x = _mm_loadu_si128((const __m128i *)q);
	__m128i t = _mm_sub_epi8(x, _mm_set1_epi8('0'));
	__m128i isd = _mm_cmpeq_epi8(_mm_min_epu8(t, _mm_set1_epi8(9)), t);
	unsigned m = (unsigned)_mm_movemask_epi8(isd);
	unsigned nd = (unsigned)__builtin_ctz(~m); /* m < 2^16, so ~m has a bit below 17 */
	if (nd == 0 || (nd > 1 && *q == '0')) return false;
	uint64_t val;
	if (nd <= 8) {
		memcpy(&val, q, 8);
		val = parse8(val << (8 * (8 - nd)));
	} else if (nd <= 16) {
		uint64_t hi, lo;
		memcpy(&hi, q, 8);
		memcpy(&lo, q + nd - 8, 8);
		unsigned nh = nd - 8;
		val = parse8(hi << (8 * (8 - nh))) * 100000000ull + parse8(lo);
	} else {
		const char *ds = q;
		uint64_t y = 0;
		const char *r = q;
		while ((unsigned char)(*r - '0') < 10u) { y = y * 10 + (uint64_t)(*r - '0'); r++; }
		if (r - ds > 18) return false;
		v = (pos_t)y; q = r;
		return true;
	}
	v = (pos_t)val;
	q += nd;
	return true;
}

/* Parse a field of known length nd as a canonical unsigned decimal. */
static inline bool knownCoord(const char *q, size_t nd, pos_t &v) {
	if (nd == 0 || nd > 18 || (nd > 1 && q[0] == '0')) return false;
	if (nd <= 16) {
		__m128i x = _mm_loadu_si128((const __m128i *)q);
		__m128i t = _mm_sub_epi8(x, _mm_set1_epi8('0'));
		unsigned m = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_min_epu8(t, _mm_set1_epi8(9)), t));
		unsigned need = (1u << nd) - 1u;
		if ((m & need) != need) return false;
		uint64_t val;
		if (nd <= 8) {
			memcpy(&val, q, 8);
			val = parse8(val << (8 * (8 - nd)));
		} else {
			uint64_t hi, lo;
			memcpy(&hi, q, 8);
			memcpy(&lo, q + nd - 8, 8);
			val = parse8(hi << (8 * (16 - nd))) * 100000000ull + parse8(lo);
		}
		v = (pos_t)val;
		return true;
	}
	uint64_t y = 0;
	for (size_t i = 0; i < nd; i++) {
		unsigned d = (unsigned char)q[i] - '0';
		if (d > 9) return false;
		y = y * 10 + d;
	}
	v = (pos_t)y;
	return true;
}

/* Does the line at p start with the previous line's "chrom\t" (prevL+1 bytes at prevC)? */
static inline bool sameChromPrefix(const char *p, const char *prevC, size_t prevL) {
	if (prevL < 16) {
		__m128i a = _mm_loadu_si128((const __m128i *)p), b = _mm_loadu_si128((const __m128i *)prevC);
		unsigned eq = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(a, b));
		unsigned need = (1u << (prevL + 1)) - 1u;
		return (eq & need) == need;
	}
	return memcmp(p, prevC, prevL + 1) == 0;
}

/* Iterate over all lines; determine the format from the first line; call
 * sink(LineRec&, chromId) for every record.  Returns false (fallback) on any
 * line the fast path does not accept. */
/* Stage 1 of the line scanner: append offsets (relative to p) of tabs to tx and
 * of newlines to nx for bytes [off, lim); false if a NUL byte occurs. */
static bool stage1_sse2(const char *p, size_t off, size_t lim, uint32_t *tx, size_t &tc, uint32_t *nx, size_t &nc) {
	const __m128i TB = _mm_set1_epi8('\t'), NL = _mm_set1_epi8('\n'), Z = _mm_setzero_si128();
	uint64_t bad = 0;
	for (size_t i = off; i < lim; i += 64) {
		const char *u = p + i;
		_mm_prefetch(u + 8192, _MM_HINT_T0);
		uint64_t mt = 0, mn = 0, mz = 0;
		for (int w = 0; w < 4; w++) {
			__m128i v = _mm_loadu_si128((const __m128i *)(u + 16 * w));
			mt |= (uint64_t)(unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, TB)) << (16 * w);
			mn |= (uint64_t)(unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, NL)) << (16 * w);
			mz |= (uint64_t)(unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, Z)) << (16 * w);
		}
		if (lim - i < 64) {
			uint64_t keep = ((uint64_t)1 << (lim - i)) - 1;
			mt &= keep; mn &= keep; mz &= keep;
		}
		bad |= mz;
		const uint32_t base = (uint32_t)i;
		while (mt) { tx[tc++] = base + (uint32_t)__builtin_ctzll(mt); mt &= mt - 1; }
		while (mn) { nx[nc++] = base + (uint32_t)__builtin_ctzll(mn); mn &= mn - 1; }
	}
	return bad == 0;
}

#if defined(__GNUC__) && defined(__x86_64__)
#include <immintrin.h>
__attribute__((target("avx2,bmi"))) static bool stage1_avx2(const char *p, size_t off, size_t lim, uint32_t *tx, size_t &tc, uint32_t *nx, size_t &nc) {
	const __m256i TB = _mm256_set1_epi8('\t'), NL = _mm256_set1_epi8('\n'), Z = _mm256_setzero_si256();
	__m256i mins = _mm256_set1_epi8((char)0xff);
	size_t i = off;
	const size_t full = off + ((lim - off) & ~(size_t)63);
	for (; i < full; i += 64) {
		_mm_prefetch(p + i + 8192, _MM_HINT_T0);
		__m256i v0 = _mm256_loadu_si256((const __m256i *)(p + i));
		__m256i v1 = _mm256_loadu_si256((const __m256i *)(p + i + 32));
		uint64_t mt = (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v0, TB)) |
		              ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v1, TB)) << 32);
		uint64_t mn = (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v0, NL)) |
		              ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v1, NL)) << 32);
		mins = _mm256_min_epu8(mins, _mm256_min_epu8(v0, v1));
		const uint32_t base = (uint32_t)i;
		while (mt) { tx[tc++] = base + (uint32_t)_tzcnt_u64(mt); mt &= mt - 1; }
		while (mn) { nx[nc++] = base + (uint32_t)_tzcnt_u64(mn); mn &= mn - 1; }
	}
	if ((unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(mins, Z)) != 0) return false;
	if (i < lim) return stage1_sse2(p, i, lim, tx, tc, nx, nc);
	return true;
}
/* SWIFTCOV_NO_AVX2 forces the SSE2 scanner (testing) */
static bool g_haveAvx2 = !getenv("SWIFTCOV_NO_AVX2") && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("bmi");
/* per-line scanner needs AVX2, BMI1/2 and POPCNT; SWIFTCOV_BLOCK_SCAN forces the block scanner (testing) */
static bool g_haveLineScan = g_haveAvx2 && !getenv("SWIFTCOV_BLOCK_SCAN") && __builtin_cpu_supports("bmi2") && __builtin_cpu_supports("popcnt");
#endif

/* number of v[i] <= x for i in [i0, i1) (window scans of the track index) */
static inline uint64_t countLE_scalar(const pos_t *v, size_t i0, size_t i1, pos_t x) {
	uint64_t c = 0;
	for (size_t i = i0; i < i1; i++) c += (v[i] <= x);
	return c;
}
#if defined(__GNUC__) && defined(__x86_64__)
__attribute__((target("avx2"))) static uint64_t countLE_avx2(const pos_t *v, size_t i0, size_t i1, pos_t x) {
	const __m256i X = _mm256_set1_epi64x((long long)x);
	__m256i acc = _mm256_setzero_si256();
	size_t i = i0;
	for (; i + 4 <= i1; i += 4) acc = _mm256_sub_epi64(acc, _mm256_cmpgt_epi64(_mm256_loadu_si256((const __m256i *)(v + i)), X));
	long long t[4];
	_mm256_storeu_si256((__m256i *)t, acc);
	uint64_t gt = (uint64_t)(t[0] + t[1] + t[2] + t[3]);
	return (uint64_t)(i - i0) - gt + countLE_scalar(v, i, i1, x);
}
#endif
static inline uint64_t countLE32_scalar(const rec_t *v, size_t i0, size_t i1, rec_t x) {
	uint64_t c = 0;
	for (size_t i = i0; i < i1; i++) c += (v[i] <= x);
	return c;
}
#if defined(__GNUC__) && defined(__x86_64__)
__attribute__((target("avx2"))) static uint64_t countLE32_avx2(const rec_t *v, size_t i0, size_t i1, rec_t x) {
	const __m256i X = _mm256_set1_epi32(x);
	__m256i acc = _mm256_setzero_si256();
	size_t i = i0;
	for (; i + 8 <= i1; i += 8) acc = _mm256_sub_epi32(acc, _mm256_cmpgt_epi32(_mm256_loadu_si256((const __m256i *)(v + i)), X));
	int32_t t[8];
	_mm256_storeu_si256((__m256i *)t, acc);
	uint64_t gt = 0;
	for (int k = 0; k < 8; k++) gt += (uint64_t)(uint32_t)t[k];
	return (uint64_t)(i - i0) - gt + countLE32_scalar(v, i, i1, x);
}
#endif
/* number of v[i] <= x for 32-bit stored coordinates and a 64-bit bound */
static inline uint64_t countLE(const rec_t *v, size_t i0, size_t i1, pos_t x) {
	if (x >= (pos_t)INT32_MAX) return (uint64_t)(i1 - i0);
	if (x < (pos_t)INT32_MIN) return 0;
#if defined(__GNUC__) && defined(__x86_64__)
	if (i1 - i0 >= 16 && g_haveAvx2) return countLE32_avx2(v, i0, i1, (rec_t)x);
#endif
	return countLE32_scalar(v, i0, i1, (rec_t)x);
}
static inline uint64_t countLE(const pos_t *v, size_t i0, size_t i1, pos_t x) {
#if defined(__GNUC__) && defined(__x86_64__)
	if (i1 - i0 >= 16 && g_haveAvx2) return countLE_avx2(v, i0, i1, x);
#endif
	return countLE_scalar(v, i0, i1, x);
}

static inline bool stage1(const char *p, size_t off, size_t lim, uint32_t *tx, size_t &tc, uint32_t *nx, size_t &nc) {
#if defined(__GNUC__) && defined(__x86_64__)
	if (g_haveAvx2) return stage1_avx2(p, off, lim, tx, tc, nx, nc);
#endif
	return stage1_sse2(p, off, lim, tx, tc, nx, nc);
}

struct BedParser {
	FileFmt fmt;
	bool first = true;
	bool needStrand = false;
	ChromTable *chroms = nullptr;
	std::vector<char> prevPad; /* previous chromosome name + '\t' + zero padding */
	size_t prevL = (size_t)-1; /* no previous chromosome yet */
	int prevId = -1;
	std::vector<uint32_t> idx, nidx;
	void setPrev(const char *c, size_t n, int id) {
		prevPad.assign(n + 40, 0);
		memcpy(prevPad.data(), c, n);
		prevPad[n] = '\t';
		prevL = n; prevId = id;
	}
#if defined(__GNUC__) && defined(__x86_64__)
	/* Per-line scanner (AVX2): one 64-byte window at each line start gives the
	 * newline (first set bit of the newline mask) and every tab before it, so the
	 * only serial dependency between lines is load -> compare -> tzcnt -> next
	 * line start; field parsing overlaps with the next line's scan.  Lines of 64
	 * bytes or more continue window by window.  Same checks as the block scanner:
	 * exactly nf-1 tabs per line, no NUL byte, no CR, canonical coordinates.
	 * Needs 64 readable bytes after end (guard page / padded copies). */
	template <class Sink>
	__attribute__((target("avx2,bmi,bmi2,popcnt"))) bool linesAvx2(const char *p, const char *end, Sink &sink) {
		const int nf = fmt.nfields;
		const unsigned ntab = (unsigned)nf - 1;
		const bool wantStrand = needStrand && fmt.strandParsed;
		const __m256i TB = _mm256_set1_epi8('\t'), NL = _mm256_set1_epi8('\n'), Z = _mm256_setzero_si256();
		uint64_t bad = 0;
		__m256i mins = _mm256_set1_epi8((char)0xff);
		/* current chromosome prefix ("name\t") as a vector and its byte mask */
		__m256i cv = _mm256_setzero_si256();
		uint32_t cneed = 0;
#define SWIFTCOV_SETVEC() do { \
			if (prevL != (size_t)-1 && prevL + 1 <= 32) { \
				cv = _mm256_loadu_si256((const __m256i *)prevPad.data()); \
				cneed = prevL + 1 == 32 ? 0xffffffffu : ((1u << (prevL + 1)) - 1u); \
			} else cneed = 0; \
		} while (0)
		SWIFTCOV_SETVEC();
		LineRec r;
		while (p < end) {
			_mm_prefetch(p + 4096, _MM_HINT_T0);
			const __m256i v0 = _mm256_loadu_si256((const __m256i *)p);
			const __m256i v1 = _mm256_loadu_si256((const __m256i *)(p + 32));
			uint64_t mn = (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v0, NL)) |
			              ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v1, NL)) << 32);
			uint64_t mt = (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v0, TB)) |
			              ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v1, TB)) << 32);
			size_t nlo, t0, t1, t2 = 0, t4 = 0, t5 = 0;
			if (__builtin_expect(mn != 0, 1)) {
				nlo = (size_t)_tzcnt_u64(mn);
				const uint64_t below = _bzhi_u64(~0ull, (unsigned)nlo);
				mt &= below;
				/* NUL check: windows inside the block fold into a running byte
				 * minimum; a window reaching past end checks only its line */
				if (__builtin_expect(p + 64 <= end, 1)) mins = _mm256_min_epu8(mins, _mm256_min_epu8(v0, v1));
				else bad |= ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v0, Z)) |
				             ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v1, Z)) << 32)) & below;
				if ((unsigned)_mm_popcnt_u64(mt) != ntab) return fallback("inconsistent field count");
				t0 = (size_t)_tzcnt_u64(mt); mt = _blsr_u64(mt);
				t1 = (size_t)_tzcnt_u64(mt); mt = _blsr_u64(mt);
				t2 = (size_t)_tzcnt_u64(mt);
				if (wantStrand) {
					mt = _blsr_u64(mt); mt = _blsr_u64(mt);
					t4 = (size_t)_tzcnt_u64(mt); mt = _blsr_u64(mt);
					t5 = (size_t)_tzcnt_u64(mt);
				}
			} else {
				/* long line: continue in 64-byte windows, keeping the first six tabs */
				uint64_t mz = (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v0, Z)) |
				              ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v1, Z)) << 32);
				size_t o[6] = {0, 0, 0, 0, 0, 0};
				unsigned cnt = 0;
				size_t wb = 0;
				for (;;) {
					if (mn) {
						const unsigned nl = (unsigned)_tzcnt_u64(mn);
						const uint64_t below = _bzhi_u64(~0ull, nl);
						mt &= below; mz &= below;
					}
					bad |= mz;
					while (mt && cnt < 6) { o[cnt++] = wb + (size_t)_tzcnt_u64(mt); mt = _blsr_u64(mt); }
					cnt += (unsigned)_mm_popcnt_u64(mt);
					if (mn) { nlo = wb + (size_t)_tzcnt_u64(mn); break; }
					wb += 64;
					const __m256i w0 = _mm256_loadu_si256((const __m256i *)(p + wb));
					const __m256i w1 = _mm256_loadu_si256((const __m256i *)(p + wb + 32));
					mn = (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(w0, NL)) |
					     ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(w1, NL)) << 32);
					mt = (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(w0, TB)) |
					     ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(w1, TB)) << 32);
					mz = (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(w0, Z)) |
					     ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(w1, Z)) << 32);
				}
				if (cnt != ntab) return fallback("inconsistent field count");
				t0 = o[0]; t1 = o[1]; t2 = o[2]; t4 = o[4]; t5 = o[5];
			}
			const size_t chromLen = t0;
			const bool same = cneed ? (((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v0, cv)) & cneed) == cneed)
			                        : (chromLen == prevL && sameChromPrefix(p, prevPad.data(), prevL));
			if (!same) {
				if (nlo < 5) return fallback("short/empty line");
				if (maybeHeader(p, nlo)) return fallback("header-like line");
				if (chromLen == 0) return fallback("missing chrom/tab");
				const int c = chroms->get(p, chromLen);
				setPrev(p, chromLen, c);
				SWIFTCOV_SETVEC();
			}
			const size_t endTerm = ntab >= 3 ? t2 : nlo;
			pos_t s, e;
			const size_t n1 = t1 - t0 - 1, n2 = endTerm - t1 - 1;
			if (__builtin_expect(n1 - 1 < 8 && n2 - 1 < 8 && endTerm <= 32, 1)) {
				/* both coordinates (1-8 digits) inside the first 32 bytes: one digit
				 * mask validates both fields, and one SIMD multiply-add chain converts
				 * them together (each number right-aligned in a 64-bit lane) */
				const __m256i dz = _mm256_sub_epi8(v0, _mm256_set1_epi8('0'));
				const uint32_t dm = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_min_epu8(dz, _mm256_set1_epi8(9)), dz));
				const uint32_t want = (uint32_t)(_bzhi_u64(~0ull, (unsigned)endTerm) & ~_bzhi_u64(~0ull, (unsigned)t0 + 1)) & ~(1u << t1);
				if ((dm & want) != want) return fallback("coordinate field");
				if ((p[t0 + 1] == '0' && n1 > 1) || (p[t1 + 1] == '0' && n2 > 1)) return fallback("coordinate field");
				uint64_t wa, wb;
				memcpy(&wa, p + t0 + 1, 8);
				memcpy(&wb, p + t1 + 1, 8);
				__m128i x = _mm_set_epi64x((long long)wb, (long long)wa);
				x = _mm_sllv_epi64(x, _mm_set_epi64x((long long)(64 - 8 * n2), (long long)(64 - 8 * n1)));
				x = _mm_and_si128(x, _mm_set1_epi8(0x0f));
				x = _mm_maddubs_epi16(x, _mm_set1_epi16(0x010a));  /* byte pairs: 10*a + b */
				x = _mm_madd_epi16(x, _mm_set1_epi32(0x00010064));  /* 16-bit pairs: 100*a + b */
				x = _mm_packus_epi32(x, x);
				x = _mm_madd_epi16(x, _mm_set1_epi32(0x00012710));  /* 10000*a + b */
				const uint64_t se = (uint64_t)_mm_cvtsi128_si64(x);
				s = (pos_t)(uint32_t)se;
				e = (pos_t)(se >> 32);
			} else {
				if (!knownCoord(p + t0 + 1, n1, s)) return fallback("start field");
				if (!knownCoord(p + t1 + 1, n2, e)) return fallback("end field");
			}
			if (p[nlo - 1] == '\r') return fallback("CR line ending");
			if (e < s) return fallback("end < start");
			r.line = p; r.len = nlo; r.chrom = p; r.chromLen = chromLen;
			r.start = s; r.end = e;
			r.strand = STR_UNKNOWN;
			if (wantStrand) {
				const size_t f6e = nf > 6 ? t5 : nlo;
				if (f6e - t4 == 2) {
					const char c = p[t4 + 1];
					if (c == '+') r.strand = STR_FWD;
					else if (c == '-') r.strand = STR_REV;
				}
			}
			if (!sink(r, prevId)) return false;
			p += nlo + 1;
		}
#undef SWIFTCOV_SETVEC
		if (bad || (unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(mins, Z)) != 0) return fallback("NUL byte");
		return true;
	}
#endif
	/* Process the complete lines in [p, end) (end[-1] == '\n'; at least 16
	 * readable bytes follow end). */
	template <class Sink>
	bool block(const char *p, const char *end, Sink &sink) {
	LineRec r;
	int cid;
	/* ---- first line: generic parse establishes the record type ---- */
	if (first) {
		const char *nl = (const char *)memchr(p, '\n', end - p);
		const char *le = nl ? nl : end;
		size_t len = le - p;
		if (len < 5) return fallback("short/empty line");
		if (le[-1] == '\r') return fallback("CR line ending");
		/* upstream's type checker re-tokenises the first data line with getline(),
		 * which drops an empty last field, and then exits with "wrong number of
		 * fields"; leave that case to upstream */
		if (le[-1] == '\t') return fallback("first line ends with a tab");
		if (maybeHeader(p, len)) return fallback("header-like line");
		if (memchr(p, '\0', len)) return fallback("NUL byte");
		const char *t1 = (const char *)memchr(p, '\t', len);
		if (!t1 || t1 == p) return fallback("missing chrom/tab");
		const char *q = t1 + 1;
		pos_t s, e;
		if (!parseCoord(q, le, s)) return fallback("start field");
		if (q >= le || *q != '\t') return fallback("start terminator");
		q++;
		if (!parseCoord(q, le, e)) return fallback("end field");
		int nf = 3;
		const char *strandP = nullptr; size_t strandN = 0;
		const char *f6 = nullptr, *f7 = nullptr, *f9 = nullptr; size_t f6n = 0, f7n = 0, f9n = 0;
		if (q < le) {
			if (*q != '\t') return fallback("end terminator");
			const char *fs = q + 1;
			for (const char *r2 = q + 1;; r2++) {
				if (r2 == le || *r2 == '\t') {
					if (nf == 5) { strandP = fs; strandN = r2 - fs; }
					else if (nf == 6) { f6 = fs; f6n = r2 - fs; }
					else if (nf == 7) { f7 = fs; f7n = r2 - fs; }
					else if (nf == 9) { f9 = fs; f9n = r2 - fs; }
					nf++;
					if (r2 == le) break;
					fs = r2 + 1;
				}
			}
		}
		if (e < s) return fallback("end < start");
		fmt.nfields = nf;
		bool strandLike = nf >= 6 && isStrandText(strandP, strandN);
		if (nf == 12 && strandLike && isNumericText(f6, f6n) && isNumericText(f7, f7n) && isNumericText(f9, f9n))
			return fallback("BED12 record type");
		fmt.strandParsed = strandLike;
		if (needStrand && !fmt.strandParsed) return fallback("strand requested but record type has no parsed strand");
		r.line = p; r.len = len; r.chrom = p; r.chromLen = t1 - p;
		r.start = s; r.end = e;
		r.strand = STR_UNKNOWN;
		if (fmt.strandParsed && strandN == 1) {
			if (strandP[0] == '+') r.strand = STR_FWD;
			else if (strandP[0] == '-') r.strand = STR_REV;
		}
		cid = chroms->get(r.chrom, r.chromLen);
		setPrev(r.chrom, r.chromLen, cid);
		first = false;
		if (!sink(r, cid)) return false;
		p = le + 1;
	}
	/* ---- remaining lines: structural index ----
	 * Stage 1 records the offsets of every tab and every newline in separate
	 * arrays (AVX2 or SSE2 compare masks over 64-byte units); stage 2 takes each
	 * line's nf-1 tabs from the tab array and checks that exactly that many lie
	 * between the previous and the current newline, so successive lines carry no
	 * serial scanning dependency. */
	if (p >= end) return true;
#if defined(__GNUC__) && defined(__x86_64__)
	if (g_haveLineScan) return linesAvx2(p, end, sink);
#endif
	const size_t n = (size_t)(end - p);
	if (n >= 0xffffff00u) return fallback("block too large");
	const int nf = fmt.nfields;
	const size_t ntab = (size_t)nf - 1;
	const bool wantStrand = needStrand && fmt.strandParsed;
	const size_t SUB = 65536;
	if (idx.size() < SUB + 128 + (size_t)nf) idx.resize(SUB + 128 + (size_t)nf);
	if (nidx.size() < SUB + 128) nidx.resize(SUB + 128);
	uint32_t *tx = idx.data(), *nx = nidx.data();
	size_t tc = 0;      /* tab entries (including carried tabs of a partial line) */
	size_t ls = 0;      /* offset of the current line start */
	for (size_t off = 0; off < n; off += SUB) {
		const size_t lim = off + SUB < n ? off + SUB : n;
		size_t nc = 0;
		if (!stage1(p, off, lim, tx, tc, nx, nc)) return fallback("NUL byte");
		size_t tk = 0;
		for (size_t j = 0; j < nc; j++) {
			const size_t nlo = nx[j];
			if (tk + ntab > tc) return fallback("inconsistent field count");
			const uint32_t *o = tx + tk;
			if (o[ntab - 1] > nlo) return fallback("inconsistent field count");
			if (tk + ntab < tc && o[ntab] < nlo) return fallback("inconsistent field count");
			const char *line = p + ls;
			const size_t chromLen = o[0] - ls;
			if (!(chromLen == prevL && sameChromPrefix(line, prevPad.data(), prevL))) {
				size_t len = nlo - ls;
				if (len < 5) return fallback("short/empty line");
				if (maybeHeader(line, len)) return fallback("header-like line");
				if (chromLen == 0 || o[0] < ls) return fallback("missing chrom/tab");
				cid = chroms->get(line, chromLen);
				setPrev(line, chromLen, cid);
			}
			const size_t endTerm = nf > 3 ? o[2] : nlo;
			pos_t s, e;
			if (!knownCoord(p + o[0] + 1, o[1] - o[0] - 1, s)) return fallback("start field");
			if (!knownCoord(p + o[1] + 1, endTerm - o[1] - 1, e)) return fallback("end field");
			if (p[nlo - 1] == '\r') return fallback("CR line ending");
			if (e < s) return fallback("end < start");
			r.line = line; r.len = nlo - ls; r.chrom = line; r.chromLen = chromLen;
			r.start = s; r.end = e;
			r.strand = STR_UNKNOWN;
			if (wantStrand) {
				const size_t f6e = nf > 6 ? o[5] : nlo;
				if (f6e - o[4] == 2) {
					char c = p[o[4] + 1];
					if (c == '+') r.strand = STR_FWD;
					else if (c == '-') r.strand = STR_REV;
				}
			}
			if (!sink(r, prevId)) return false;
			tk += ntab;
			ls = nlo + 1;
		}
		/* carry the partial line's tabs; more than nf-1 means too many fields */
		size_t carry = tc - tk;
		if (carry > ntab) return fallback("inconsistent field count");
		for (size_t j = 0; j < carry; j++) tx[j] = tx[tk + j];
		tc = carry;
	}
	if (tc != 0 || ls != n) return fallback("inconsistent field count");
	return true;
}
};

/* Sort-order validation for -sorted (see header comment).  Each parser records
 * its chromosome runs (first/last original and widened starts) and checks start
 * order inside a run; validateSorted() checks the concatenated runs of one file:
 * starts across part boundaries, contiguity and chromosome order. */
struct RunInfo { int cid; pos_t firstS, firstAdj, lastS, lastAdj; };
struct SortTrace {
	bool enabled = false;
	std::vector<RunInfo> runs;
	inline bool on(int cid, pos_t s, pos_t adj) {
		if (!enabled) return true;
		if (runs.empty() || runs.back().cid != cid) { runs.push_back(RunInfo{cid, s, adj, s, adj}); return true; }
		RunInfo &r = runs.back();
		if (s < r.lastS || adj < r.lastAdj) return fallback("-sorted: start order");
		r.lastS = s; r.lastAdj = adj;
		return true;
	}
};
static bool validateSorted(const std::vector<RunInfo> &in, const ChromTable &tab,
                           const std::unordered_map<std::string, int> *genome) {
	std::vector<RunInfo> runs;
	for (const RunInfo &r : in) {
		if (!runs.empty() && runs.back().cid == r.cid) {
			RunInfo &q = runs.back();
			if (r.firstS < q.lastS || r.firstAdj < q.lastAdj) return fallback("-sorted: start order");
			q.lastS = r.lastS; q.lastAdj = r.lastAdj;
		} else runs.push_back(r);
	}
	std::vector<char> seen(tab.names.size(), 0);
	int prevGid = -1;
	for (size_t i = 0; i < runs.size(); i++) {
		int c = runs[i].cid;
		if (seen[c]) return fallback("-sorted: chromosome not contiguous");
		seen[c] = 1;
		if (genome) {
			auto it = genome->find(tab.names[c]);
			if (it == genome->end()) return fallback("-sorted: chromosome missing from genome file");
			if (i > 0 && it->second <= prevGid) return fallback("-sorted: genome order");
			prevGid = it->second;
		} else if (i > 0) {
			if (!(tab.names[runs[i - 1].cid] < tab.names[c])) return fallback("-sorted: not lexicographic");
		}
	}
	return true;
}

/* Unsorted (bin tree) path: coordinates comfortably inside the bin range. */
static const pos_t MAX_UNSORTED_COORD = (pos_t)1 << 33;

/* ------------------------------------------------------------------ */
/* Database tracks                                                     */
/* ------------------------------------------------------------------ */

template <class T>
static void radixSortPos(std::vector<T> &v) {
	size_t n = v.size();
	if (n < 2) return;
	bool sorted = true;
	for (size_t i = 1; i < n; i++) if (v[i] < v[i - 1]) { sorted = false; break; }
	if (sorted) return;
	if (n < 4096) { std::sort(v.begin(), v.end()); return; }
	pos_t mn = v[0], mx = v[0];
	for (pos_t x : v) { if (x < mn) mn = x; if (x > mx) mx = x; }
	uint64_t range = (uint64_t)mx - (uint64_t)mn;
	std::vector<uint64_t> a(n), b(n);
	for (size_t i = 0; i < n; i++) a[i] = (uint64_t)(pos_t)v[i] - (uint64_t)mn;
	int shift = 0;
	std::vector<size_t> cnt(1 << 16);
	while (shift < 64 && (range >> shift) != 0) {
		std::fill(cnt.begin(), cnt.end(), 0);
		for (size_t i = 0; i < n; i++) cnt[(a[i] >> shift) & 0xffff]++;
		size_t sum = 0;
		for (size_t k = 0; k < cnt.size(); k++) { size_t t = cnt[k]; cnt[k] = sum; sum += t; }
		for (size_t i = 0; i < n; i++) b[cnt[(a[i] >> shift) & 0xffff]++] = a[i];
		a.swap(b);
		shift += 16;
	}
	for (size_t i = 0; i < n; i++) v[i] = (T)(pos_t)(a[i] + (uint64_t)mn);
}

/* Branch-free binary searches (Khuong & Morin): the loop compiles to a
 * conditional move per step, so queries cost no branch mispredictions.  Arrays
 * beyond L2 keep std:: searches, whose speculation hides memory latency
 * (measured: branch-free +0.3 ms on transcript_default, with or without prefetch). */
static const size_t BL_MAX = (size_t)1 << 18; /* 2 MiB of keys: larger arrays use the branchy search */
template <class T>
static inline size_t lbBL(const T *a, size_t n, pos_t x) {
	if (n > BL_MAX) return (size_t)(std::lower_bound(a, a + n, x, [](pos_t u, pos_t v) { return u < v; }) - a);
	if (n == 0) return 0;
	const T *b = a;
	while (n > 1) { size_t h = n / 2; b = (b[h] < x) ? b + h : b; n -= h; }
	return (size_t)(b - a) + (*b < x);
}
template <class T>
static inline size_t ubBL(const T *a, size_t n, pos_t x) {
	if (n > BL_MAX) return (size_t)(std::upper_bound(a, a + n, x, [](pos_t u, pos_t v) { return u < v; }) - a);
	if (n == 0) return 0;
	const T *b = a;
	while (n > 1) { size_t h = n / 2; b = (b[h] <= x) ? b + h : b; n -= h; }
	return (size_t)(b - a) + (*b <= x);
}

struct Track {
	rec_t *S = nullptr, *E = nullptr; /* paired adjusted coordinates (slice of a flat array or own storage) */
	size_t n = 0;
	std::vector<rec_t> ownS, ownE;
	bool sByStart = true;
	/* derived */
	bool prepared = false;
	const rec_t *Ss = nullptr;      /* sorted starts (== S after sortPairs) */
	std::vector<rec_t> Es;          /* sorted ends */
	std::vector<pos_t> uS, uE;      /* union of intervals */
	std::vector<pos_t> uPre;        /* uPre[k] = total union length before uS[k] */
	bool havePrefix = false;
	std::vector<pos_t> preS, preE;  /* prefix sums of Ss / Es (length n+1) */
	/* interval tree over pairs sorted by start (fraction mode) */
	bool haveTree = false;
	const rec_t *tS = nullptr, *tE = nullptr;
	std::vector<rec_t> tMax;
	int maxLevel = -1;

	void sortPairs() {
		if (sByStart) {
			for (size_t i = 1; i < n; i++) if (S[i] < S[i - 1]) { sByStart = false; break; }
			if (sByStart) return;
		}
		std::vector<std::pair<rec_t, rec_t>> v(n);
		for (size_t i = 0; i < n; i++) v[i] = std::make_pair(S[i], E[i]);
		std::sort(v.begin(), v.end());
		for (size_t i = 0; i < n; i++) { S[i] = v[i].first; E[i] = v[i].second; }
		sByStart = true;
	}

	/* Window mode: with maxSpan M, a record ending at or before x must start at
	 * or before x - M unless it is one of the few starting in (x - M, x]; this
	 * answers end-rank queries from the start-sorted pairs without sorting ends. */
	pos_t M = 0;
	bool windowMode = false;
	static const pos_t WINDOW_MAX_SPAN = (pos_t)1 << 16;
	static const pos_t DIFF_MAX_LEN = (pos_t)1 << 24;

	void prepareSortedEnds() {
		if (!Es.empty() || n == 0) return;
		Es.assign(E, E + n);
		radixSortPos(Es);
	}

	/* One pass over [a, b): start order, max span and (for sorted input) the
	 * interval union; with several threads the chunks are merged afterwards. */
	struct ChunkStat { bool sorted = true; pos_t maxSpan = 0; std::vector<pos_t> us, ue; };
#if defined(__GNUC__) && defined(__x86_64__)
	/* AVX2 version of scanChunk (32-bit records, eight per step).  Start order
	 * compares each vector of starts with the previous-element vector; spans use
	 * a vector max.  Union: a new interval starts at record j iff S[j] exceeds the
	 * maximum of all earlier ends (that running maximum is the current interval's
	 * end), so each block computes the exclusive prefix maximum of its ends, seeded
	 * with the carried maximum, and only a block containing such a gap is handled
	 * by the scalar code.  Checked against scanChunk on 20,000 random ranges
	 * (sorted and unsorted) in a scratch program. */
	__attribute__((target("avx2"))) void scanChunkAvx2(size_t a, size_t b, bool wantUnion, ChunkStat &c) const {
		if (a >= b) return;
		pos_t m = (pos_t)E[a] - S[a], cs = S[a], ce = E[a];
		bool sorted = true;
		size_t i = a + 1;
		const __m256i NEG = _mm256_set1_epi32(INT32_MIN);
		const __m256i I1 = _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6), I2 = _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5);
		const __m256i I4 = _mm256_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3), I7 = _mm256_set1_epi32(7);
		__m256i vm = _mm256_setzero_si256(), bad = _mm256_setzero_si256(), carry = _mm256_set1_epi32((int)ce);
		for (; i + 8 <= b; i += 8) {
			const __m256i vs = _mm256_loadu_si256((const __m256i *)(S + i));
			const __m256i ve = _mm256_loadu_si256((const __m256i *)(E + i));
			const __m256i vp = _mm256_loadu_si256((const __m256i *)(S + i - 1));
			bad = _mm256_or_si256(bad, _mm256_cmpgt_epi32(vp, vs));
			vm = _mm256_max_epi32(vm, _mm256_sub_epi32(ve, vs));
			if (wantUnion) {
				__m256i x = _mm256_max_epi32(ve, _mm256_blend_epi32(_mm256_permutevar8x32_epi32(ve, I1), NEG, 0x01));
				x = _mm256_max_epi32(x, _mm256_blend_epi32(_mm256_permutevar8x32_epi32(x, I2), NEG, 0x03));
				x = _mm256_max_epi32(x, _mm256_blend_epi32(_mm256_permutevar8x32_epi32(x, I4), NEG, 0x0f));
				const __m256i excl = _mm256_blend_epi32(_mm256_permutevar8x32_epi32(x, I1), NEG, 0x01);
				const __m256i pm = _mm256_max_epi32(excl, carry);
				const __m256i gaps = _mm256_cmpgt_epi32(vs, pm);
				if (__builtin_expect(!_mm256_testz_si256(gaps, gaps), 0)) {
					ce = (pos_t)_mm256_extract_epi32(carry, 0);
					for (size_t j = i; j < i + 8; j++) {
						const pos_t si = S[j], ei = E[j];
						if (si <= ce) { if (ei > ce) ce = ei; }
						else { c.us.push_back(cs); c.ue.push_back(ce); cs = si; ce = ei; }
					}
					carry = _mm256_set1_epi32((int)ce);
				} else carry = _mm256_max_epi32(carry, _mm256_permutevar8x32_epi32(x, I7));
			}
		}
		if (wantUnion) ce = (pos_t)_mm256_extract_epi32(carry, 0);
		if (!_mm256_testz_si256(bad, bad)) sorted = false;
		{
			int32_t t[8];
			_mm256_storeu_si256((__m256i *)t, vm);
			for (int k = 0; k < 8; k++) if ((pos_t)t[k] > m) m = t[k];
		}
		for (; i < b; i++) {
			const pos_t si = S[i], ei = E[i];
			if (si < S[i - 1]) sorted = false;
			if (ei - si > m) m = ei - si;
			if (wantUnion) {
				if (si <= ce) { if (ei > ce) ce = ei; }
				else { c.us.push_back(cs); c.ue.push_back(ce); cs = si; ce = ei; }
			}
		}
		if (wantUnion) { c.us.push_back(cs); c.ue.push_back(ce); }
		c.sorted = sorted; c.maxSpan = m;
	}
#endif
	void scanChunk(size_t a, size_t b, bool wantUnion, ChunkStat &c) const {
#if defined(__GNUC__) && defined(__x86_64__)
		if (g_haveAvx2) { scanChunkAvx2(a, b, wantUnion, c); return; }
#endif
		if (a >= b) return;
		pos_t m = 0, prev = S[a], cs = S[a], ce = E[a];
		bool sorted = true;
		for (size_t i = a; i < b; i++) {
			const pos_t si = S[i], ei = E[i];
			if (si < prev) sorted = false;
			prev = si;
			if (ei - si > m) m = ei - si;
			if (wantUnion) {
				if (si <= ce) { if (ei > ce) ce = ei; }
				else { c.us.push_back(cs); c.ue.push_back(ce); cs = si; ce = ei; }
			}
		}
		if (wantUnion) { c.us.push_back(cs); c.ue.push_back(ce); }
		c.sorted = sorted; c.maxSpan = m;
	}
	/* Fast preparation: returns false if the records are not start-sorted (the
	 * caller then uses the general path). */
	bool prepareFast(bool wantUnion, int threads) {
		if (prepared || n == 0) return false;
		int T = (threads > 1 && n >= ((size_t)1 << 20)) ? threads : 1;
		std::vector<ChunkStat> cs((size_t)T);
		auto work = [&](int t) { scanChunk(n * (size_t)t / (size_t)T, n * (size_t)(t + 1) / (size_t)T, wantUnion, cs[(size_t)t]); };
		std::vector<std::thread> th;
		for (int t = 1; t < T; t++) th.emplace_back(work, t);
		work(0);
		for (auto &x : th) x.join();
		pos_t m = 0;
		for (int t = 0; t < T; t++) {
			if (!cs[(size_t)t].sorted) return false;
			size_t b = n * (size_t)t / (size_t)T;
			if (t > 0 && b < n && b > 0 && S[b] < S[b - 1]) return false;
			if (cs[(size_t)t].maxSpan > m) m = cs[(size_t)t].maxSpan;
		}
		sByStart = true;
		Ss = S;
		M = m;
		windowMode = M <= WINDOW_MAX_SPAN && n < ((size_t)1 << 31);
		if (!windowMode) prepareSortedEnds();
		prepared = true;
		if (wantUnion) {
			for (int t = 0; t < T; t++) {
				ChunkStat &c = cs[(size_t)t];
				for (size_t j = 0; j < c.us.size(); j++) {
					if (!uE.empty() && c.us[j] <= uE.back()) { if (c.ue[j] > uE.back()) uE.back() = c.ue[j]; }
					else { uS.push_back(c.us[j]); uE.push_back(c.ue[j]); }
				}
			}
			uPre.resize(uS.size());
			pos_t acc = 0;
			for (size_t k = 0; k < uS.size(); k++) { uPre[k] = acc; acc += uE[k] - uS[k]; }
		}
		return true;
	}

	void prepare(bool wantUnion, bool wantPrefix) {
		if (!prepared) {
			sortPairs();
			Ss = S;
			pos_t m = 0;
			for (size_t i = 0; i < n; i++) { pos_t d = E[i] - S[i]; if (d > m) m = d; }
			M = m;
			windowMode = M <= WINDOW_MAX_SPAN && n < ((size_t)1 << 31);
			if (!windowMode) prepareSortedEnds();
			prepared = true;
		}
		if (wantUnion && uS.empty() && n) {
			pos_t cs = S[0], ce = E[0];
			pos_t acc = 0;
			for (size_t i = 1; i < n; i++) {
				if (S[i] <= ce) { if (E[i] > ce) ce = E[i]; }
				else { uS.push_back(cs); uE.push_back(ce); uPre.push_back(acc); acc += ce - cs; cs = S[i]; ce = E[i]; }
			}
			uS.push_back(cs); uE.push_back(ce); uPre.push_back(acc);
		}
		if (wantPrefix && !havePrefix && !windowMode) {
			preS.assign(n + 1, 0); preE.assign(n + 1, 0);
			for (size_t i = 0; i < n; i++) { preS[i + 1] = preS[i] + Ss[i]; preE[i + 1] = preE[i] + Es[i]; }
			havePrefix = true;
		}
	}

	/* first index with S > x */
	inline size_t ubS(pos_t x) const { return ubBL(Ss, n, x); }
	inline size_t lbS(pos_t x) const { return lbBL(Ss, n, x); }
	/* number of records with E <= x */
	uint64_t endsLE(pos_t x) const {
		if (!windowMode) return (uint64_t)ubBL(Es.data(), Es.size(), x);
		size_t i0 = ubS(x - M), i1 = ubS(x);
		return (uint64_t)i0 + countLE(E, i0, i1, x);
	}

	/* number of records with s < ae and e > as */
	uint64_t count(pos_t as, pos_t ae) const {
		return (uint64_t)lbS(ae) - endsLE(as);
	}
	/* union length in (-inf, x) */
	pos_t unionBefore(pos_t x) const {
		size_t k = ubBL(uS.data(), uS.size(), x);
		if (k == 0) return 0;
		k--;
		pos_t r = uPre[k];
		pos_t e = uE[k] < x ? uE[k] : x;
		if (e > uS[k]) r += e - uS[k];
		return r;
	}
	uint64_t covered(pos_t as, pos_t ae) const {
		if (uS.empty()) return 0;
		return (uint64_t)(unionBefore(ae) - unionBefore(as));
	}
	/* sum over x < X of depth(x) (sorted-ends mode) */
	__int128 depthBefore(pos_t x) const {
		size_t cs = lbBL(Ss, n, x);
		size_t ce = lbBL(Es.data(), Es.size(), x);
		return (__int128)x * (__int128)cs - (__int128)preS[cs] - (__int128)x * (__int128)ce + (__int128)preE[ce];
	}
	uint64_t depthSum(pos_t as, pos_t ae) const {
		if (!windowMode) return (uint64_t)(depthBefore(ae) - depthBefore(as));
		size_t i0 = ubS(as - M), i1 = lbS(ae);
		uint64_t sum = 0;
		for (size_t i = i0; i < i1; i++) {
			if (E[i] <= as) continue;
			pos_t a = S[i] > as ? S[i] : as, b = E[i] < ae ? E[i] : ae;
			sum += (uint64_t)(b - a);
		}
		return sum;
	}
	/* window mode: add +1/-1 events of the records overlapping [as, ae) to dd[0..L] */
	void addDiff(pos_t as, pos_t ae, int32_t *dd) const {
		const size_t L = (size_t)(ae - as);
		size_t i0 = ubS(as - M), i1 = lbS(ae);
		for (size_t i = i0; i < i1; i++) {
			if (E[i] <= as) continue;
			size_t a = S[i] > as ? (size_t)(S[i] - as) : 0;
			size_t b = E[i] < ae ? (size_t)(E[i] - as) : L;
			dd[a]++; dd[b]--;
		}
	}
	template <class CB>
	static void emitDiff(pos_t as, pos_t ae, const int32_t *dd, CB &cb) {
		const size_t L = (size_t)(ae - as);
		int64_t depth = dd[0], cur = depth;
		size_t start = 0;
		for (size_t x = 1; x < L; x++) {
			depth += dd[x];
			if (depth != cur) { cb(as + (pos_t)start, as + (pos_t)x, (uint64_t)cur); start = x; cur = depth; }
		}
		cb(as + (pos_t)start, ae, (uint64_t)cur);
	}
	/* Iterate depth runs over [as, ae): cb(runStart, runEnd, depth). */
	template <class CB>
	void runs(pos_t as, pos_t ae, CB &cb, std::vector<int32_t> &diff) const {
		if (windowMode && ae - as <= DIFF_MAX_LEN) {
			/* local difference array over A */
			const size_t L = (size_t)(ae - as);
			diff.assign(L + 1, 0);
			addDiff(as, ae, diff.data());
			emitDiff(as, ae, diff.data(), cb);
			return;
		}
		size_t i = ubBL(Ss, n, as);
		size_t j = ubBL(Es.data(), Es.size(), as);
		uint64_t depth = (uint64_t)(i - j);
		pos_t cur = as;
		for (;;) {
			pos_t ns = i < n ? Ss[i] : INT64_MAX;
			pos_t ne = j < n ? Es[j] : INT64_MAX;
			pos_t x = ns < ne ? ns : ne;
			if (x >= ae) break;
			uint64_t nd = depth;
			while (i < n && Ss[i] == x) { nd++; i++; }
			while (j < n && Es[j] == x) { nd--; j++; }
			if (nd != depth) {
				if (x > cur) cb(cur, x, depth);
				cur = x;
				depth = nd;
			}
		}
		cb(cur, ae, depth);
	}

	/* Implicit augmented interval tree (cgranges layout, Li 2019). */
	void buildTree() {
		if (haveTree) return;
		sortPairs();
		tS = S; tE = E;
		tMax.assign(n, 0);
		haveTree = true;
		if (n == 0) { maxLevel = -1; return; }
		int64_t i, last_i = 0;
		pos_t last = 0;
		for (i = 0; i < (int64_t)n; i += 2) { last_i = i; last = tMax[i] = tE[i]; }
		int k;
		for (k = 1; (int64_t)1 << k <= (int64_t)n; ++k) {
			int64_t x = (int64_t)1 << (k - 1), i0 = (x << 1) - 1, step = x << 2;
			for (i = i0; i < (int64_t)n; i += step) {
				pos_t el = tMax[i - x];
				pos_t er = i + x < (int64_t)n ? tMax[i + x] : last;
				pos_t e = tE[i];
				if (el > e) e = el;
				if (er > e) e = er;
				tMax[i] = (rec_t)e;
			}
			last_i = (last_i >> k & 1) ? last_i - x : last_i + x;
			if (last_i < (int64_t)n && tMax[last_i] > last) last = tMax[last_i];
		}
		maxLevel = k - 1;
	}
	/* Visit records with s < ae and e > as. */
	template <class CB>
	void overlaps(pos_t as, pos_t ae, CB &cb) const {
		const int64_t n = (int64_t)this->n;
		if (n == 0) return;
		struct Item { int64_t x; int k, w; } stack[128];
		int t = 0;
		stack[t++] = Item{((int64_t)1 << maxLevel) - 1, maxLevel, 0};
		while (t) {
			Item z = stack[--t];
			if (z.k <= 3) {
				int64_t i0 = z.x >> z.k << z.k, i1 = i0 + ((int64_t)1 << (z.k + 1)) - 1;
				if (i1 >= n) i1 = n;
				for (int64_t i = i0; i < i1 && tS[i] < ae; ++i)
					if (as < tE[i]) cb(tS[i], tE[i]);
			} else if (z.w == 0) {
				int64_t y = z.x - ((int64_t)1 << (z.k - 1));
				stack[t++] = Item{z.x, z.k, 1};
				if (y >= n || tMax[y] > as) stack[t++] = Item{y, z.k - 1, 0};
			} else if (z.x < n && tS[z.x] < ae) {
				if (as < tE[z.x]) cb(tS[z.x], tE[z.x]);
				stack[t++] = Item{z.x + ((int64_t)1 << (z.k - 1)), z.k - 1, 0};
			}
		}
	}
};

/* The records of one chromosome (and strand set) as one or more segments, each
 * a contiguous slice of one parse part's arrays, so that multi-part parses need
 * no copy into a single array.  Counts, depth sums and overlap enumeration add
 * up over segments; covered bases use the merged union; depth runs accumulate
 * all segments into one difference array (segments that would need sorted end
 * arrays for -hist/-d are first concatenated into one track). */
struct TGroup {
	std::vector<Track> seg;
	size_t n = 0;
	std::vector<pos_t> uS, uE, uPre; /* merged union (several segments) */
	uint64_t count(pos_t as, pos_t ae) const {
		if (seg.size() == 1) return seg[0].count(as, ae);
		uint64_t c = 0;
		for (const Track &t : seg) c += t.count(as, ae);
		return c;
	}
	pos_t unionBefore(pos_t x) const {
		size_t k = ubBL(uS.data(), uS.size(), x);
		if (k == 0) return 0;
		k--;
		pos_t r = uPre[k];
		pos_t e = uE[k] < x ? uE[k] : x;
		if (e > uS[k]) r += e - uS[k];
		return r;
	}
	uint64_t covered(pos_t as, pos_t ae) const {
		if (seg.size() == 1) return seg[0].covered(as, ae);
		if (uS.empty()) return 0;
		return (uint64_t)(unionBefore(ae) - unionBefore(as));
	}
	uint64_t depthSum(pos_t as, pos_t ae) const {
		if (seg.size() == 1) return seg[0].depthSum(as, ae);
		uint64_t d = 0;
		for (const Track &t : seg) d += t.depthSum(as, ae);
		return d;
	}
	template <class CB>
	void runs(pos_t as, pos_t ae, CB &cb, std::vector<int32_t> &diff) const {
		if (seg.size() == 1) { seg[0].runs(as, ae, cb, diff); return; }
		const size_t L = (size_t)(ae - as);
		diff.assign(L + 1, 0);
		for (const Track &t : seg) t.addDiff(as, ae, diff.data());
		Track::emitDiff(as, ae, diff.data(), cb);
	}
	template <class CB>
	void overlaps(pos_t as, pos_t ae, CB &cb) const {
		for (const Track &t : seg) t.overlaps(as, ae, cb);
	}
	void prepare(bool fraction, bool wantUnion, bool wantPrefix, bool needRuns, pos_t maxALen, int threads) {
		if (fraction) { for (Track &t : seg) t.buildTree(); return; }
		for (Track &t : seg) { t.prepareFast(wantUnion, threads); t.prepare(wantUnion, wantPrefix); }
		bool needEnds = false;
		if (needRuns) for (const Track &t : seg) if (!t.windowMode || maxALen > Track::DIFF_MAX_LEN) needEnds = true;
		if (needEnds && seg.size() > 1) {
			/* concatenate into one track (rare: long spans or very long queries) */
			Track m;
			m.ownS.reserve(n); m.ownE.reserve(n);
			for (const Track &t : seg) {
				m.ownS.insert(m.ownS.end(), t.S, t.S + t.n);
				m.ownE.insert(m.ownE.end(), t.E, t.E + t.n);
			}
			m.S = m.ownS.data(); m.E = m.ownE.data(); m.n = m.ownS.size();
			m.prepareFast(wantUnion, threads);
			m.prepare(wantUnion, wantPrefix);
			seg.clear();
			seg.push_back(std::move(m));
		}
		if (needEnds) seg[0].prepareSortedEnds();
		if (wantUnion && seg.size() > 1) {
			std::vector<std::pair<pos_t, pos_t>> iv;
			for (const Track &t : seg)
				for (size_t j = 0; j < t.uS.size(); j++) iv.push_back(std::make_pair(t.uS[j], t.uE[j]));
			std::sort(iv.begin(), iv.end());
			pos_t acc = 0;
			for (const auto &q : iv) {
				if (!uE.empty() && q.first <= uE.back()) { if (q.second > uE.back()) { acc += q.second - uE.back(); uE.back() = q.second; } }
				else { uS.push_back(q.first); uE.push_back(q.second); uPre.push_back(acc); acc += q.second - q.first; }
			}
		}
	}
};

/* ------------------------------------------------------------------ */
/* Output                                                              */
/* ------------------------------------------------------------------ */

struct Out {
	char *buf;
	size_t cap, len = 0;
	bool failed = false;
	explicit Out(size_t c) : cap(c) { buf = (char *)malloc(cap); }
	~Out() { flush(); free(buf); }
	void writeAll(const char *p, size_t n) {
		while (n && !failed) {
			ssize_t w = ::write(1, p, n);
			if (w < 0) { if (errno == EINTR) continue; failed = true; break; }
			p += w; n -= (size_t)w;
		}
	}
	void flush() { if (len) { writeAll(buf, len); len = 0; } }
	inline char *reserve(size_t n) {
		if (len + n > cap) {
			flush();
			if (n > cap) { free(buf); cap = n * 2; buf = (char *)malloc(cap); }
		}
		return buf + len;
	}
	inline void commit(char *p) { len = (size_t)(p - buf); }
};

/* Growable private output buffer for one block of queries. */
struct Buf {
	char *p = nullptr;
	size_t cap = 0, len = 0;
	inline char *reserve(size_t n) {
		if (len + n > cap) {
			size_t nc = cap ? cap * 2 : ((size_t)1 << 20);
			while (nc < len + n) nc *= 2;
			p = (char *)realloc(p, nc);
			cap = nc;
		}
		return p + len;
	}
	inline void commit(char *q) { len = (size_t)(q - p); }
};

static inline char *putU64(char *p, uint64_t v) {
	char tmp[24];
	int n = 0;
	do { tmp[n++] = (char)('0' + v % 10); v /= 10; } while (v);
	while (n) *p++ = tmp[--n];
	return p;
}

/* printf("%0.7f", (double)f) for finite non-negative floats. */
static char *putF7(char *p, float f) {
	double d = (double)f;
	if (!(d >= 0.0) || d >= 1099511627776.0 /* 2^40 */) {
		char tmp[512];
		int n = snprintf(tmp, sizeof(tmp), "%0.7f", d);
		memcpy(p, tmp, n);
		return p + n;
	}
	uint32_t bits;
	memcpy(&bits, &f, 4);
	uint32_t expo = (bits >> 23) & 0xff;
	uint64_t mant = bits & 0x7fffff;
	int e2;
	if (expo == 0) { e2 = -149; }
	else { mant |= 0x800000; e2 = (int)expo - 150; }
	/* value = mant * 2^e2; want q = round_half_even(value * 1e7) */
	uint64_t q;
	if (mant == 0) q = 0;
	else if (e2 >= 0) {
		q = (mant << e2) * 10000000ull; /* < 2^40 * 1e7 < 2^64 */
	} else {
		int k = -e2;
		unsigned __int128 num = (unsigned __int128)mant * 10000000ull;
		if (k >= 120) q = 0; /* num < 2^48, far below half an ulp */
		else {
			unsigned __int128 qq = num >> k;
			unsigned __int128 rem = num - (qq << k);
			unsigned __int128 half = (unsigned __int128)1 << (k - 1);
			if (rem > half || (rem == half && (qq & 1))) qq++;
			q = (uint64_t)qq;
		}
	}
	uint64_t ip = q / 10000000ull, fp = q % 10000000ull;
	p = putU64(p, ip);
	*p++ = '.';
	for (int i = 6; i >= 0; i--) { p[i] = (char)('0' + fp % 10); fp /= 10; }
	return p + 7;
}

/* ------------------------------------------------------------------ */
/* Main                                                                */
/* ------------------------------------------------------------------ */

struct ARec {
	const char *line;
	uint32_t len;
	int32_t chrom;
	pos_t s, e; /* adjusted */
	int8_t strand;
};

struct Hist {
	std::vector<uint64_t> dense; /* depth < dense.size() */
	std::map<uint64_t, uint64_t> sparse;
	inline void add(uint64_t d, uint64_t n) {
		if (d < 65536) {
			if (d >= dense.size()) dense.resize(d + 1, 0);
			dense[d] += n;
		} else sparse[d] += n;
	}
};

static int run(int argc, char **argv) {
	g_debug = getenv("SWIFTCOV_DEBUG") != nullptr;
	g_prof = getenv("SWIFTCOV_PROFILE") != nullptr;
	g_t0 = g_tl = nowSec();
	Config cfg;
	if (!parseArgs(argc, argv, cfg)) return -1;

	std::unordered_map<std::string, int> genome;
	if (cfg.hasG && !loadGenome(cfg.gPath, genome)) return -1;

	InFile fa;
	if (!fa.open(cfg.aPath)) return -1;
	std::vector<InFile> fbs(cfg.bPaths.size());
	for (size_t i = 0; i < fbs.size(); i++)
		if (!fbs[i].open(cfg.bPaths[i])) return -1;

	prof("open+mmap");
	const bool stranded = cfg.same || cfg.diff;
	const bool fraction = cfg.haveFA || cfg.haveFB;
	ChromTable chroms;

	/* --- query --- */
	std::vector<ARec> A;
	char *abuf = nullptr;
	{
		SortTrace sc;
		sc.enabled = cfg.sorted;
		abuf = (char *)malloc(fa.size + 64);
		if (!abuf) return fallback("allocation"), -1;
		if (fa.readAt(abuf, fa.size, 0) != fa.size) return fallback("short read"), -1;
		memset(abuf + fa.size, 0, 64);
		if (!textPrefixOk(abuf, fa.size)) return -1;
		size_t aend = fa.size;
		if (abuf[aend - 1] != '\n') abuf[aend++] = '\n';
		A.reserve(fa.size / 24 + 16);
		auto sink = [&](LineRec &r, int cid) -> bool {
			pos_t s = r.start, e = r.end;
			if (s == e) { s--; e++; }
			if (!sc.on(cid, r.start, s)) return false;
			if (!cfg.sorted && e > MAX_UNSORTED_COORD) return fallback("unsorted: query coordinate beyond bin range");
			if (r.len > 0x7fffffffu) return fallback("line too long");
			ARec a;
			a.line = r.line; a.len = (uint32_t)r.len; a.chrom = cid; a.s = s; a.e = e; a.strand = (int8_t)r.strand;
			A.push_back(a);
			return true;
		};
		BedParser bp;
		bp.needStrand = stranded; bp.chroms = &chroms;
		if (!bp.block(abuf, abuf + aend, sink)) return -1;
		if (bp.first) return fallback("no records"), -1;
		if (cfg.sorted && !validateSorted(sc.runs, chroms, cfg.hasG ? &genome : nullptr)) return -1;
	}
	prof("parse A");

	/* Union of query intervals per chromosome (widened coordinates): a database
	 * record that overlaps no query interval can never be a hit, so it is
	 * validated like every other line but not stored. */
	struct QUnion { std::vector<pos_t> s, e; };
	std::vector<QUnion> qunion(chroms.names.size());
	{
		std::vector<std::vector<std::pair<pos_t, pos_t>>> iv(chroms.names.size());
		for (const ARec &a : A) iv[(size_t)a.chrom].push_back(std::make_pair(a.s, a.e));
		for (size_t c = 0; c < iv.size(); c++) {
			auto &v = iv[c];
			if (v.empty()) continue;
			bool sorted = true;
			for (size_t i = 1; i < v.size(); i++) if (v[i].first < v[i - 1].first) { sorted = false; break; }
			if (!sorted) std::sort(v.begin(), v.end());
			pos_t cs = v[0].first, ce = v[0].second;
			for (size_t i = 1; i < v.size(); i++) {
				if (v[i].first <= ce) { if (v[i].second > ce) ce = v[i].second; }
				else { qunion[c].s.push_back(cs); qunion[c].e.push_back(ce); cs = v[i].first; ce = v[i].second; }
			}
			qunion[c].s.push_back(cs); qunion[c].e.push_back(ce);
		}
	}
	const size_t nQueryChroms = chroms.names.size();
	prof("query union");

	/* --- database --- */
	/* Records go to flat arrays in file order (one set, or one per strand for
	 * -s/-S); a chromosome stored as one contiguous run is used in place.  Large
	 * files are split at line boundaries into one byte range per thread; part 0
	 * writes straight into the final arrays, later parts are copied after it. */
	struct Run { int cid; size_t b, e; };
	struct FlatSet { BigArray<rec_t> S, E; std::vector<Run> runs; };
	const int nsets = stranded ? 2 : 1;
	struct BPart {
		size_t file = 0, begin = 0, end = 0;
		bool first = true;
		ChromTable chroms;
		FlatSet sets[2];
		SortTrace trace;
		FileFmt fmt;
		bool ok = false;
	};
	std::vector<const char *> fileMap(fbs.size(), nullptr);
	for (size_t fi = 0; fi < fbs.size(); fi++) {
		fileMap[fi] = mapWithGuard(fbs[fi].fd, fbs[fi].size);
		if (!fileMap[fi]) return fallback("mmap failed"), -1;
	}
	std::vector<FileFmt> fileFmt(fbs.size());
	std::vector<BPart> parts;
	{
		/* SWIFTCOV_SPLIT_MIN (bytes) lowers the split threshold for testing only */
		const size_t SPLIT_MIN = getenv("SWIFTCOV_SPLIT_MIN") ? (size_t)atol(getenv("SWIFTCOV_SPLIT_MIN")) : ((size_t)8 << 20);
		for (size_t fi = 0; fi < fbs.size(); fi++) {
			const InFile &f = fbs[fi];
			std::vector<size_t> cuts;
			cuts.push_back(0);
			if (cfg.threads > 1 && f.size >= SPLIT_MIN) {
				/* record type from the first line, exactly as part 0 will parse it */
				size_t hn = f.size < ((size_t)1 << 20) ? f.size : ((size_t)1 << 20);
				std::vector<char> head(hn + 64, 0);
				if (f.readAt(head.data(), hn, 0) != hn) return fallback("short read"), -1;
				const char *nl = (const char *)memchr(head.data(), '\n', hn);
				if (nl) {
					ChromTable tmp;
					BedParser hp; hp.needStrand = stranded; hp.chroms = &tmp;
					auto nosink = [](LineRec &, int) -> bool { return true; };
					if (!textPrefixOk(head.data(), hn) || !hp.block(head.data(), nl + 1, nosink)) return -1;
					fileFmt[fi] = hp.fmt;
					char probe[65536 + 64];
					for (int k = 1; k < cfg.threads; k++) {
						size_t pos = f.size / (size_t)cfg.threads * (size_t)k;
						if (pos <= cuts.back()) continue;
						size_t got = f.readAt(probe, 65536, pos - 1);
						const char *q = (const char *)memchr(probe, '\n', got);
						if (!q) continue;
						size_t cut = pos - 1 + (size_t)(q - probe) + 1;
						if (cut > cuts.back() && cut < f.size) cuts.push_back(cut);
					}
				}
			}
			cuts.push_back(f.size);
			for (size_t k = 0; k + 1 < cuts.size(); k++) {
				parts.emplace_back();
				BPart &pt = parts.back();
				pt.file = fi; pt.begin = cuts[k]; pt.end = cuts[k + 1]; pt.first = (k == 0);
			}
		}
		size_t capB = 2;
		for (auto &fb : fbs) capB += fb.size / 5 + 2;
		for (size_t pi = 0; pi < parts.size(); pi++) {
			size_t cap = pi == 0 ? capB : (parts[pi].end - parts[pi].begin) / 5 + 2;
			for (int k = 0; k < nsets; k++)
				if (!parts[pi].sets[k].S.alloc(cap) || !parts[pi].sets[k].E.alloc(cap)) return fallback("allocation"), -1;
		}
	}
	auto parsePart = [&](BPart &pt) {
		const InFile &f = fbs[pt.file];
		pt.trace.enabled = cfg.sorted;
		BedParser bp;
		bp.needStrand = stranded; bp.chroms = &pt.chroms;
		if (!pt.first) { bp.first = false; bp.fmt = fileFmt[pt.file]; }
		int curC[2] = {-1, -1};
		/* per local chromosome: query union (or none) and a forward cursor */
		struct QRef { const QUnion *u; size_t cur; pos_t lastS; };
		std::vector<QRef> qref;
		int qc = -1;
		QRef *q = nullptr;
		/* The record sink is instantiated per (sorted, stranded) combination so
		 * that the per-record flag tests compile away. */
		auto makeSink = [&](auto SORTED, auto STRANDED) {
		return [&, SORTED, STRANDED](LineRec &r, int cid) -> bool {
			constexpr bool kSorted = decltype(SORTED)::value;
			constexpr bool kStranded = decltype(STRANDED)::value;
			pos_t s = r.start, e = r.end;
			if (s == e) { s--; e++; }
			if constexpr (kSorted) {
				if (!pt.trace.on(cid, r.start, s)) return false;
			} else {
				if (s < 0) return fallback("unsorted: zero-length database record at 0");
				if (e > MAX_UNSORTED_COORD) return fallback("unsorted: coordinate beyond bin range");
			}
			if (cid != qc) {
				while ((int)qref.size() <= cid) {
					const std::string &nm = pt.chroms.names[qref.size()];
					auto it = chroms.ids.find(nm); /* read-only during the parallel parse */
					const QUnion *u = (it != chroms.ids.end() && (size_t)it->second < nQueryChroms && !qunion[(size_t)it->second].s.empty()) ? &qunion[(size_t)it->second] : nullptr;
					qref.push_back(QRef{u, 0, INT64_MIN});
				}
				qc = cid; q = &qref[(size_t)cid];
			}
			if (!q->u) return true;
			{
				const QUnion &u = *q->u;
				if (s < q->lastS) q->cur = (size_t)(std::upper_bound(u.e.begin(), u.e.end(), s) - u.e.begin());
				q->lastS = s;
				size_t k = q->cur;
				while (k < u.e.size() && u.e[k] <= s) k++;
				q->cur = k;
				if (k >= u.e.size() || u.s[k] >= e) return true; /* overlaps no query interval */
			}
			int k = 0;
			if constexpr (kStranded) {
				if (r.strand == STR_UNKNOWN) return true; /* can never satisfy -s/-S */
				k = r.strand == STR_REV ? 1 : 0;
			}
			FlatSet &fs = pt.sets[k];
			if (cid != curC[k]) {
				if (!fs.runs.empty()) fs.runs.back().e = fs.S.n;
				fs.runs.push_back(Run{cid, fs.S.n, fs.S.n});
				curC[k] = cid;
			}
			if (e > (pos_t)INT32_MAX) return fallback("database coordinate beyond 32-bit storage");
			fs.S.p[fs.S.n++] = (rec_t)s;
			fs.E.p[fs.E.n++] = (rec_t)e;
			return true;
		};
		};
		auto doParse = [&](auto &sink) {
		/* Parse the range straight from a MAP_POPULATE mapping of the file with
		 * software prefetch in the scanner (measured on the target VM: mapped scan
		 * with prefetch 0.09 s vs pread copy + scan 0.21 s for 600 MB; without
		 * prefetch the mapped scan is slowest, 0.28 s).  Blocks end at newlines and
		 * stay below 1 GiB (32-bit offsets); a final line without a newline is
		 * parsed from a small private copy. */
		const char *base = fileMap[pt.file];
		if (pt.first && !textPrefixOk(base, f.size)) return;
		size_t endc = pt.end;
		if (pt.end == f.size && base[f.size - 1] != '\n') {
			const char *lastNL = (const char *)memrchr(base + pt.begin, '\n', pt.end - pt.begin);
			endc = lastNL ? (size_t)(lastNL + 1 - base) : pt.begin;
		}
		size_t pos = pt.begin;
		while (pos < endc) {
			size_t lim = endc - pos > ((size_t)1 << 30) ? pos + ((size_t)1 << 30) : endc;
			if (lim < endc) {
				const char *q = (const char *)memrchr(base + pos, '\n', lim - pos);
				if (!q) { fallback("line longer than 1 GiB"); return; }
				lim = (size_t)(q + 1 - base);
			}
			if (!bp.block(base + pos, base + lim, sink)) return;
			pos = lim;
		}
		if (endc < pt.end) {
			size_t tl = pt.end - endc;
			std::vector<char> tail(tl + 64, 0);
			memcpy(tail.data(), base + endc, tl);
			tail[tl] = '\n';
			if (!bp.block(tail.data(), tail.data() + tl + 1, sink)) return;
		}
		if (bp.first) { fallback("no records"); return; }
		for (int k = 0; k < nsets; k++) if (!pt.sets[k].runs.empty()) pt.sets[k].runs.back().e = pt.sets[k].S.n;
		pt.fmt = bp.fmt;
		pt.ok = true;
		};
		if (cfg.sorted) {
			if (stranded) { auto sk = makeSink(std::true_type{}, std::true_type{}); doParse(sk); }
			else { auto sk = makeSink(std::true_type{}, std::false_type{}); doParse(sk); }
		} else {
			if (stranded) { auto sk = makeSink(std::false_type{}, std::true_type{}); doParse(sk); }
			else { auto sk = makeSink(std::false_type{}, std::false_type{}); doParse(sk); }
		}
	};
	{
		std::atomic<size_t> nextPart(0);
		auto worker = [&]() {
			for (;;) {
				size_t i = nextPart.fetch_add(1);
				if (i >= parts.size()) break;
				parsePart(parts[i]);
			}
		};
		int nt = cfg.threads < (int)parts.size() ? cfg.threads : (int)parts.size();
		std::vector<std::thread> th;
		for (int t = 1; t < nt; t++) th.emplace_back(worker);
		worker();
		for (auto &x : th) x.join();
	}
	for (auto &pt : parts) if (!pt.ok) return -1;
	/* merge: global chromosome ids, per-file sort validation */
	std::vector<std::vector<int>> remap(parts.size());
	for (size_t pi = 0; pi < parts.size(); pi++) {
		BPart &pt = parts[pi];
		for (auto &nm : pt.chroms.names) remap[pi].push_back(chroms.get(nm.data(), nm.size()));
	}
	if (cfg.sorted) {
		for (size_t fi = 0; fi < fbs.size(); fi++) {
			std::vector<RunInfo> runs;
			for (size_t pi = 0; pi < parts.size(); pi++) {
				if (parts[pi].file != fi) continue;
				for (RunInfo r : parts[pi].trace.runs) { r.cid = remap[pi][r.cid]; runs.push_back(r); }
			}
			if (!validateSorted(runs, chroms, cfg.hasG ? &genome : nullptr)) return -1;
		}
	}
	/* groups indexed by chrom (unstranded) or 2*chrom + (strand==REV); each run of
	 * each part becomes one segment, in file order (no copies) */
	std::vector<TGroup> tracks(chroms.names.size() * (size_t)nsets);
	for (size_t pi = 0; pi < parts.size(); pi++)
		for (int k = 0; k < nsets; k++) {
			FlatSet &ps = parts[pi].sets[k];
			for (const Run &r : ps.runs) {
				if (r.e <= r.b) continue;
				const size_t g = (size_t)remap[pi][(size_t)r.cid];
				TGroup &G = tracks[stranded ? 2 * g + (size_t)k : g];
				G.seg.emplace_back();
				Track &t = G.seg.back();
				t.S = ps.S.p + r.b; t.E = ps.E.p + r.b; t.n = r.e - r.b;
				G.n += t.n;
			}
		}
	prof("parse B");


	/* ---------------- compute + output (no fallback past this point) ---------------- */
	const Mode mode = cfg.mode;
	const bool wantUnion = (mode == M_DEFAULT) && !fraction;
	const bool wantPrefix = (mode == M_MEAN) && !fraction;

	auto trackIndex = [&](const ARec &a) -> size_t {
		if (stranded) {
			if (a.strand == STR_UNKNOWN) return (size_t)-1;
			int want = cfg.same ? a.strand : (a.strand == STR_FWD ? STR_REV : STR_FWD);
			return 2 * (size_t)a.chrom + (want == STR_REV ? 1 : 0);
		}
		return (size_t)a.chrom;
	};
	/* prepare every track a query needs before any (possibly parallel) work */
	uint64_t totalLen = 0;
	{
		std::vector<pos_t> maxALen(tracks.size(), -1);
		for (const ARec &a : A) {
			totalLen += (uint64_t)(a.e - a.s);
			size_t ti = trackIndex(a);
			if (ti < tracks.size() && tracks[ti].n && a.e - a.s > maxALen[ti]) maxALen[ti] = a.e - a.s;
		}
		for (size_t ti = 0; ti < tracks.size(); ti++) {
			if (maxALen[ti] < 0) continue;
			tracks[ti].prepare(fraction, wantUnion, wantPrefix, mode == M_HIST || mode == M_PERBASE, maxALen[ti], cfg.threads);
		}
	}
	prof("prepare tracks");
	auto trackFor = [&](const ARec &a) -> const TGroup * {
		size_t ti = trackIndex(a);
		if (ti >= tracks.size() || tracks[ti].n == 0) return nullptr;
		return &tracks[ti];
	};

	struct Worker {
		std::vector<char> tmpl, blk;
		Hist allHist;
		std::vector<uint64_t> hdense, htouched;
		std::map<uint64_t, uint64_t> hsparse;
		std::vector<std::pair<pos_t, pos_t>> hits;
		std::vector<std::pair<pos_t, int>> events;
		std::vector<int32_t> diff;
	};

	auto processA = [&](const ARec &a, auto &out, Worker &w) {
		const uint64_t L = (uint64_t)(a.e - a.s);
		const TGroup *t = trackFor(a);

		/* Hit set: either the whole track (unfiltered) or a filtered list. */
		uint64_t hitCount = 0;
		bool filtered = fraction && t;
		if (filtered) {
			w.hits.clear();
			const pos_t as = a.s, ae = a.e;
			const float aLen = (float)(ae - as);
			auto cb = [&](pos_t s, pos_t e) {
				pos_t ms = s > as ? s : as;
				pos_t me = e < ae ? e : ae;
				pos_t ov = me - ms;
				float oa = (float)ov / aLen;
				float ob = (float)ov / (float)(e - s);
				bool okA = oa >= cfg.fa, okB = ob >= cfg.fb;
				bool ok = cfg.either ? (okA || okB) : (okA && okB);
				if (ok) w.hits.push_back(std::make_pair(ms, me));
			};
			t->overlaps(as, ae, cb);
			hitCount = w.hits.size();
		} else if (t) {
			hitCount = t->count(a.s, a.e);
		}

		/* depth runs over A, for hist / per-base */
		auto forRuns = [&](auto &&cb) {
			if (!t || (filtered && w.hits.empty()) || (!filtered && hitCount == 0)) { cb(a.s, a.e, (uint64_t)0); return; }
			if (!filtered) { t->runs(a.s, a.e, cb, w.diff); return; }
			w.events.clear();
			for (auto &h : w.hits) { w.events.push_back(std::make_pair(h.first, 1)); w.events.push_back(std::make_pair(h.second, -1)); }
			std::sort(w.events.begin(), w.events.end());
			pos_t cur = a.s;
			int64_t depth = 0;
			size_t k = 0;
			while (k < w.events.size()) {
				pos_t x = w.events[k].first;
				int64_t nd = depth;
				while (k < w.events.size() && w.events[k].first == x) { nd += w.events[k].second; k++; }
				if (nd != depth) {
					if (x > cur) cb(cur, x, (uint64_t)depth);
					cur = x; depth = nd;
				}
			}
			if (cur < a.e) cb(cur, a.e, (uint64_t)depth);
		};

		char *p;
		switch (mode) {
		case M_COUNTS: {
			p = out.reserve(a.len + 32);
			memcpy(p, a.line, a.len); p += a.len; *p++ = '\t';
			p = putU64(p, hitCount); *p++ = '\n';
			out.commit(p);
			break;
		}
		case M_DEFAULT: {
			uint64_t cov = 0;
			if (t && hitCount) {
				if (!filtered) cov = t->covered(a.s, a.e);
				else {
					std::sort(w.hits.begin(), w.hits.end());
					pos_t cs = w.hits[0].first, ce = w.hits[0].second;
					for (size_t k = 1; k < w.hits.size(); k++) {
						if (w.hits[k].first <= ce) { if (w.hits[k].second > ce) ce = w.hits[k].second; }
						else { cov += (uint64_t)(ce - cs); cs = w.hits[k].first; ce = w.hits[k].second; }
					}
					cov += (uint64_t)(ce - cs);
				}
			}
			float frac = (float)cov / (float)L;
			p = out.reserve(a.len + 96);
			memcpy(p, a.line, a.len); p += a.len; *p++ = '\t';
			p = putU64(p, hitCount); *p++ = '\t';
			p = putU64(p, cov); *p++ = '\t';
			p = putU64(p, L); *p++ = '\t';
			p = putF7(p, frac); *p++ = '\n';
			out.commit(p);
			break;
		}
		case M_MEAN: {
			uint64_t sum = 0;
			if (t && hitCount) {
				if (!filtered) sum = t->depthSum(a.s, a.e);
				else for (auto &h : w.hits) sum += (uint64_t)(h.second - h.first);
			}
			float mean = (float)sum / (float)L;
			p = out.reserve(a.len + 600);
			memcpy(p, a.line, a.len); p += a.len; *p++ = '\t';
			p = putF7(p, mean); *p++ = '\n';
			out.commit(p);
			break;
		}
		case M_HIST: {
			w.htouched.clear();
			w.hsparse.clear();
			auto cb = [&](pos_t rs, pos_t re, uint64_t d) {
				uint64_t n = (uint64_t)(re - rs);
				if (d < 65536) {
					if (d >= w.hdense.size()) w.hdense.resize(d + 1, 0);
					if (w.hdense[d] == 0) w.htouched.push_back(d);
					w.hdense[d] += n;
				} else w.hsparse[d] += n;
			};
			forRuns(cb);
			std::sort(w.htouched.begin(), w.htouched.end());
			auto emit = [&](uint64_t d, uint64_t n) {
				w.allHist.add(d, n);
				float fr = (float)n / (float)L;
				char *q = out.reserve(a.len + 96);
				memcpy(q, a.line, a.len); q += a.len; *q++ = '\t';
				q = putU64(q, d); *q++ = '\t';
				q = putU64(q, n); *q++ = '\t';
				q = putU64(q, L); *q++ = '\t';
				q = putF7(q, fr); *q++ = '\n';
				out.commit(q);
			};
			for (uint64_t d : w.htouched) { emit(d, w.hdense[d]); w.hdense[d] = 0; }
			for (auto &kv : w.hsparse) emit(kv.first, kv.second);
			break;
		}
		case M_PERBASE: {
			/* line = A line + '\t' + 1-based position + '\t' + depth + '\n'.  The
			 * position is an ASCII counter incremented in place.  Within a run, whole
			 * decades (positions x0..x9) are emitted as one prebuilt 10-line block in
			 * which only the changed higher digits are patched between decades. */
			const size_t pl = (size_t)a.len + 1;
			if (w.tmpl.size() < pl + 64) w.tmpl.resize(pl + 64);
			char *T = w.tmpl.data();
			memcpy(T, a.line, a.len);
			T[a.len] = '\t';
			char *dig = T + pl;
			size_t nd = 1;
			dig[0] = '1';
			auto incr = [&](size_t from) -> size_t { /* increment digits [0, from) ; returns #changed low digits, or SIZE_MAX if length grew */
				size_t k = from;
				while (k > 0 && dig[k - 1] == '9') { dig[k - 1] = '0'; k--; }
				if (k > 0) { dig[k - 1]++; return from - k + 1; }
				/* all nines: 99..9 -> 100..0 */
				dig[0] = '1';
				for (size_t j = 1; j <= nd; j++) dig[j] = '0';
				nd++;
				return (size_t)-1;
			};
			auto cb = [&](pos_t rs, pos_t re, uint64_t d) {
				char dstr[48];
				dstr[0] = '\t';
				char *de = putU64(dstr + 1, d);
				*de++ = '\n';
				const size_t dl = (size_t)(de - dstr);
				uint64_t left = (uint64_t)(re - rs);
				bool blockValid = false;
				size_t bll = 0;
				while (left) {
					if (nd >= 2 && dig[nd - 1] == '0' && left >= 10) {
						const size_t ll = pl + nd + dl;
						if (!blockValid || bll != ll) {
							if (w.blk.size() < 10 * ll + 64) w.blk.resize(10 * ll + 64);
							char *b = w.blk.data();
							for (int j = 0; j < 10; j++) {
								memcpy(b + j * ll, T, pl + nd);
								b[j * ll + pl + nd - 1] = (char)('0' + j);
								memcpy(b + j * ll + pl + nd, dstr, dl);
							}
							blockValid = true; bll = ll;
						}
						char *q = out.reserve(10 * ll + 64);
						memcpy(q, w.blk.data(), 10 * ll);
						out.commit(q + 10 * ll);
						left -= 10;
						/* advance by ten: increment the tens part; the last digit stays '0' */
						size_t ch = incr(nd - 1);
						if (ch == (size_t)-1) blockValid = false;
						else {
							char *b = w.blk.data();
							const size_t lo = pl + nd - 1 - ch; /* first changed digit offset */
							if (ch == 1) {
								const char c = T[lo];
								for (int j = 0; j < 10; j++) b[j * ll + lo] = c;
							} else for (int j = 0; j < 10; j++) memcpy(b + j * ll + lo, T + lo, ch);
						}
						continue;
					}
					const size_t ll = pl + nd;
					char *q = out.reserve(ll + 64);
					if (ll >= 16 && ll <= 32) {
						__m128i v0 = _mm_loadu_si128((const __m128i *)T);
						__m128i v1 = _mm_loadu_si128((const __m128i *)(T + ll - 16));
						_mm_storeu_si128((__m128i *)q, v0);
						_mm_storeu_si128((__m128i *)(q + ll - 16), v1);
					} else memcpy(q, T, ll);
					_mm_storeu_si128((__m128i *)(q + ll), _mm_loadu_si128((const __m128i *)dstr));
					if (dl > 16) memcpy(q + ll + 16, dstr + 16, dl - 16);
					out.commit(q + ll + dl);
					left--;
					if (incr(nd) == (size_t)-1) blockValid = false;
					else if (dig[nd - 1] == '0') blockValid = false; /* higher digits changed */
				}
			};
			forRuns(cb);
			break;
		}
		}
		};

	const int nth = cfg.threads;
	std::vector<Worker> workers((size_t)(nth > 1 ? nth : 1));
	/* 512 KiB: the formatted bytes are still in L2 (2 MiB per core here) when write() copies them */
	Out out((size_t)512 << 10);
	/* Blocks of consecutive queries (~1 MiB of estimated output) are formatted
	 * by all threads into private buffers; the main thread writes them in order
	 * and formats blocks itself whenever the next block is not ready.  Outputs
	 * estimated below one block are formatted by the main thread alone. */
	std::vector<std::pair<size_t, size_t>> blocks;
	if (nth > 1 && A.size() >= 2) {
		const double TARGET = 1048576.0;
		double acc = 0;
		size_t b0 = 0;
		for (size_t i = 0; i < A.size(); i++) {
			const ARec &a = A[i];
			double L = (double)(a.e - a.s);
			double est = mode == M_PERBASE ? L * ((double)a.len + 14.0) : mode == M_HIST ? ((double)a.len + 40.0) * 4.0 : (double)a.len + 40.0;
			acc += est;
			if (acc >= TARGET) { blocks.push_back(std::make_pair(b0, i + 1)); b0 = i + 1; acc = 0; }
		}
		if (b0 < A.size()) blocks.push_back(std::make_pair(b0, A.size()));
	}
	if (blocks.size() < 2) {
		for (const ARec &a : A) processA(a, out, workers[0]);
	} else {
		const size_t nb = blocks.size();
		std::vector<Buf> bufs(nb);
		std::unique_ptr<std::atomic<int>[]> done(new std::atomic<int>[nb]);
		for (size_t i = 0; i < nb; i++) done[i].store(0);
		std::atomic<size_t> nextBlock(0), written(0);
		std::mutex mu;
		std::condition_variable cv;
		const size_t LIMIT = (size_t)nth * 6;
		auto runBlock = [&](size_t bi, Worker &w) {
			for (size_t i = blocks[bi].first; i < blocks[bi].second; i++) processA(A[i], bufs[bi], w);
			{ std::lock_guard<std::mutex> lk(mu); done[bi].store(1); }
			cv.notify_all();
		};
		auto workerLoop = [&](int wi) {
			for (;;) {
				size_t bi;
				{
					std::unique_lock<std::mutex> lk(mu);
					cv.wait(lk, [&] { return nextBlock.load() >= nb || nextBlock.load() < written.load() + LIMIT; });
					bi = nextBlock.load();
					if (bi >= nb) return;
					nextBlock.store(bi + 1);
				}
				runBlock(bi, workers[(size_t)wi]);
			}
		};
		std::vector<std::thread> th;
		for (int t = 1; t < nth; t++) th.emplace_back(workerLoop, t);
		out.flush();
		for (size_t wb = 0; wb < nb; wb++) {
			for (;;) {
				if (done[wb].load()) break;
				size_t bi = (size_t)-1;
				{
					std::unique_lock<std::mutex> lk(mu);
					if (done[wb].load()) break;
					if (nextBlock.load() < nb && nextBlock.load() < written.load() + LIMIT) { bi = nextBlock.load(); nextBlock.store(bi + 1); }
					else { cv.wait(lk, [&] { return done[wb].load() != 0; }); break; }
				}
				runBlock(bi, workers[0]);
			}
			out.writeAll(bufs[wb].p, bufs[wb].len);
			free(bufs[wb].p); bufs[wb].p = nullptr;
			{ std::lock_guard<std::mutex> lk(mu); written.store(wb + 1); }
			cv.notify_all();
		}
		for (auto &x : th) x.join();
		for (size_t wi = 1; wi < workers.size(); wi++) {
			Hist &h = workers[wi].allHist;
			for (size_t d = 0; d < h.dense.size(); d++) if (h.dense[d]) workers[0].allHist.add(d, h.dense[d]);
			for (auto &kv : h.sparse) workers[0].allHist.add(kv.first, kv.second);
		}
	}
	if (mode == M_HIST) {
		Hist &allHist = workers[0].allHist;
		auto emitAll = [&](uint64_t d, uint64_t n) {
			float pct = (float)n / (float)totalLen;
			char *q = out.reserve(128);
			memcpy(q, "all\t", 4); q += 4;
			q = putU64(q, d); *q++ = '\t';
			q = putU64(q, n); *q++ = '\t';
			q = putU64(q, totalLen); *q++ = '\t';
			q = putF7(q, pct); *q++ = '\n';
			out.commit(q);
		};
		for (size_t d = 0; d < allHist.dense.size(); d++) if (allHist.dense[d]) emitAll(d, allHist.dense[d]);
		for (auto &kv : allHist.sparse) emitAll(kv.first, kv.second);
	}
	out.flush();
	prof("compute+output");
	return 0;
}

} /* namespace swiftcov */

extern char **environ;

/* Returns the exit status of a completed fast-path run, or -1 if the upstream
 * implementation must handle the invocation in this process.  -1 is returned
 * only before any heap activity; once the fast path has done work, fallback
 * re-executes the binary (marking the child with SWIFTCOV_FALLBACK=1) because
 * upstream reads the uninitialised ContextCoverage::_mean member, whose value
 * depends on prior heap state (observed: spurious "mutually exclusive" errors
 * after an in-process fallback). */
int swiftcov_fast_coverage(int argc, char **argv) {
	if (getenv("SWIFTCOV_FALLBACK") || getenv("SWIFTCOV_NOFAST")) return -1;
	for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "--cram-ref")) return -1;
	int rc = swiftcov::run(argc, argv);
	if (rc != -1) return rc;
	size_t n = 0;
	while (environ && environ[n]) n++;
	char **env = (char **)malloc((n + 2) * sizeof(char *));
	if (env) {
		for (size_t i = 0; i < n; i++) env[i] = environ[i];
		env[n] = (char *)"SWIFTCOV_FALLBACK=1";
		env[n + 1] = nullptr;
		execve("/proc/self/exe", argv, env);
	}
	if (swiftcov::g_debug) fprintf(stderr, "swiftcov: re-exec failed; in-process fallback\n");
	return -1;
}
