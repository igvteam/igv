package org.igv.figwig;

import org.json.JSONObject;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * The regional summary index (header optional bit 1; "The regional summary index" in the figwig format
 * specification): where a group's summary and zoom streams were restarted, so that a view of one region
 * inflates only the windows of bins it covers instead of a whole level.  Its level catalogs also give
 * the levels' places in a group's record without reading the zoom section.
 *
 * Every check the specification asks of a reader is made before a value is used.  A stale index (its
 * copy of the header differs from the file's) is ignored; one that fails a check makes the file corrupt.
 */
final class FigWigSummaryIndex {

    static final String KEY = "org.figwig.summary-restarts";
    private static final int VERSION = 2, ROOT_BYTES = 256, PAGE_ENTRIES = 102, PAGE_STRIDE = 4096;
    private static final int TARGET_BYTES = 1 << 15, MAX_WINDOW_BINS = 2048, MAX_LEVELS = 13;
    private static final long INDEX_BYTES = 1 << 24;
    // Decoded windows kept for panning back and forth, in bytes.
    private static final long WINDOW_CACHE_BYTES = Long.getLong("figwig.summaryWindowBytes", 64L << 20);

    /** The file's bytes, from FigWigFile: one range, or several at once (concurrent requests over HTTP). */
    interface Source {
        byte[] readAt(long offset, int n) throws IOException;

        byte[][] readEach(long[] offsets, int[] lengths) throws IOException;
    }

    /** One level's entry in a group's catalog. */
    static final class LevelEntry {
        long sourceOffset, indexOffset, sourceCrc;
        int sourceRows, sourceBytes, windowBins;
        long bins;            // the level's bins over the genome
        volatile boolean streamChecked;

        boolean restarted() {
            return windowBins != 0;
        }

        long windows() {
            return restarted() ? (bins + windowBins - 1) / windowBins : 0;
        }
    }

    /** A window's rows: bin numbers and the sums of the group's tracks per bin. */
    static final class Window {
        final int[] bins;
        final byte[] sums;

        Window(int[] bins, byte[] sums) {
            this.bins = bins;
            this.sums = sums;
        }
    }

    private final Source file;
    private final byte[] header;
    private final long root, chromsOffset;
    private final int nGroups;
    private final boolean hasZooms;
    private volatile Boolean usable;
    private byte[] directory;   // the groups' locators, read with the root
    private byte[] catalogBytes;   // every group's catalog when they lie together, read with the first
    private long catalogFrom;
    private final Object[] groupLocks;
    private final LevelEntry[][] catalogs;
    private final LinkedHashMap<Long, Window> windows = new LinkedHashMap<>(256, 0.75f, true);
    private long windowBytes = 0;

    private FigWigSummaryIndex(Source file, byte[] header, long root, long chromsOffset, int nGroups, boolean hasZooms) {
        this.file = file;
        this.header = header;
        this.root = root;
        this.chromsOffset = chromsOffset;
        this.nGroups = nGroups;
        this.hasZooms = hasZooms;
        this.catalogs = new LevelEntry[nGroups][];
        this.groupLocks = new Object[nGroups];
        for (int g = 0; g < nGroups; g++) {
            groupLocks[g] = new Object();
        }
    }

    /**
     * The index the header and metadata announce, or null: without bit 1 or the key, or with a locator
     * of another version.  A locator of this version that fails its checks makes the file corrupt.
     */
    static FigWigSummaryIndex of(Source file, byte[] header, long optional, JSONObject metadata, long chromsOffset,
                                 int nGroups) {
        if ((optional & 2) == 0 || !metadata.has(KEY) || !(metadata.get(KEY) instanceof JSONObject)) {
            return null;
        }
        JSONObject locator = metadata.getJSONObject(KEY);
        if (locator.optInt("version", -1) != VERSION) {
            return null;
        }
        long offset = locator.optLong("offset", -1), length = locator.optLong("length", -1);
        check(length == ROOT_BYTES && offset >= 128 && offset + ROOT_BYTES <= chromsOffset
                && chromsOffset - offset <= INDEX_BYTES, "the summary index's locator is out of place");
        return new FigWigSummaryIndex(file, header, offset, chromsOffset, nGroups, (optional & 1) != 0);
    }

    /** Whether the index holds for this file: its root is read and checked the first time. */
    boolean usable() throws IOException {
        Boolean u = usable;
        if (u != null) {
            return u;
        }
        synchronized (this) {
            if (usable == null) {
                check(root + ROOT_BYTES + 32L * nGroups <= chromsOffset, "the summary index's directory passes its region");
                byte[] rd = file.readAt(root, ROOT_BYTES + 32 * nGroups);
                byte[] r = Arrays.copyOf(rd, ROOT_BYTES);
                directory = Arrays.copyOfRange(rd, ROOT_BYTES, rd.length);
                ByteBuffer b = le(r);
                check(Arrays.equals(Arrays.copyOf(r, 8), "FWRSUM1\0".getBytes(java.nio.charset.StandardCharsets.US_ASCII))
                        && (b.getShort(8) & 0xffff) == VERSION && (b.getShort(10) & 0xffff) == ROOT_BYTES,
                        "the summary index's root is not one");
                check(crcWithout(r, 0, ROOT_BYTES, 12) == (b.getInt(12) & 0xffffffffL), "the summary index's root does not match its CRC-32");
                check(zeros(r, 36, 64) && zeros(r, 192, 256), "the summary index's root has reserved bytes set");
                if (!Arrays.equals(Arrays.copyOfRange(r, 64, 192), header)) {
                    usable = false;   // stale: the file changed after the index was written
                } else {
                    check(b.getLong(16) == root + ROOT_BYTES && (b.getInt(24) & 0xffffffffL) == nGroups
                            && b.getInt(28) == PAGE_ENTRIES && b.getInt(32) == TARGET_BYTES,
                            "the summary index's root does not match the file");
                    usable = true;
                }
            }
            return usable;
        }
    }

    /**
     * Group g's level catalog, checked against its record: the summary's header at summaryAt (n, bytes,
     * CRC-32 as the record holds them), the zoom section from zoomAt to the record's end.
     */
    LevelEntry[] catalog(int g, int nTracks, long recordStart, long recordLength, int summaryAt, int summaryN,
                         int summaryBytes, long summaryCrc, int zoomAt, long[] levelBins) throws IOException {
        LevelEntry[] c = catalogs[g];
        if (c != null) {
            return c;
        }
        synchronized (groupLocks[g]) {
            if (catalogs[g] != null) {
                return catalogs[g];
            }
            byte[] loc = Arrays.copyOfRange(directory, 32 * g, 32 * g + 32);
            ByteBuffer lb = le(loc);
            check(crcWithout(loc, 0, 32, 24) == (lb.getInt(24) & 0xffffffffL) && lb.getInt(28) == 0,
                    "group " + g + "'s summary index locator does not match its CRC-32");
            long levelsOffset = lb.getLong(0);
            int nLevels = lb.getInt(8);
            check((lb.getInt(16) & 0xffffffffL) == g && (lb.getInt(20) & 0xffffffffL) == nTracks,
                    "group " + g + "'s summary index locator names another group");
            check(nLevels >= 1 && nLevels <= MAX_LEVELS && (nLevels == 1 || hasZooms) && nLevels <= levelBins.length,
                    "group " + g + "'s summary index has " + nLevels + " levels");
            inData(levelsOffset, 32L * nLevels, "group " + g + "'s level catalog");
            byte[] cat = catalogBytes(levelsOffset, 32 * nLevels);
            check(crc(cat, 0, cat.length) == (lb.getInt(12) & 0xffffffffL), "group " + g + "'s level catalog does not match its CRC-32");
            ByteBuffer cb = le(cat);
            LevelEntry[] out = new LevelEntry[nLevels];
            List<long[]> spans = new ArrayList<>();
            spans.add(new long[]{levelsOffset, levelsOffset + 32L * nLevels});
            long expectAt = recordStart + summaryAt;
            for (int k = 0; k < nLevels; k++) {
                LevelEntry e = new LevelEntry();
                int p = 32 * k;
                e.sourceOffset = cb.getLong(p);
                e.indexOffset = cb.getLong(p + 8);
                e.sourceRows = cb.getInt(p + 16);
                e.sourceBytes = cb.getInt(p + 20);
                e.sourceCrc = cb.getInt(p + 24) & 0xffffffffL;
                e.windowBins = cb.getInt(p + 28);
                e.bins = levelBins[k];
                String what = "group " + g + "'s summary index level " + k;
                if (k == 0) {
                    check(e.sourceRows == summaryN && e.sourceBytes == summaryBytes && e.sourceCrc == summaryCrc,
                            what + " does not match the record's summary");
                } else if (k == 1) {
                    expectAt = recordStart + zoomAt + 4 + 16;   // the zoom count, then level 1's header
                } else {
                    expectAt += 16;                             // level k's header
                }
                check(e.sourceOffset == expectAt, what + " is not where the record puts it");
                check(e.sourceRows >= 0 && e.sourceBytes >= 0 && e.sourceRows <= e.bins, what + " lists more bins than it has");
                check(e.sourceBytes != 0 || e.sourceRows == 0, what + " lists rows but holds no stream");
                check(e.sourceRows != 0 || e.sourceBytes == 0 || e.sourceBytes == 8, what + " is empty but holds a stream");
                if (e.indexOffset == 0 && e.windowBins == 0) {
                    // ordinary
                } else {
                    check(e.sourceRows > 0 && e.windowBins > 0 && Integer.bitCount(e.windowBins) == 1
                                    && e.windowBins <= MAX_WINDOW_BINS
                                    && (long) e.windowBins * (4 + (long) (k == 0 ? 4 : 8) * nTracks) <= TARGET_BYTES,
                            what + "'s windows are out of bounds");
                    long w = e.windows(), pages = (w + PAGE_ENTRIES - 1) / PAGE_ENTRIES;
                    long last = 16 + 40 * (w - PAGE_ENTRIES * (pages - 1));
                    long span = PAGE_STRIDE * (pages - 1) + last;
                    inData(e.indexOffset, span, what + "'s pages");
                    spans.add(new long[]{e.indexOffset, e.indexOffset + span});
                }
                expectAt = e.sourceOffset + e.sourceBytes;
                out[k] = e;
            }
            long end = nLevels > 1 ? expectAt : recordStart + summaryAt + summaryBytes;
            check(nLevels == 1 ? end <= recordStart + recordLength : end == recordStart + recordLength,
                    "group " + g + "'s summary index levels do not end where the record does");
            spans.sort(Comparator.comparingLong(s -> s[0]));
            for (int i = 1; i < spans.size(); i++) {
                check(spans.get(i)[0] >= spans.get(i - 1)[1], "group " + g + "'s summary index structures overlap");
            }
            catalogs[g] = out;
            return out;
        }
    }

    /**
     * A catalog's bytes: the catalogs of all groups are read in one request the first time when they lie
     * within 1 MiB of one another, as a writer lays them; otherwise each alone.
     */
    private byte[] catalogBytes(long offset, int length) throws IOException {
        synchronized (this) {
            if (catalogBytes == null) {
                long lo = Long.MAX_VALUE, hi = 0;
                ByteBuffer d = le(directory);
                for (int g = 0; g < nGroups; g++) {
                    long at = d.getLong(32 * g), n = d.getInt(32 * g + 8) & 0xffffffffL;
                    if (n < 1 || n > MAX_LEVELS || at < 128 || at + 32 * n > chromsOffset) {
                        lo = Long.MAX_VALUE;   // a bad locator: checked when its group is read
                        break;
                    }
                    lo = Math.min(lo, at);
                    hi = Math.max(hi, at + 32 * n);
                }
                if (lo != Long.MAX_VALUE && hi - lo <= 1 << 20) {
                    catalogFrom = lo;
                    catalogBytes = file.readAt(lo, (int) (hi - lo));
                } else {
                    catalogBytes = new byte[0];
                }
            }
        }
        if (offset >= catalogFrom && offset + length <= catalogFrom + catalogBytes.length) {
            int at = (int) (offset - catalogFrom);
            return Arrays.copyOfRange(catalogBytes, at, at + length);
        }
        return file.readAt(offset, length);
    }

    /**
     * Before reading level k's stream or fragments the first time: the record's header for the level
     * (the summary's is checked by catalog) and the stream's zlib header.  For level 1 the record's zoom
     * count, just before its header, is read too.  The bytes to read for it: {offset, length}.
     */
    static long[] streamHead(LevelEntry e, int k) {
        int before = k == 0 ? 0 : k == 1 ? 20 : 16;
        return new long[]{e.sourceOffset - before, before + Math.min(2, e.sourceBytes)};
    }

    /** Level k's stream, whole, read with the headers before it and checked; null if it holds none. */
    byte[] stream(int g, LevelEntry[] cat, int k) throws IOException {
        LevelEntry e = cat[k];
        if (e.sourceBytes == 0) {
            return null;
        }
        long[] head = streamHead(e, k);
        int before = (int) (e.sourceOffset - head[0]);
        byte[] b = file.readAt(head[0], before + e.sourceBytes);
        checkStream(g, cat, k, cat.length, b);
        return Arrays.copyOfRange(b, before, b.length);
    }

    private void checkStream(int g, LevelEntry[] cat, int k, int levels, byte[] b) {
        LevelEntry e = cat[k];
        int before = k == 0 ? 0 : k == 1 ? 20 : 16;
        ByteBuffer h = le(b);
        String what = "group " + g + "'s level " + k;
        if (k >= 1) {
            int at = before - 16;
            if (k == 1) {
                check((h.getInt(0) & 0xffffffffL) == levels - 1, "group " + g + "'s zoom count does not match its summary index");
            }
            check(h.getInt(at) == (1 << (2 * k)) && h.getInt(at + 4) == e.sourceRows && h.getInt(at + 8) == e.sourceBytes
                    && (h.getInt(at + 12) & 0xffffffffL) == e.sourceCrc, what + "'s header does not match its summary index");
        }
        if (e.sourceBytes >= 2) {
            int cmf = b[before] & 0xff, flg = b[before + 1] & 0xff;
            check((cmf & 0x0f) == 8 && (flg & 0x20) == 0 && (cmf * 256 + flg) % 31 == 0
                            && (!e.restarted() || (cmf >> 4) == 7),
                    what + "'s stream does not start as zlib");
        }
        e.streamChecked = true;
    }

    /**
     * The rows of level k of group g in windows [j0, j1], from the cache or the file; w is 4 for the
     * summary, 8 for a zoom level.  Pages are read once for the range, then the bin numbers and the
     * sums of the windows missing from the cache as one read each.
     */
    Window[] windows(int g, LevelEntry[] cat, int k, int nTracks, long j0, long j1) throws IOException {
        Window[] out = cached(g, k, j0, j1);
        if (out != null) {
            return out;
        }
        // One read of a level's windows at a time: the tracks of a group ask for the same windows at once.
        synchronized (cat[k]) {
            return fetch(g, cat, k, nTracks, j0, j1);
        }
    }

    /** Windows [j0, j1] if all are cached, else null. */
    private Window[] cached(int g, int k, long j0, long j1) {
        Window[] out = new Window[(int) (j1 - j0 + 1)];
        synchronized (windows) {
            for (int i = 0; i < out.length; i++) {
                out[i] = windows.get(key(g, k, j0 + i));
                if (out[i] == null) {
                    return null;
                }
            }
        }
        return out;
    }

    private Window[] fetch(int g, LevelEntry[] cat, int k, int nTracks, long j0, long j1) throws IOException {
        LevelEntry e = cat[k];
        int n = (int) (j1 - j0 + 1);
        Window[] out = new Window[n];
        long firstMissing = -1, lastMissing = -1;
        synchronized (windows) {
            for (int i = 0; i < n; i++) {
                out[i] = windows.get(key(g, k, j0 + i));
                if (out[i] == null) {
                    if (firstMissing < 0) firstMissing = j0 + i;
                    lastMissing = j0 + i;
                }
            }
        }
        if (firstMissing < 0) {
            return out;
        }
        int w = k == 0 ? 4 : 8;
        int[][] d = descriptors(g, cat, k, nTracks, firstMissing, lastMissing);
        // Fragments of the windows with rows lie in window order: one read for the ids, one for the sums.
        int a = -1, z = -1;
        for (int i = 0; i < d.length; i++) {
            if (d[i][0] > 0) {
                if (a < 0) a = i;
                z = i;
            }
        }
        byte[] ids = null, sums = null;
        long idsFrom = 0, sumsFrom = 0;
        if (a >= 0) {
            idsFrom = d[a][2];
            sumsFrom = d[a][6];
            byte[][] both = file.readEach(new long[]{e.sourceOffset + idsFrom, e.sourceOffset + sumsFrom},
                    new int[]{(int) (d[z][2] + (long) d[z][3] - idsFrom), (int) (d[z][6] + (long) d[z][7] - sumsFrom)});
            ids = both[0];
            sums = both[1];
        }
        for (int i = 0; i < d.length; i++) {
            long j = firstMissing + i;
            int rows = d[i][0];
            Window win;
            if (rows == 0) {
                win = new Window(new int[0], new byte[0]);
            } else {
                byte[] idBytes = fragment(ids, (int) (d[i][2] - idsFrom), d[i][3], 4 * rows, d[i][4], d[i][5]);
                long sumLen = (long) rows * nTracks * w;
                check(sumLen < Integer.MAX_VALUE, "a summary window is too large");
                byte[] sumBytes = fragment(sums, (int) (d[i][6] - sumsFrom), d[i][7], (int) sumLen, d[i][8], d[i][9]);
                int[] bins = new int[rows];
                ByteBuffer ib = le(idBytes);
                long lo = j * e.windowBins, hi = Math.min((j + 1) * e.windowBins, e.bins);
                long prev = -1;
                for (int r = 0; r < rows; r++) {
                    long bin = ib.getInt(4 * r) & 0xffffffffL;
                    check(bin > prev && bin >= lo && bin < hi, "a summary window's bins do not ascend inside it");
                    bins[r] = (int) bin;
                    prev = bin;
                }
                win = new Window(bins, sumBytes);
            }
            synchronized (windows) {
                windows.put(key(g, k, j), win);
                windowBytes += 64 + 4L * win.bins.length + win.sums.length;
                Iterator<Map.Entry<Long, Window>> it = windows.entrySet().iterator();
                while (windowBytes > WINDOW_CACHE_BYTES && it.hasNext()) {
                    Window old = it.next().getValue();
                    windowBytes -= 64 + 4L * old.bins.length + old.sums.length;
                    it.remove();
                }
            }
            out[(int) (j - j0)] = win;
        }
        return out;
    }

    /** Decoded windows kept now; for tests. */
    int windowsKept() {
        synchronized (windows) {
            return windows.size();
        }
    }

    void clear() {
        synchronized (windows) {
            windows.clear();
            windowBytes = 0;
        }
    }

    /** The ten fields of the descriptors of windows [j0, j1], from their pages, checked. */
    private int[][] descriptors(int g, LevelEntry[] cat, int k, int nTracks, long j0, long j1) throws IOException {
        LevelEntry e = cat[k];
        long total = e.windows();
        check(j0 >= 0 && j1 < total, "a summary window past the level");
        long p0 = j0 / PAGE_ENTRIES, p1 = j1 / PAGE_ENTRIES;
        long lastCount = Math.min(PAGE_ENTRIES, total - PAGE_ENTRIES * p1);
        long pagesAt = e.indexOffset + PAGE_STRIDE * p0;
        int pagesLength = (int) (PAGE_STRIDE * (p1 - p0) + 16 + 40 * lastCount);
        byte[] pages;
        if (e.streamChecked) {
            pages = file.readAt(pagesAt, pagesLength);
        } else {
            // The level's first read: its headers with the pages.
            long[] head = streamHead(e, k);
            byte[][] both = file.readEach(new long[]{pagesAt, head[0]}, new int[]{pagesLength, (int) head[1]});
            pages = both[0];
            checkStream(g, cat, k, cat.length, both[1]);
        }
        int[][] out = new int[(int) (j1 - j0 + 1)][];
        String what = "group " + g + "'s level " + k + " summary page";
        for (long p = p0; p <= p1; p++) {
            int at = (int) (PAGE_STRIDE * (p - p0));
            ByteBuffer b = le(pages);
            long count = Math.min(PAGE_ENTRIES, total - PAGE_ENTRIES * p);
            check(pages[at] == 'F' && pages[at + 1] == 'W' && pages[at + 2] == 'R' && pages[at + 3] == 'P'
                    && (b.getInt(at + 4) & 0xffffffffL) == PAGE_ENTRIES * p && b.getInt(at + 8) == count, what + " is not one");
            check(crcWithout(pages, at, 16 + 40 * (int) count, 12) == (b.getInt(at + 12) & 0xffffffffL), what + " does not match its CRC-32");
            long nextRow = -1, nextIds = -1, nextSums = -1;
            for (int i = 0; i < count; i++) {
                int q = at + 16 + 40 * i;
                int[] f = new int[10];
                for (int x = 0; x < 10; x++) f[x] = b.getInt(q + 4 * x);
                long j = PAGE_ENTRIES * p + i;
                if (f[0] == 0) {
                    for (int x = 1; x < 10; x++) check(f[x] == 0, what + " has a window with no rows but fragments");
                } else {
                    long rows = f[0] & 0xffffffffL, first = f[1] & 0xffffffffL;
                    long winBins = Math.min((j + 1) * e.windowBins, e.bins) - j * e.windowBins;
                    long idsOff = f[2] & 0xffffffffL, idsLen = f[3] & 0xffffffffL;
                    long sumsOff = f[6] & 0xffffffffL, sumsLen = f[7] & 0xffffffffL;
                    long sumRaw = rows * nTracks * (k == 0 ? 4 : 8);
                    check(rows <= e.windowBins && rows <= winBins && first + rows <= e.sourceRows, what + " lists too many rows");
                    check(idsOff >= 2 && idsLen > 0 && idsLen <= 2 * 4 * rows + 64 && idsOff + idsLen <= sumsOff
                                    && sumsLen > 0 && sumsLen <= 2 * sumRaw + 64 && sumsOff + sumsLen <= e.sourceBytes - 4,
                            what + " has fragments out of place");
                    check(nextRow < 0 || (first >= nextRow && idsOff >= nextIds && sumsOff >= nextSums),
                            what + " has windows out of order");
                    nextRow = first + rows;
                    nextIds = idsOff + idsLen;
                    nextSums = sumsOff + sumsLen;
                }
                if (j >= j0 && j <= j1) {
                    out[(int) (j - j0)] = f;
                }
            }
        }
        return out;
    }

    /**
     * One fragment inflated by itself: raw deflate ending in the full-flush marker, followed by an
     * empty final block; exactly {@code raw} bytes out, with both CRC-32 values.
     */
    private static byte[] fragment(byte[] buf, int off, int len, int raw, int crcCompressed, int crcRaw) throws IOException {
        check(off >= 0 && len > 0 && off + (long) len <= buf.length && len <= 2L * raw + 64, "a summary fragment is out of place");
        check(crc(buf, off, len) == (crcCompressed & 0xffffffffL), "a summary fragment does not match its CRC-32");
        check(len >= 4 && buf[off + len - 4] == 0 && buf[off + len - 3] == 0 && (buf[off + len - 2] & 0xff) == 0xff
                && (buf[off + len - 1] & 0xff) == 0xff, "a summary fragment does not end in a full flush");
        byte[] in = Arrays.copyOfRange(buf, off, off + len + 2);
        in[len] = 3;
        in[len + 1] = 0;
        Inflater inf = new Inflater(true);
        try {
            inf.setInput(in);
            byte[] out = new byte[raw];
            int done = 0;
            while (done < raw) {
                int r = inf.inflate(out, done, raw - done);
                if (r == 0 && (inf.finished() || inf.needsInput() || inf.needsDictionary())) {
                    break;
                }
                done += r;
            }
            check(done == raw, "a summary fragment inflates to " + done + " bytes, not " + raw);
            if (!inf.finished()) {
                check(inf.inflate(new byte[1]) == 0 && inf.finished(), "a summary fragment inflates past its size");
            }
            check(inf.getRemaining() == 0, "bytes after a summary fragment");
            check(crc(out, 0, raw) == (crcRaw & 0xffffffffL), "a summary fragment's bytes do not match their CRC-32");
            return out;
        } catch (DataFormatException ex) {
            throw new IOException("figwig: a corrupt summary fragment: " + ex.getMessage(), ex);
        } finally {
            inf.end();
        }
    }

    private void inData(long offset, long length, String what) {
        check(offset >= 128 && length >= 0 && offset + length <= chromsOffset, what + " lies outside the data region");
    }

    private static long key(int g, int k, long j) {
        return ((long) g << 40) | ((long) k << 36) | j;
    }

    private static ByteBuffer le(byte[] b) {
        return ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN);
    }

    private static long crc(byte[] b, int off, int len) {
        CRC32 c = new CRC32();
        c.update(b, off, len);
        return c.getValue();
    }

    private static long crcWithout(byte[] b, int off, int len, int field) {
        byte[] copy = Arrays.copyOfRange(b, off, off + len);
        Arrays.fill(copy, field, field + 4, (byte) 0);
        return crc(copy, 0, len);
    }

    private static boolean zeros(byte[] b, int from, int to) {
        for (int i = from; i < to; i++) {
            if (b[i] != 0) return false;
        }
        return true;
    }

    private static void check(boolean ok, String what) {
        if (!ok) {
            throw new IllegalStateException("figwig: " + what);
        }
    }
}
