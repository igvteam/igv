package org.igv.figwig;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.igv.util.stream.IGVSeekableStreamFactory;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import java.util.function.IntToLongFunction;
import java.util.zip.CRC32;

import static org.igv.figwig.FigWigChunks.*;
import static org.igv.figwig.FigWigPool.*;

/**
 * A reader of figwig files, version 1, as specified in docs/format.rst of the figwig repository
 * (https://github.com/jmschrei/figwig).  Opening reads the header and the tables; a group's index record
 * is read the first time one of its tracks is read (or in the background when a track is created), and a
 * summary or zoom level the first time it is asked for, by windows when the file has a regional summary
 * index ({@link FigWigSummaryIndex}).
 *
 * A read loads the chunks under its range once for every track of the group (a chunk holds all of the
 * group's tracks), caches them under one budget shared by every open file, and decodes the track's
 * values into runs ({@link FigWigChunks}).  Long reads are split over the decoding threads
 * ({@link FigWigPool}); reads of 16 or more tracks of a group over one range are decoded in one pass.
 * Local files are read with positional reads; http(s), s3 and gs paths with one range request per read,
 * the requests of a read sent at once.
 */
public class FigWigFile {

    public static final int MAGIC = 0x57474946;  // "FIGW"

    private static final long CACHE_BYTES = Long.getLong("figwig.cacheMB", 256) << 20;
    private static final boolean POSITIONAL = !Boolean.getBoolean("figwig.noPositional");
    private static final long MAX_SPAN = 8 << 20;
    // Reads split into parts of at least this many chunks.
    private static final int PART_CHUNKS = Integer.getInteger("figwig.partChunks", 8);
    // Chunks one thread claims at a time; the tracks of a group then share a view's decoding.
    private static final int SPAN_CHUNKS = Integer.getInteger("figwig.spanChunks", 4);
    // Over HTTP each read is a round trip, so a thread claims longer stretches.
    private static final int HTTP_SPAN_CHUNKS = Integer.getInteger("figwig.httpSpanChunks", 1 << 20);
    // Over HTTP the parts of a load are requests sent at once; Java keeps 5 idle connections per host.
    private static final int HTTP_PARTS = Integer.getInteger("figwig.httpParts", 4);

    public final String path;
    private final SeekableStream stream;

    public final int tile;
    public final int codec;
    public final boolean trackOrder;
    public final String genome;

    public final String[] chromNames;
    public final long[] chromLengths;
    private final long[] chromTileStart;
    private final Map<String, Integer> chromIds = new HashMap<>();

    public final String[] trackNames;
    public final String[] trackDtypes;
    public final int[] trackGroup;
    public final int[] trackSlot;
    public final List<LinkedHashMap<String, String>> trackAttributes = new ArrayList<>();
    public final double[] statSum, statMin, statMax;
    public final long[] statNonzero;

    private final int nGroups;
    private final int[] groupSize;
    private final int[] groupTrack;   // a track of each group
    private final long[] recordStart, recordEnd;
    private final boolean hasZooms;
    private final FigWigSummaryIndex summaryIndex;   // the regional summary index, or null
    // Reading levels by windows through the regional summary index; false reads them whole.
    private static final boolean SUMMARY_INDEX = !"false".equals(System.getProperty("figwig.summaryIndex"));
    private final AtomicReferenceArray<Group> groups;
    private final Object[] groupLocks;
    private final AtomicIntegerArray groupUsers;   // tracks of each group in use
    private final AtomicIntegerArray trackUsers;   // users of each track
    private volatile int[] zoomFactors;       // level 0 is the summary, factor 1; from the first group with one
    private volatile long[][] chromBinStart;  // per level, per chromosome

    public final long openNanos;
    public final long openBytesRead;
    public final int openReads;


    /**
     * Loaded chunks and blocks by key.  IGV loads the tracks of a view on several threads at once, and
     * the tracks of a group share their chunks: the first thread to want a chunk loads it, the others
     * wait for it.  Bounded by bytes, evicted first in, first out.
     */
    private final ConcurrentHashMap<Long, CompletableFuture<Loaded>> cache =
            new ConcurrentHashMap<>();

    /** One entry of the cache budget every open figwig file shares, evicted first in, first out. */
    private static final class Cached {
        final FigWigFile file;
        final long key;
        final int bytes;

        Cached(FigWigFile file, long key, int bytes) {
            this.file = file;
            this.key = key;
            this.bytes = bytes;
        }
    }

    private static final ConcurrentLinkedQueue<Cached> CACHE_ORDER = new ConcurrentLinkedQueue<>();
    private static final AtomicLong CACHED_BYTES = new AtomicLong();

    private final AtomicInteger users = new AtomicInteger();
    private volatile boolean closed = false;

    private final AtomicLong chunkBytesRead = new AtomicLong();
    private final AtomicLong chunkReads = new AtomicLong();
    private final FileChannel channel;
    private final boolean HTTP;   // http(s), s3 and gs paths: IGV reads them over HTTP range requests
    private final int spanChunks;

    private static class Group {
        int nTracks;
        // position order
        int[] liveTiles;
        int[] kinds;
        int[] sizes;
        int[] raws;
        long[] crcs;
        long[] offsets;
        // track order
        int[] bChrom, bFirst, bLast, bN, bSize, bRaw, bKind;
        long[] bOffset;
        // summary (level 0) and zoom levels: where each stream lies in the record, and once inflated,
        // its bin numbers and sums (tracks per bin).  The summary is read from the file when first
        // asked for; the zoom levels' headers and streams (zoomBytes, from zoomAt to the record's end)
        // when a zoom level or the levels' factors are first asked for.
        int index;
        long recordLength;
        int summaryAt, summaryN, summaryBytes, zoomAt;
        long summaryCrc;
        byte[] zoomBytes;
        volatile int[] factor;
        int[] levelN, levelAt, levelBytes;
        long[] levelCrc;
        AtomicReferenceArray<Level> levelData;
        FigWigSummaryIndex.LevelEntry[] catalog;   // from the regional summary index, when the file has one
    }

    /** One summary or zoom level of a group, inflated: bin numbers, and the sums of the group's tracks per bin. */
    private static final class Level {
        final int[] bins;
        final byte[] raw;     // the inflated stream; sums (float32 for the summary, float64 for zooms) at sumsAt
        final int sumsAt;
        final int width;

        Level(int[] bins, byte[] raw, int sumsAt, int width) {
            this.bins = bins;
            this.raw = raw;
            this.sumsAt = sumsAt;
            this.width = width;
        }

        double sum(int i) {
            return width == 4 ? (float) F32.get(raw, sumsAt + 4 * i) : (double) F64.get(raw, sumsAt + 8 * i);
        }
    }

    public static FigWigFile open(String path) throws IOException {
        return new FigWigFile(path);
    }

    private FigWigFile(String path) throws IOException {
        long t0 = System.nanoTime();
        this.path = path;
        this.stream = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        java.io.File local = new java.io.File(path);
        this.channel = POSITIONAL && local.isFile() ?
                FileChannel.open(local.toPath(), java.nio.file.StandardOpenOption.READ) : null;
        this.HTTP = channel == null && org.igv.util.FileUtils.isRemote(path) && !path.startsWith("ftp:");
        this.spanChunks = HTTP ? HTTP_SPAN_CHUNKS : SPAN_CHUNKS;

        byte[] headerBytes = readAt(0, 128);
        ByteBuffer h = ByteBuffer.wrap(headerBytes).order(ByteOrder.LITTLE_ENDIAN);
        if (h.getInt(0) != MAGIC) {
            throw new IOException("figwig: " + path + " is not a figwig file");
        }
        byte[] zeroed = headerBytes.clone();
        Arrays.fill(zeroed, 104, 108, (byte) 0);
        if (crc32(zeroed, 0, 128) != (h.getInt(104) & 0xffffffffL)) {
            throw new IOException("figwig: the header's CRC-32 does not match");
        }
        long fileSize = stream.length() > 0 ? stream.length() : h.getLong(88);

        int version = h.getShort(4) & 0xffff, flags = h.getShort(6) & 0xffff;
        long tl = h.getInt(8) & 0xffffffffL;
        int codecField = h.getShort(12) & 0xffff;
        long nChroms = h.getInt(16) & 0xffffffffL, nTracks = h.getInt(20) & 0xffffffffL, ng = h.getInt(24) & 0xffffffffL;
        long tablesCrc = h.getInt(28) & 0xffffffffL;
        long chromsOffset = h.getLong(32), tracksOffset = h.getLong(40), indexOffset = h.getLong(48);
        long metadataOffset = h.getLong(64), metadataLength = h.getLong(72), length = h.getLong(88);
        long required = h.getInt(96) & 0xffffffffL, optional = h.getInt(100) & 0xffffffffL;
        check(version == 1, "version " + version + " is not 1");
        check(required == 0, "required feature bits " + required + " are unknown");
        check((flags & ~1) == 0, "unknown flags " + flags);
        check(tl >= 16 && tl <= 65536 && (tl & (tl - 1)) == 0, "tile " + tl);
        check(codecField <= 1, "codec " + codecField);
        check(length >= 128 && length <= fileSize, "length " + length);
        check(128 <= chromsOffset && chromsOffset <= tracksOffset && tracksOffset <= metadataOffset
                && metadataOffset <= indexOffset && indexOffset <= length, "table offsets out of order");
        check(nChroms > 0 && nTracks > 0 && ng > 0, "a count is 0");
        long tablesEnd = indexOffset + 8 * ng;
        check(tablesEnd <= length && tablesEnd - chromsOffset < Integer.MAX_VALUE, "the index directory passes the end");
        check(metadataOffset + metadataLength <= indexOffset, "the metadata passes the index");
        byte[] tables = readAt(chromsOffset, (int) (tablesEnd - chromsOffset));
        check(crc32(tables, 0, tables.length) == tablesCrc, "the tables' CRC-32 does not match");
        ByteBuffer tb = ByteBuffer.wrap(tables).order(ByteOrder.LITTLE_ENDIAN);

        this.tile = (int) tl;
        this.codec = codecField;
        this.trackOrder = (flags & 1) != 0;
        this.hasZooms = (optional & 1) != 0;

        // Chromosomes
        int nc = (int) nChroms;
        chromNames = new String[nc];
        chromLengths = new long[nc];
        chromTileStart = new long[nc + 1];
        int at = 0;
        int chromsEnd = (int) (tracksOffset - chromsOffset);
        for (int c = 0; c < nc; c++) {
            check(at + 2 <= chromsEnd, "the chromosome table passes its end");
            int len = tb.getShort(at) & 0xffff;
            check(at + 2 + len + 8 <= chromsEnd, "the chromosome table passes its end");
            chromNames[c] = new String(tables, at + 2, len, StandardCharsets.UTF_8);
            chromLengths[c] = tb.getLong(at + 2 + len);
            at += 2 + len + 8;
            chromIds.put(chromNames[c], c);
            chromTileStart[c + 1] = chromTileStart[c] + (chromLengths[c] + tile - 1) / tile;
        }

        // Statistics
        int nt = (int) nTracks;
        check((long) nt * 32 <= metadataOffset - tracksOffset, "the statistics table passes its end");
        statSum = new double[nt];
        statMin = new double[nt];
        statMax = new double[nt];
        statNonzero = new long[nt];
        for (int t = 0; t < nt; t++) {
            int o = (int) (tracksOffset - chromsOffset) + 32 * t;
            statSum[t] = tb.getDouble(o);
            statNonzero[t] = tb.getLong(o + 8);
            statMin[t] = tb.getDouble(o + 16);
            statMax[t] = tb.getDouble(o + 24);
        }

        // Metadata
        JSONObject metadata = new JSONObject(new String(tables, (int) (metadataOffset - chromsOffset),
                (int) metadataLength, StandardCharsets.UTF_8));
        genome = metadata.optString("genome", null);
        JSONArray tracks = metadata.getJSONArray("tracks");
        check(tracks.length() == nt, "metadata describes " + tracks.length() + " tracks, the header " + nt);
        trackNames = new String[nt];
        trackDtypes = new String[nt];
        trackGroup = new int[nt];
        trackSlot = new int[nt];
        nGroups = (int) ng;
        groupSize = new int[nGroups];
        for (int t = 0; t < nt; t++) {
            JSONObject record = tracks.getJSONObject(t);
            trackNames[t] = record.getString("name");
            trackDtypes[t] = record.getString("dtype");
            trackGroup[t] = record.getInt("group");
            check(trackGroup[t] >= 0 && trackGroup[t] < nGroups, "track " + t + " has group " + trackGroup[t]);
            trackSlot[t] = groupSize[trackGroup[t]]++;
            LinkedHashMap<String, String> attributes = new LinkedHashMap<>();
            for (String key : record.keySet()) {
                if (!key.equals("name") && !key.equals("dtype") && !key.equals("group")) {
                    attributes.put(key, String.valueOf(record.get(key)));
                }
            }
            trackAttributes.add(attributes);
        }

        groupTrack = new int[nGroups];
        for (int t = nt - 1; t >= 0; t--) {
            groupTrack[trackGroup[t]] = t;
        }
        summaryIndex = SUMMARY_INDEX ? FigWigSummaryIndex.of(new FigWigSummaryIndex.Source() {
            public byte[] readAt(long offset, int n) throws IOException {
                return FigWigFile.this.readAt(offset, n);
            }

            public byte[][] readEach(long[] offsets, int[] lengths) throws IOException {
                if (!HTTP) {
                    byte[][] out = new byte[offsets.length][];
                    for (int i = 0; i < offsets.length; i++) {
                        out[i] = FigWigFile.this.readAt(offsets[i], lengths[i]);
                    }
                    return out;
                }
                List<Callable<byte[]>> tasks = new ArrayList<>();
                for (int i = 0; i < offsets.length; i++) {
                    final long o = offsets[i];
                    final int n = lengths[i];
                    tasks.add(() -> FigWigFile.this.readAt(o, n));
                }
                return parallel(tasks, PREFETCHER).toArray(new byte[0][]);
            }
        }, headerBytes, optional, metadata, chromsOffset, nGroups) : null;

        // The index directory: records lie back to back from the directory's end to the file's length.
        recordStart = new long[nGroups];
        recordEnd = new long[nGroups];
        int dir = (int) (indexOffset - chromsOffset);
        for (int g = 0; g < nGroups; g++) {
            recordStart[g] = tb.getLong(dir + 8 * g);
        }
        for (int g = 0; g < nGroups; g++) {
            recordEnd[g] = g + 1 < nGroups ? recordStart[g + 1] : length;
            check(recordStart[g] >= tablesEnd && recordStart[g] <= recordEnd[g] && recordEnd[g] <= length
                    && recordEnd[g] - recordStart[g] < Integer.MAX_VALUE, "group " + g + "'s index record is out of place");
        }
        groups = new AtomicReferenceArray<>(nGroups);
        groupUsers = new AtomicIntegerArray(nGroups);
        trackUsers = new AtomicIntegerArray(nt);
        groupLocks = new Object[nGroups];
        for (int g = 0; g < nGroups; g++) {
            groupLocks[g] = new Object();
        }

        openBytesRead = 128 + tables.length;
        openReads = 2;
        openNanos = System.nanoTime() - t0;
    }

    /** A group's index record, read and parsed the first time one of its tracks is read. */
    private Group group(int g) throws IOException {
        Group x = groups.get(g);
        if (x != null) {
            return x;
        }
        synchronized (groupLocks[g]) {
            x = groups.get(g);
            if (x == null) {
                x = parseGroup(g);
                groups.set(g, x);
            }
            return x;
        }
    }

    // The first read of a group's record: its runs and live tiles usually fit, else the rest follows.
    private static final int RECORD_HEAD = Integer.getInteger("figwig.recordHead", 1 << 16);

    private Group parseGroup(int g) throws IOException {
        long len = recordEnd[g] - recordStart[g];
        if (trackOrder) {
            return parseBlockGroup(g, readAt(recordStart[g], (int) len));
        }
        Group group = new Group();
        group.nTracks = groupSize[g];
        group.index = g;
        group.recordLength = len;
        check(len >= 8, "group " + g + "'s index record is cut short");
        // Live tiles take about 10 bytes each compressed: one read usually brings all of them.
        long guess = RECORD_HEAD + 10 * chromTileStart[chromNames.length];
        byte[] rb = readAt(recordStart[g], (int) Math.min(len, Math.min(guess, 64L << 20)));
        ByteBuffer r = ByteBuffer.wrap(rb).order(ByteOrder.LITTLE_ENDIAN);
        long n = r.getInt(0) & 0xffffffffL;
        long crc = r.getInt(4) & 0xffffffffL;
        check(8 + 32 * n + 8 <= len, "group " + g + "'s runs pass its record");
        int q = (int) (8 + 32 * n);
        if (q + 8 > rb.length) {
            rb = extend(g, rb, q + 8);
            r = ByteBuffer.wrap(rb).order(ByteOrder.LITTLE_ENDIAN);
        }
        long nLive = r.getInt(q) & 0xffffffffL;
        long nBytes = r.getInt(q + 4) & 0xffffffffL;
        check(q + 8 + nBytes <= len, "group " + g + "'s live tiles pass its record");
        int p = (int) (q + 8 + nBytes);
        // The live tiles and the summary's header
        int need = (int) Math.min(len, (long) p + 12);
        if (need > rb.length) {
            rb = extend(g, rb, need);
            r = ByteBuffer.wrap(rb).order(ByteOrder.LITTLE_ENDIAN);
        }
        check(crc32(rb, 8, (int) (32 * n + 8 + nBytes)) == crc, "group " + g + "'s index record does not match its CRC-32");
        check(nLive <= chromTileStart[chromNames.length], "group " + g + " lists more live tiles than the genome has");
        check(n <= chromTileStart[chromNames.length], "group " + g + " lists more runs than the genome has tiles");
        byte[] live = inflate(Arrays.copyOfRange(rb, q + 8, p), (int) (17 * nLive));
        readLiveTiles(r, (int) n, live, (int) nLive, group);
        check(p + 12 <= len, "group " + g + "'s summary passes its record");
        group.summaryN = r.getInt(p);
        group.summaryBytes = r.getInt(p + 4);
        group.summaryCrc = r.getInt(p + 8) & 0xffffffffL;
        group.summaryAt = p + 12;
        long zoomAt = (long) p + 12 + group.summaryBytes;
        check(zoomAt <= len && zoomAt >= p + 12, "group " + g + "'s summary passes its record");
        group.zoomAt = (int) zoomAt;
        return group;
    }

    /** The record's bytes up to {@code need}, reading what {@code have} lacks. */
    private byte[] extend(int g, byte[] have, int need) throws IOException {
        byte[] more = readAt(recordStart[g] + have.length, need - have.length);
        byte[] all = Arrays.copyOf(have, need);
        System.arraycopy(more, 0, all, have.length, more.length);
        return all;
    }

    /** A track-order group's record, read whole. */
    private Group parseBlockGroup(int g, byte[] rb) throws IOException {
        ByteBuffer r = ByteBuffer.wrap(rb).order(ByteOrder.LITTLE_ENDIAN);
        Group group = new Group();
        group.nTracks = groupSize[g];
        group.index = g;
        group.recordLength = rb.length;
        check(rb.length >= 8, "group " + g + "'s index record is cut short");
        long n = r.getInt(0) & 0xffffffffL;
        long crc = r.getInt(4) & 0xffffffffL;
        check(8 + 40 * n <= rb.length, "group " + g + "'s block records pass its record");
        check(crc32(rb, 8, (int) (40 * n)) == crc, "group " + g + "'s block records do not match their CRC-32");
        readBlockRecords(r, (int) n, group);
        int p = (int) (8 + 40 * n);
        check(p + 12 <= rb.length, "group " + g + "'s summary passes its record");
        group.summaryN = r.getInt(p);
        group.summaryBytes = r.getInt(p + 4);
        group.summaryCrc = r.getInt(p + 8) & 0xffffffffL;
        group.summaryAt = p + 12;
        long zoomAt = (long) p + 12 + group.summaryBytes;
        check(zoomAt <= rb.length && zoomAt >= p + 12, "group " + g + "'s summary passes its record");
        group.zoomAt = (int) zoomAt;
        return group;
    }

    /** A group's zoom levels' headers, read the first time a level or the factors are asked for. */
    private int[] zooms(Group group) throws IOException {
        int[] f = group.factor;
        if (f != null) {
            return f;
        }
        synchronized (group) {
            if (group.factor != null) {
                return group.factor;
            }
            int g = group.index;
            if (summaryIndex != null && summaryIndex.usable()) {
                // The catalog places every level in the record: the zoom section is not read whole.
                long[] bins = new long[13];
                for (int k = 0; k < bins.length; k++) {
                    for (int c = 0; c < chromNames.length; c++) {
                        long t = chromTileStart[c + 1] - chromTileStart[c];
                        bins[k] += (t + (1L << (2 * k)) - 1) >> (2 * k);
                    }
                }
                FigWigSummaryIndex.LevelEntry[] cat = summaryIndex.catalog(g, group.nTracks, recordStart[g],
                        group.recordLength, group.summaryAt, group.summaryN, group.summaryBytes, group.summaryCrc,
                        group.zoomAt, bins);
                int levels = cat.length;
                int[] factor = new int[levels];
                group.levelN = new int[levels];
                group.levelAt = new int[levels];
                group.levelBytes = new int[levels];
                group.levelCrc = new long[levels];
                group.levelData = new AtomicReferenceArray<>(levels);
                for (int k = 0; k < levels; k++) {
                    factor[k] = 1 << (2 * k);
                    group.levelN[k] = cat[k].sourceRows;
                    group.levelBytes[k] = cat[k].sourceBytes;
                    group.levelCrc[k] = cat[k].sourceCrc;
                }
                group.levelAt[0] = group.summaryAt;
                group.catalog = cat;
                group.factor = factor;
                return factor;
            }
            int levels = 1;
            byte[] zb = new byte[0];
            if (hasZooms) {
                check(group.zoomAt + 4 <= group.recordLength, "group " + g + "'s zoom section passes its record");
                zb = readAt(recordStart[g] + group.zoomAt, (int) (group.recordLength - group.zoomAt));
                levels += (int) I32.get(zb, 0);
                check(levels >= 1 && levels <= 13, "group " + g + " has " + (levels - 1) + " zoom levels");
            }
            int[] factor = new int[levels];
            group.levelN = new int[levels];
            group.levelAt = new int[levels];
            group.levelBytes = new int[levels];
            group.levelCrc = new long[levels];
            group.levelData = new AtomicReferenceArray<>(levels);
            factor[0] = 1;
            group.levelN[0] = group.summaryN;
            group.levelBytes[0] = group.summaryBytes;
            group.levelCrc[0] = group.summaryCrc;
            group.levelAt[0] = group.summaryAt;   // in the record; levels 1 and up: in zoomBytes
            int z = 4;
            for (int k = 1; k < levels; k++) {
                check(z + 16 <= zb.length, "group " + g + "'s zoom level " + k + " passes its record");
                factor[k] = (int) I32.get(zb, z);
                group.levelN[k] = (int) I32.get(zb, z + 4);
                group.levelBytes[k] = (int) I32.get(zb, z + 8);
                group.levelCrc[k] = (int) I32.get(zb, z + 12) & 0xffffffffL;
                group.levelAt[k] = z + 16;
                check(group.levelBytes[k] >= 0 && (long) z + 16 + group.levelBytes[k] <= zb.length,
                        "group " + g + "'s zoom level " + k + " passes its record");
                z += 16 + group.levelBytes[k];
            }
            group.zoomBytes = zb;
            group.factor = factor;
            return factor;
        }
    }

    /**
     * The zoom factors and the bins per chromosome.  The factors are those of the first group that keeps
     * a summary, which every group keeping one must match; group 0's when none does.
     */
    private void zoomGrid() throws IOException {
        if (zoomFactors != null) {
            return;
        }
        synchronized (this) {
            if (zoomFactors == null) {
                int first = 0;
                for (int g = 0; g < nGroups; g++) {
                    if (group(g).summaryBytes > 0) {
                        first = g;
                        break;
                    }
                }
                int[] f = zooms(group(first));
                long[][] start = new long[f.length][chromNames.length + 1];
                for (int k = 0; k < f.length; k++) {
                    for (int c = 0; c < chromNames.length; c++) {
                        long tiles = chromTileStart[c + 1] - chromTileStart[c];
                        start[k][c + 1] = start[k][c] + (tiles + f[k] - 1) / f[k];
                    }
                }
                chromBinStart = start;
                zoomFactors = f;
            }
        }
    }

    /** A group's summary (level 0) or zoom level, inflated the first time it is asked for. */
    private Level level(Group group, int k) throws IOException {
        Level l = group.levelData.get(k);
        if (l != null) {
            return l;
        }
        synchronized (group) {
            l = group.levelData.get(k);
            if (l != null) {
                return l;
            }
            int n = group.levelN[k];
            int tracks = group.nTracks;
            int w = k == 0 ? 4 : 8;
            long levelBins = 0;
            for (int c = 0; c < chromNames.length; c++) {
                long t = chromTileStart[c + 1] - chromTileStart[c];
                levelBins += (t + group.factor[k] - 1) / Math.max(1, group.factor[k]);
            }
            check(n >= 0 && n <= levelBins && group.factor[k] >= 1, "a summary or zoom level lists more bins than it has");
            check(group.levelBytes[k] >= 0, "a summary or zoom level's size is negative");
            int[] bins = new int[n];
            byte[] b = new byte[0];
            if (group.levelBytes[k] > 0) {
                byte[] stream;
                if (k == 0) {
                    stream = readAt(recordStart[group.index] + group.levelAt[0], group.levelBytes[0]);
                } else if (group.catalog != null) {
                    stream = summaryIndex.stream(group.index, group.catalog, k);
                } else {
                    stream = Arrays.copyOfRange(group.zoomBytes, group.levelAt[k], group.levelAt[k] + group.levelBytes[k]);
                }
                check(crc32(stream, 0, stream.length) == group.levelCrc[k],
                        "a summary or zoom level does not match its CRC-32");
                long raw = (long) n * (4 + (long) w * tracks);
                check(raw < Integer.MAX_VALUE, "a summary or zoom level is too large");
                b = inflate(stream, (int) raw);
                for (int i = 0; i < n; i++) {
                    bins[i] = (int) I32.get(b, 4 * i);
                }
            } else {
                check(n == 0, "a summary or zoom level lists bins but holds no stream");
            }
            l = new Level(bins, b, 4 * n, w);
            group.levelData.set(k, l);
            return l;
        }
    }

    private void readLiveTiles(ByteBuffer r, int nRuns, byte[] live, int n, Group group) {
        group.liveTiles = new int[n];
        group.kinds = new int[n];
        group.sizes = new int[n];
        group.raws = new int[n];
        group.crcs = new long[n];
        group.offsets = new long[n];
        // The largest a chunk of this group can be, decoded and in the file (docs/format.rst).
        long most = 8L * group.nTracks * tile + 16L * group.nTracks + 64;
        long tiles = chromTileStart[chromNames.length];
        for (int i = 0; i < n; i++) {
            group.liveTiles[i] = (int) I32.get(live, 4 * i);
            group.kinds[i] = live[4 * n + i] & 0xff;
            group.sizes[i] = (int) I32.get(live, 5 * n + 4 * i);
            group.raws[i] = (int) I32.get(live, 9 * n + 4 * i);
            group.crcs[i] = (int) I32.get(live, 13 * n + 4 * i) & 0xffffffffL;
            long t = group.liveTiles[i] & 0xffffffffL;
            check(t < tiles && (i == 0 || t > (group.liveTiles[i - 1] & 0xffffffffL)), "live tiles out of order or range");
            check(group.raws[i] > 0 && group.raws[i] <= most && group.sizes[i] > 0 && group.sizes[i] <= most + 64,
                    "a chunk's size is out of range");
        }
        // Each live tile's offset: walk the runs in tile order.
        long[][] runs = new long[nRuns][];
        for (int j = 0; j < nRuns; j++) {
            int o = 8 + 32 * j;
            runs[j] = new long[]{r.getLong(o), r.getLong(o + 8), r.getLong(o + 16), r.getLong(o + 24)};
        }
        Arrays.sort(runs, Comparator.comparingLong(x -> x[0]));
        int k = 0;
        for (long[] run : runs) {
            if (k < n && group.liveTiles[k] < run[0]) {
                throw new IllegalStateException("figwig: live tile " + group.liveTiles[k] + " lies in no run");
            }
            long offset = run[2];
            while (k < n && group.liveTiles[k] <= run[1]) {
                group.offsets[k] = offset;
                offset += group.sizes[k];
                k++;
            }
            check(offset == run[2] + run[3], "a run's bytes do not match its chunks");
        }
        check(k == n, "live tiles past the last run");
    }

    private void readBlockRecords(ByteBuffer r, int n, Group group) {
        group.bChrom = new int[n];
        group.bFirst = new int[n];
        group.bLast = new int[n];
        group.bN = new int[n];
        group.bSize = new int[n];
        group.bRaw = new int[n];
        group.bKind = new int[n];
        group.bOffset = new long[n];
        for (int i = 0; i < n; i++) {
            int o = 8 + 40 * i;
            group.bChrom[i] = r.getInt(o);
            group.bFirst[i] = r.getInt(o + 4);
            group.bLast[i] = r.getInt(o + 8);
            group.bN[i] = r.getInt(o + 12);
            group.bOffset[i] = r.getLong(o + 16);
            group.bSize[i] = r.getInt(o + 24);
            group.bRaw[i] = r.getInt(o + 28);
            group.bKind[i] = r.getInt(o + 32);
            check((group.bKind[i] & 7) == 1, "a block whose kind is not 1");
            check(group.bN[i] >= 1 && group.bN[i] <= 65536 && group.bRaw[i] > 0 && group.bRaw[i] <= 12L * group.bN[i] + 64
                    && group.bSize[i] > 0 && group.bSize[i] <= group.bRaw[i] + 1024 && group.bChrom[i] >= 0
                    && group.bChrom[i] < chromNames.length, "a block record is out of range");
        }
    }

    // ------------------------------------------------------------------------------------------
    // Values

    /**
     * Runs of one value, nonzero, in ascending position; NaN is kept.
     */
    public static class Runs {
        public int n = 0;
        public int[] starts;
        public int[] ends;
        public double[] values;

        Runs() {
            this(256);
        }

        Runs(int capacity) {
            capacity = Math.max(16, capacity);
            starts = new int[capacity];
            ends = new int[capacity];
            values = new double[capacity];
        }

        void reserve(int more) {
            if (n + more > starts.length) {
                int c = Math.max(n + more, 2 * starts.length);
                starts = Arrays.copyOf(starts, c);
                ends = Arrays.copyOf(ends, c);
                values = Arrays.copyOf(values, c);
            }
        }

        /** Append the runs of a later part of the same read, merging across the seam as add would. */
        void append(Runs other) {
            if (other.n == 0) {
                return;
            }
            add(other.starts[0], other.ends[0], other.values[0]);
            int m = other.n - 1;
            reserve(m);
            System.arraycopy(other.starts, 1, starts, n, m);
            System.arraycopy(other.ends, 1, ends, n, m);
            System.arraycopy(other.values, 1, values, n, m);
            n += m;
        }

        void add(int start, int end, double value) {
            if (value == 0 || end <= start) {
                return;
            }
            if (n > 0 && ends[n - 1] == start && Double.compare(values[n - 1], value) == 0) {
                ends[n - 1] = end;
                return;
            }
            if (n == starts.length) {
                starts = Arrays.copyOf(starts, 2 * n);
                ends = Arrays.copyOf(ends, 2 * n);
                values = Arrays.copyOf(values, 2 * n);
            }
            starts[n] = start;
            ends[n] = end;
            values[n] = value;
            n++;
        }
    }

    public Integer chromId(String name) {
        return chromIds.get(name);
    }

    public int nTracks() {
        return trackNames.length;
    }

    /**
     * The nonzero values of a track over [start, end) of a chromosome, as runs.
     */
    public Runs read(int track, int chrom, int start, int end) throws IOException {
        Runs runs = new Runs();
        end = (int) Math.min(end, chromLengths[chrom]);
        start = Math.max(0, start);
        if (end <= start) {
            return runs;
        }
        try {
            if (trackOrder) {
                readBlocks(track, chrom, start, end, runs);
            } else {
                Runs batched = batched(track, chrom, start, end);
                if (batched != null) {
                    return batched;
                }
                readTiles(track, chrom, start, end, runs);
            }
        } catch (IndexOutOfBoundsException | NegativeArraySizeException | IllegalArgumentException e) {
            throw new IOException("figwig: " + path + " is corrupt: " + e, e);
        }
        return runs;
    }

    /**
     * Load, without decoding, the chunks or blocks a read of [start, end) of a track will need, so
     * that reads of parts of the range that follow find them cached.  Over HTTP, also start loading
     * the same range for the file's other groups that have been read before, on background threads,
     * so that the groups' round trips overlap instead of following one another.
     */
    public void prefetch(int track, int chrom, int start, int end) throws IOException {
        end = (int) Math.min(end, chromLengths[chrom]);
        start = Math.max(0, start);
        if (end <= start) {
            return;
        }
        int g = trackGroup[track];
        if (HTTP) {
            activeGroups.add(g);
            for (int x = 0; x < nGroups; x++) {
                if (groupUsers.get(x) > 0) {
                    activeGroups.add(x);
                }
            }
            String range = chrom + ":" + start + ":" + end;
            for (int other : activeGroups) {
                // One background load per group and range: every track of every group asks.
                if (other != g && !range.equals(prefetched.put(other, range))) {
                    final int og = other, s0 = start, e0 = end;
                    PREFETCHER.submit(() -> {
                        try {
                            loadRange(og, groupTrack[og], chrom, s0, e0);
                        } catch (Exception ignored) {
                            // A reader of that group loads the range again and reports the error.
                        }
                    });
                }
            }
        }
        try {
            loadRange(g, track, chrom, start, end);
        } catch (IndexOutOfBoundsException | NegativeArraySizeException | IllegalArgumentException e) {
            throw new IOException("figwig: " + path + " is corrupt: " + e, e);
        }
    }

    private void loadRange(int g, int track, int chrom, int start, int end) throws IOException {
        Group group = group(g);
        if (trackOrder) {
            int[] r = blockRange(group, chrom, start, end);
            loadBlocks(group, track, r[0], r[1]);
        } else {
            int[] r = tileRange(group, chrom, start, end);
            int k = r[0], k1 = r[1];
            int parts = parts(group, k, k1);
            if (parts <= 1) {
                loadTiles(group, g, k, k1);
                return;
            }
            List<Callable<Object>> tasks = new ArrayList<>(parts);
            for (int p = 0; p < parts; p++) {
                final int a = k + (int) ((long) (k1 - k) * p / parts), b = k + (int) ((long) (k1 - k) * (p + 1) / parts);
                tasks.add(() -> loadTiles(group, g, a, b));
            }
            parallel(tasks, HTTP ? PREFETCHER : DECODER);
        }
    }

    // Over HTTP, a part is one request: at least HTTP_PART_BYTES each, at most HTTP_PARTS.
    private static final long HTTP_PART_BYTES = Long.getLong("figwig.httpPartBytes", 256 << 10);

    /** Into how many parts to split the load of live tiles [k, k1) of a position-order group. */
    private int parts(Group group, int k, int k1) {
        if (k1 - k < 2) {
            return 1;
        }
        if (HTTP) {
            long bytes = group.offsets[k1 - 1] + group.sizes[k1 - 1] - group.offsets[k];
            return (int) Math.max(1, Math.min(Math.min(HTTP_PARTS, k1 - k), bytes / HTTP_PART_BYTES));
        }
        if (DECODER == null) {
            return 1;
        }
        // The decoding threads are shared among the reads running at once.
        int active = Math.max(1, ACTIVE_READS.get());
        int share = Math.max(1, (SHARE * DECODE_THREADS + active - 1) / active);
        return Math.min(share, (k1 - k) / PART_CHUNKS);
    }

    // How many times the decoding threads a lone read splits into: parts are claimed cooperatively, so
    // extra parts keep every thread busy when parts take unequal time.
    private static final int SHARE = Integer.getInteger("figwig.share", 2);

    /** Reads running now, in every open file. */
    static final AtomicInteger ACTIVE_READS = new AtomicInteger();

    /** Groups read so far or with tracks in use, which an HTTP prefetch loads too, and the range last prefetched for each. */
    private final java.util.Set<Integer> activeGroups = ConcurrentHashMap.newKeySet();
    private final ConcurrentHashMap<Integer, String> prefetched = new ConcurrentHashMap<>();


    /**
     * Every base of [start, end), 0 where the track has no value.
     */
    public double[] readDense(int track, int chrom, int start, int end) throws IOException {
        double[] out = new double[end - start];
        Runs runs = read(track, chrom, start, end);
        for (int i = 0; i < runs.n; i++) {
            Arrays.fill(out, runs.starts[i] - start, runs.ends[i] - start, runs.values[i]);
        }
        return out;
    }

    /** The live tiles of a group over [start, end) of a chromosome: [k, k1) of its live list. */
    private int[] tileRange(Group group, int chrom, int start, int end) {
        long first = chromTileStart[chrom] + start / tile;
        long last = chromTileStart[chrom] + (end - 1) / tile;
        int k = lowerBound(group.liveTiles, first);
        int k1 = k;
        while (k1 < group.liveTiles.length && group.liveTiles[k1] <= last) {
            k1++;
        }
        return new int[]{k, k1};
    }

    private Loaded[] loadTiles(Group group, int g, int k, int k1) throws IOException {
        return load(k, k1, i -> ((long) g << 32) | i, group.offsets, group.sizes,
                (i, bytes) -> loadChunk(group, i, bytes));
    }

    private void readTiles(int track, int chrom, int start, int end, Runs runs) throws IOException {
        Group group = group(trackGroup[track]);
        int[] range = tileRange(group, chrom, start, end);
        int k = range[0], k1 = range[1];
        if (HTTP || k1 - k < 2 * PART_CHUNKS) {
            readTilesSplit(group, track, chrom, start, end, k, k1, runs);
            return;
        }
        // Reads large enough to split are counted, so that reads running at once share the threads.
        ACTIVE_READS.incrementAndGet();
        try {
            readTilesSplit(group, track, chrom, start, end, k, k1, runs);
        } finally {
            ACTIVE_READS.decrementAndGet();
        }
    }

    // ------------------------------------------------------------------------------------------
    // Reads of many tracks of a group over one range (a view of many cells): the first read decodes
    // every track of the group in use in one pass over the chunks, and the others take their runs.
    // Read one by one, each track walks every chunk again, and the walk, not the data, costs the time.

    private static final int BATCH_TRACKS = Integer.getInteger("figwig.batchTracks", 16);
    // At most this many runs in one batch (a run takes 16 bytes); a larger one is read track by track.
    private static final long BATCH_RUNS = Long.getLong("figwig.batchRuns", 1L << 23);
    // Runs in the batches kept, as made; past this the oldest batches are dropped.
    private static final long BATCH_HELD = Long.getLong("figwig.batchHeld", 1L << 23);

    /**
     * Recent batches by group, chromosome and range, oldest first: each track's runs, removed when the
     * track takes them.  A view of many tiles has a batch per tile in use at once.  A batch every track
     * has taken stays, so that a track reading the range again reads it alone instead of making the
     * whole batch again.  The runs each batch was made with, and their total, under the same lock.
     */
    private final LinkedHashMap<String, CompletableFuture<ConcurrentHashMap<Integer, Runs>>> batches =
            new LinkedHashMap<>();
    private final Map<String, Long> batchRuns = new HashMap<>();
    private long heldRuns = 0;

    /** Summary windows read through the regional summary index and kept now; for tests. */
    int summaryWindowsKept() {
        return summaryIndex == null ? 0 : summaryIndex.windowsKept();
    }

    /** Batches kept now; for tests. */
    int batchesKept() {
        synchronized (batches) {
            return batchRuns.size();
        }
    }

    /** The track's runs from a batch of its group's tracks in use, or null to read it alone. */
    private Runs batched(int track, int chrom, int start, int end) throws IOException {
        int g = trackGroup[track];
        if (BATCH_TRACKS <= 0 || groupUsers.get(g) < BATCH_TRACKS || trackUsers.get(track) <= 0) {
            return null;
        }
        Group group = group(g);
        int[] range = tileRange(group, chrom, start, end);
        if (range[1] - range[0] < 2 * PART_CHUNKS) {
            return null;   // a short read alone costs less than a batch
        }
        String key = g + ":" + chrom + ":" + start + ":" + end;
        CompletableFuture<ConcurrentHashMap<Integer, Runs>> mine =
                new CompletableFuture<>(), held;
        synchronized (batches) {
            held = batches.get(key);
            if (held == null) {
                batches.put(key, mine);
            }
        }
        if (held == null) {
            ConcurrentHashMap<Integer, Runs> made;
            try {
                made = readBatch(g, chrom, start, end);
            } catch (IOException | RuntimeException e) {
                // Every track then reads alone, and a damaged chunk fails only the reads it fails alone.
                made = null;
            }
            long n = 0;
            if (made != null) {
                for (Runs r : made.values()) {
                    n += r.n + 1;
                }
            }
            synchronized (batches) {
                if (batches.get(key) == mine) {
                    batchRuns.put(key, n);
                    heldRuns += n;
                    Iterator<Map.Entry<String, CompletableFuture<ConcurrentHashMap<Integer, Runs>>>> it =
                            batches.entrySet().iterator();
                    while ((heldRuns > BATCH_HELD || batches.size() > 256) && it.hasNext()) {
                        Map.Entry<String, CompletableFuture<ConcurrentHashMap<Integer, Runs>>> e = it.next();
                        if (e.getValue() == mine || !batchRuns.containsKey(e.getKey())) {
                            continue;   // the one just made, or one still being made
                        }
                        heldRuns -= batchRuns.remove(e.getKey());
                        it.remove();
                    }
                }
            }
            mine.complete(made);
            held = mine;
        }
        ConcurrentHashMap<Integer, Runs> got = held.join();
        return got == null ? null : got.remove(track);
    }

    /** Every track of group g in use over [start, end), in one pass over the chunks; null if too many runs. */
    private ConcurrentHashMap<Integer, Runs> readBatch(int g, int chrom, int start, int end) throws IOException {
        Group group = group(g);
        int[] range = tileRange(group, chrom, start, end);
        int k = range[0], k1 = range[1];
        int[] tracks = new int[group.nTracks];
        int n = 0;
        for (int t = 0; t < trackGroup.length; t++) {
            if (trackGroup[t] == g && trackUsers.get(t) > 0) {
                tracks[n++] = t;
            }
        }
        tracks = Arrays.copyOf(tracks, n);
        int parts = Math.max(1, parts(group, k, k1));
        List<Callable<Runs[]>> tasks = new ArrayList<>(parts);
        final int[] members = tracks;
        for (int p = 0; p < parts; p++) {
            final int a = k + (int) ((long) (k1 - k) * p / parts), b = k + (int) ((long) (k1 - k) * (p + 1) / parts);
            tasks.add(() -> decodeBatch(group, members, chrom, start, end, a, b));
        }
        List<Runs[]> done = parallel(tasks, HTTP ? PREFETCHER : DECODER);
        long total = 0;
        for (Runs[] r : done) {
            if (r == null) {
                return null;
            }
            for (Runs x : r) {
                total += x.n;
            }
        }
        if (total > BATCH_RUNS) {
            return null;
        }
        ConcurrentHashMap<Integer, Runs> out = new ConcurrentHashMap<>(2 * n);
        for (int j = 0; j < n; j++) {
            Runs runs = done.get(0)[j];
            for (int p = 1; p < done.size(); p++) {
                runs.append(done.get(p)[j]);
            }
            out.put(members[j], runs);
        }
        return out;
    }

    /**
     * The runs of each of the tracks over chunks [a, b), chunk by chunk; null when the chunks' entries
     * (dense: bases) pass BATCH_RUNS.  readBatch checks the runs made against the bound too.
     */
    private Runs[] decodeBatch(Group group, int[] tracks, int chrom, int start, int end, int a, int b) throws IOException {
        Loaded[] loaded = loadTiles(group, trackGroup[tracks[0]], a, b);
        long most = 0;
        for (int i = a; i < b; i++) {
            Loaded l = loaded[i - a];
            int kind = group.kinds[i] & 7;
            if (kind == 2 || kind == 4) {
                most += group.nTracks == 1 ? l.slotEntry[1] : l.slotEntry[group.nTracks];
            } else if (kind == 3) {
                int nPresent = (int) I32.get(l.raw, 0);
                most += (int) I32.get(l.raw, 4 + 8 * nPresent);
            } else if (kind != 0) {
                most += (long) Math.min(tile, end - start) * tracks.length;
            }
            if (most > BATCH_RUNS) {
                return null;
            }
        }
        Runs[] out = new Runs[tracks.length];
        int[] slots = new int[tracks.length];
        for (int j = 0; j < tracks.length; j++) {
            out[j] = new Runs(16);
            slots[j] = trackSlot[tracks[j]];
        }
        for (int i = a; i < b; i++) {
            int tileStart = (int) ((group.liveTiles[i] - chromTileStart[chrom]) * tile);
            int tileEnd = (int) Math.min((long) tileStart + tile, chromLengths[chrom]);
            int from = Math.max(start, tileStart), to = Math.min(end, tileEnd);
            Loaded l = loaded[i - a];
            int kindByte = group.kinds[i];
            for (int j = 0; j < tracks.length; j++) {
                decodeChunk(l, kindByte, group.nTracks, slots[j], tile, tileStart, tileEnd, from, to, out[j]);
            }
        }
        return out;
    }

    private void readTilesSplit(Group group, int track, int chrom, int start, int end, int k, int k1, Runs runs) throws IOException {
        int parts = parts(group, k, k1);
        if (parts <= 1) {
            readTiles(group, track, chrom, start, end, k, k1, runs);
            return;
        }
        // Split the chunks into parts loaded and decoded on the decoding pool; the runs of the parts,
        // joined in order, are those one pass would make.
        List<Callable<Runs>> tasks = new ArrayList<>(parts);
        for (int p = 0; p < parts; p++) {
            final int a = k + (int) ((long) (k1 - k) * p / parts), b = k + (int) ((long) (k1 - k) * (p + 1) / parts);
            tasks.add(() -> {
                Runs r = new Runs();
                readTiles(group, track, chrom, start, end, a, b, r);
                return r;
            });
        }
        for (Runs r : parallel(tasks, HTTP ? PREFETCHER : DECODER)) {
            runs.append(r);
        }
    }

    private void readTiles(Group group, int track, int chrom, int start, int end, int k, int k1, Runs runs) throws IOException {
        int slot = trackSlot[track];
        final int g = trackGroup[track];
        Loaded[] loaded = loadTiles(group, g, k, k1);
        // The most runs this read can make: a slot's entries in sparse, runs and matrix chunks, its
        // bases in dense ones.
        long most = 0;
        for (int i = k; i < k1; i++) {
            Loaded l = loaded[i - k];
            int kind = group.kinds[i] & 7;
            if (kind == 2 || kind == 4) {
                most += (group.nTracks == 1 ? l.slotEntry[1] : l.slotEntry[slot + 1] - l.slotEntry[slot]) + 1;
            } else {
                most += Math.min(tile, end - start);
            }
        }
        runs.reserve((int) Math.min(most, 1 << 24));
        for (int i = k; i < k1; i++) {
            int tileStart = (int) ((group.liveTiles[i] - chromTileStart[chrom]) * tile);
            int tileEnd = (int) Math.min((long) tileStart + tile, chromLengths[chrom]);
            decodeChunk(loaded[i - k], group.kinds[i], group.nTracks, slot, tile, tileStart, tileEnd,
                    Math.max(start, tileStart), Math.min(end, tileEnd), runs);
        }
    }





    @FunctionalInterface
    private interface Single {
        /** Load chunk or block {@code index}; {@code bytes} are its bytes in the file if already read, else null. */
        Loaded load(int index, byte[] bytes) throws IOException;
    }

    /**
     * The chunks or blocks [a, b), each loaded once however many threads want it.  Walking the range,
     * this thread claims a short stretch of chunks nobody holds that lie back to back in the file,
     * reads it with one read and decodes it, then goes on; a chunk another thread claimed is skipped
     * and waited for at the end.  Threads reading tracks of one group so share the decoding of its
     * chunks instead of all waiting on the first one to arrive.
     */
    private Loaded[] load(int a, int b, IntToLongFunction key, long[] offsets, int[] sizes,
                          Single single) throws IOException {
        int n = b - a;
        @SuppressWarnings("unchecked")
        CompletableFuture<Loaded>[] futures = new CompletableFuture[n];
        int i = 0;
        while (i < n) {
            long k0 = key.applyAsLong(a + i);
            CompletableFuture<Loaded> held = cache.get(k0);
            if (held != null) {
                futures[i++] = held;
                continue;
            }
            CompletableFuture<Loaded> f = new CompletableFuture<>();
            held = cache.putIfAbsent(k0, f);
            if (held != null) {
                futures[i++] = held;
                continue;
            }
            futures[i] = f;
            int j = i + 1;
            long end = offsets[a + i] + sizes[a + i];
            while (j < n && j - i < spanChunks && offsets[a + j] == end && end - offsets[a + i] < MAX_SPAN) {
                CompletableFuture<Loaded> g = claim(key.applyAsLong(a + j));
                if (g == null) {
                    break;
                }
                futures[j] = g;
                end += sizes[a + j];
                j++;
            }
            try {
                if (j - i == 1) {
                    complete(key.applyAsLong(a + i), futures[i], single.load(a + i, null));
                } else {
                    byte[] span = readAt(offsets[a + i], (int) (end - offsets[a + i]));
                    for (int x = i; x < j; x++) {
                        int at = (int) (offsets[a + x] - offsets[a + i]);
                        complete(key.applyAsLong(a + x), futures[x],
                                single.load(a + x, Arrays.copyOfRange(span, at, at + sizes[a + x])));
                    }
                }
            } catch (Throwable t) {
                for (int x = i; x < j; x++) {
                    if (!futures[x].isDone()) {
                        futures[x].completeExceptionally(t);
                        cache.remove(key.applyAsLong(a + x), futures[x]);
                    }
                }
                if (t instanceof IOException) throw (IOException) t;
                if (t instanceof RuntimeException) throw (RuntimeException) t;
                throw new IOException(t);
            }
            i = j;
        }
        Loaded[] out = new Loaded[n];
        for (int x = 0; x < n; x++) {
            try {
                out[x] = futures[x].join();
            } catch (CompletionException e) {
                Throwable c = e.getCause();
                if (c instanceof IOException) throw (IOException) c;
                if (c instanceof RuntimeException) throw (RuntimeException) c;
                throw new IOException(c);
            }
        }
        return out;
    }

    /** A new future for the key if nobody holds one, which this thread must then complete; else null. */
    private CompletableFuture<Loaded> claim(long key) {
        if (cache.containsKey(key)) {
            return null;
        }
        CompletableFuture<Loaded> f = new CompletableFuture<>();
        return cache.putIfAbsent(key, f) == null ? f : null;
    }

    private void complete(long key, CompletableFuture<Loaded> f, Loaded loaded) {
        f.complete(loaded);
        CACHE_ORDER.add(new Cached(this, key, loaded.raw.length));
        long total = CACHED_BYTES.addAndGet(loaded.raw.length);
        while (total > CACHE_BYTES) {
            Cached old = CACHE_ORDER.poll();
            if (old == null) {
                break;
            }
            old.file.cache.remove(old.key);
            total = CACHED_BYTES.addAndGet(-old.bytes);
        }
    }

    private Loaded loadChunk(Group group, int i, byte[] bytes) throws IOException {
        if (bytes == null) {
            bytes = readAt(group.offsets[i], group.sizes[i]);
        }
        boolean stored = codec == 0 || (group.kinds[i] & 0x80) != 0;
        byte[] raw;
        if (stored) {
            if (codec == 1) {
                check(crc32(bytes, 0, bytes.length) == group.crcs[i], "a stored chunk's CRC-32 does not match");
            }
            check(bytes.length == group.raws[i], "a stored chunk's size is not its raw size");
            raw = bytes;
        } else {
            raw = inflate(bytes, group.raws[i]);
        }
        int kind = group.kinds[i] & 7;
        int code = (group.kinds[i] >> 3) & 15;
        int tracks = group.nTracks;
        check(kind <= 4 && code <= 8, "a chunk's kind byte is unknown");
        if (kind == 1) {
            check(code != 7 && raw.length >= (long) tracks * tile * WIDTHS[code], "a dense chunk is cut short");
        } else if (kind == 3) {
            checkMatrix(raw, tracks, code);
        }
        if (kind != 2 && kind != 4) {
            return new Loaded(raw, null, null, 0);
        }
        check(raw.length >= 4 + 4L * tracks && !(kind == 4 && code == 7), "a sparse or runs chunk is corrupt");
        // Sparse and runs: where each slot's offsets and values start, found once for every track.
        ByteBuffer b = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        int[] slotByte = new int[tracks + 1];
        int[] slotEntry = new int[tracks + 1];
        int p = 4 + 4 * tracks;
        int entries = 0;
        for (int j = 0; j < tracks; j++) {
            int length = b.getInt(4 + 4 * j);
            check(length >= 0 && (long) p + length <= raw.length, "a chunk's offsets pass its end");
            slotByte[j] = p;
            slotEntry[j] = entries;
            if (tracks > 1) {
                for (int q = p; q < p + length; q++) {
                    if (raw[q] >= 0) {
                        entries++;
                    }
                }
            }
            p += length;
        }
        slotByte[tracks] = p;
        check(tracks == 1 || entries == b.getInt(0), "a chunk's offsets do not add up to its total");
        slotEntry[tracks] = tracks > 1 ? entries : b.getInt(0);
        int valuesStart = (p + 3) & ~3;
        check(slotEntry[tracks] >= 0 && valuesStart + (long) slotEntry[tracks] * WIDTHS[code] <= raw.length,
                "a sparse or runs chunk's values pass its end");
        return new Loaded(raw, slotByte, slotEntry, valuesStart);
    }



    /** The blocks of a track-order group over [start, end) of a chromosome: [lo, k1) of its records. */
    private int[] blockRange(Group group, int chrom, int start, int end) {
        int n = group.bChrom.length;
        // first block on chrom whose last >= start
        int lo = 0, hi = n;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (group.bChrom[mid] < chrom || (group.bChrom[mid] == chrom && group.bLast[mid] < start)) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        int k1 = lo;
        while (k1 < n && group.bChrom[k1] == chrom && group.bFirst[k1] < end) {
            k1++;
        }
        return new int[]{lo, k1};
    }

    private Loaded[] loadBlocks(Group group, int track, int lo, int k1) throws IOException {
        return load(lo, k1, i -> (1L << 62) | ((long) track << 32) | i, group.bOffset, group.bSize,
                (i, bytes) -> {
                    if (bytes == null) {
                        bytes = readAt(group.bOffset[i], group.bSize[i]);
                    }
                    byte[] raw = codec == 0 ? bytes : inflate(bytes, group.bRaw[i]);
                    check(raw.length == group.bRaw[i], "a block's size is not its raw size");
                    int m = group.bN[i];
                    int code = (group.bKind[i] >> 3) & 15;
                    check(m > 0 && code <= 8 && raw.length >= ((4L * m + 7) & ~7L) + (long) m * WIDTHS[code],
                            "a block is cut short");
                    long prev = -1;
                    for (int e = 0; e < m; e++) {
                        long pos = (int) I32.get(raw, 4 * e) & 0xffffffffL;
                        check(pos > prev, "a block's positions do not ascend");
                        prev = pos;
                    }
                    check(((int) I32.get(raw, 0) & 0xffffffffL) == (group.bFirst[i] & 0xffffffffL)
                            && prev == (group.bLast[i] & 0xffffffffL), "a block's first or last position is not its record's");
                    return new Loaded(raw, null, null, 0);
                });
    }

    private void readBlocks(int track, int chrom, int start, int end, Runs runs) throws IOException {
        Group group = group(trackGroup[track]);
        int[] range = blockRange(group, chrom, start, end);
        int lo = range[0], k1 = range[1];
        Loaded[] loaded = loadBlocks(group, track, lo, k1);
        for (int i = lo; i < k1; i++) {
            byte[] raw = loaded[i - lo].raw;
            byte[] b = raw;
            int m = group.bN[i];
            int code = (group.bKind[i] >> 3) & 15;
            int w = WIDTHS[code];
            int valuesStart = (4 * m + 7) & ~7;
            // first entry at or past start
            int a = 0, z = m;
            while (a < z) {
                int mid = (a + z) >>> 1;
                if (((int) I32.get(b, 4 * mid) & 0xffffffffL) < start) {
                    a = mid + 1;
                } else {
                    z = mid;
                }
            }
            for (int e = a; e < m; e++) {
                int pos = (int) I32.get(b, 4 * e);
                if (pos >= end) {
                    break;
                }
                runs.add(pos, pos + 1, code == 7 ? 1.0 : value(b, code, valuesStart + e * w));
            }
        }
    }


    // ------------------------------------------------------------------------------------------
    // Summary and zoom levels

    public int levels() {
        try {
            zoomGrid();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        return zoomFactors.length;
    }

    /**
     * Bases per bin of a level; level 0 is the tile sums.
     */
    public long levelBases(int level) {
        levels();
        return (long) tile * zoomFactors[level];
    }

    /**
     * Over HTTP, the same level and range for the file's other groups with tracks in use, on background
     * threads: their tracks ask next, and the round trips then overlap instead of following one another.
     */
    private void prefetchLevel(int g, int level, int chrom, int start, int end) {
        String range = level + ":" + chrom + ":" + start + ":" + end;
        for (int x = 0; x < nGroups; x++) {
            if (x != g && groupUsers.get(x) > 0 && !range.equals(levelPrefetched.put(x, range))) {
                final int track = groupTrack[x];
                PREFETCHER.submit(() -> {
                    try {
                        levelMeans(track, level, chrom, start, end);
                    } catch (RuntimeException ignored) {
                        // The group's own reader reads the level again and reports the error.
                    }
                });
            }
        }
    }

    private final ConcurrentHashMap<Integer, String> levelPrefetched = new ConcurrentHashMap<>();

    // Over HTTP a level is read by windows only from this size: below it, one request for the whole
    // level costs less than the two for its pages and fragments.
    private static final long HTTP_WINDOWED_LEVEL_BYTES = Long.getLong("figwig.httpWindowedLevelBytes", 1L << 20);

    /**
     * The rows of a level over [start, end) of a chromosome from the regional summary index: only the
     * windows of bins the range covers.  Null when the level was not restarted, or the range covers half
     * its windows or more, where reading the level whole costs less.
     */
    private Level levelWindows(Group group, int k, int chrom, int start, int end) throws IOException {
        if (group.catalog == null || !group.catalog[k].restarted()) {
            return null;
        }
        FigWigSummaryIndex.LevelEntry e = group.catalog[k];
        long binBases = levelBases(k);
        long first = chromBinStart[k][chrom] + Math.max(0, start) / binBases;
        long last = chromBinStart[k][chrom] + (Math.min(end, chromLengths[chrom]) - 1) / binBases;
        if (last < first) {
            return new Level(new int[0], new byte[0], 0, k == 0 ? 4 : 8);
        }
        long j0 = first / e.windowBins, j1 = last / e.windowBins;
        if (2 * (j1 - j0 + 1) >= e.windows() || (HTTP && e.sourceBytes < HTTP_WINDOWED_LEVEL_BYTES)) {
            return null;
        }
        FigWigSummaryIndex.Window[] ws = summaryIndex.windows(group.index, group.catalog, k, group.nTracks, j0, j1);
        int rows = 0, bytes = 0;
        for (FigWigSummaryIndex.Window w : ws) {
            rows += w.bins.length;
            bytes += w.sums.length;
        }
        int[] bins = new int[rows];
        byte[] sums = new byte[bytes];
        rows = 0;
        bytes = 0;
        for (FigWigSummaryIndex.Window w : ws) {
            System.arraycopy(w.bins, 0, bins, rows, w.bins.length);
            System.arraycopy(w.sums, 0, sums, bytes, w.sums.length);
            rows += w.bins.length;
            bytes += w.sums.length;
        }
        return new Level(bins, sums, 0, k == 0 ? 4 : 8);
    }

    /** Whether the track's group keeps a summary (a writer leaves it out of groups of many tracks). */
    public boolean hasSummary(int track) {
        try {
            return group(trackGroup[track]).summaryBytes > 0;
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * The bins of a level over [start, end) that hold a value, with each bin's mean over its bases
     * within the chromosome (bases without a value count as 0); null when the track's group keeps no
     * summary.
     */
    public Runs levelMeans(int track, int level, int chrom, int start, int end) {
        Group group;
        Level data;
        try {
            zoomGrid();
            group = group(trackGroup[track]);
            if (group.summaryBytes <= 0) {
                return null;
            }
            check(Arrays.equals(zooms(group), zoomFactors), "groups have different zoom levels");
            data = group.levelData.get(level);
            if (data == null && HTTP) {
                prefetchLevel(group.index, level, chrom, start, end);
            }
            if (data == null) {
                data = levelWindows(group, level, chrom, start, end);
            }
            if (data == null) {
                data = level(group, level);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        int slot = trackSlot[track];
        long binBases = levelBases(level);
        long first = chromBinStart[level][chrom] + start / binBases;
        long last = chromBinStart[level][chrom] + (Math.min(end, chromLengths[chrom]) - 1) / binBases;
        int[] bins = data.bins;
        Runs runs = new Runs();
        for (int i = lowerBound(bins, first); i < bins.length && bins[i] <= last; i++) {
            long binStart = (bins[i] - chromBinStart[level][chrom]) * binBases;
            long binEnd = Math.min(binStart + binBases, chromLengths[chrom]);
            double sum = data.sum(i * group.nTracks + slot);
            runs.add((int) binStart, (int) binEnd, sum / (binEnd - binStart));
        }
        return runs;
    }

    // ------------------------------------------------------------------------------------------
    // Bytes

    public long chunkBytesRead() {
        return chunkBytesRead.get();
    }

    public int chunkReads() {
        return (int) chunkReads.get();
    }

    public void clearCache() {
        cache.clear();
        synchronized (batches) {
            batches.clear();
            batchRuns.clear();
            heldRuns = 0;
        }
        long[] freed = {0};
        CACHE_ORDER.removeIf(c -> {
            if (c.file == this) {
                freed[0] += c.bytes;
                return true;
            }
            return false;
        });
        CACHED_BYTES.addAndGet(-freed[0]);
    }

    /** A track now reads this file. */
    public void retain() {
        users.incrementAndGet();
    }

    /**
     * Track {@code track} now reads this file.  Its group's index record is read and parsed at once in
     * the background, so that the first view does not wait for it and groups are parsed in parallel.
     */
    public void retain(int track) {
        retain();
        trackUsers.incrementAndGet(track);
        int g = trackGroup[track];
        ExecutorService pool = HTTP ? PREFETCHER : BACKGROUND_PARSE ? DECODER : null;
        if (groupUsers.getAndIncrement(g) == 0 && pool != null && groups.get(g) == null) {
            pool.submit(() -> {
                try {
                    Group group = group(g);
                    if (summaryIndex != null && group.summaryBytes > 0) {
                        zooms(group);   // the group's level catalog: small, and the first wide view needs it
                    }
                } catch (Exception ignored) {
                    // The first read of the group parses it again and reports the error.
                }
            });
        }
    }

    private static final boolean BACKGROUND_PARSE = !"false".equals(System.getProperty("figwig.backgroundParse"));

    /** A track no longer reads this file; the last one to let go closes it. */
    public void release() throws IOException {
        if (users.decrementAndGet() == 0) {
            close();
        }
    }

    public void release(int track) throws IOException {
        groupUsers.decrementAndGet(trackGroup[track]);
        trackUsers.decrementAndGet(track);
        release();
    }

    /** Read over HTTP range requests (http, https, s3, gs paths). */
    public boolean isHttp() {
        return HTTP;
    }

    public boolean isClosed() {
        return closed;
    }

    /**
     * A local file is read with positional reads, which share no file position and take no lock;
     * anything else (HTTP, S3) through IGV's SeekableStream under one lock.
     */
    byte[] readAt(long offset, int n) throws IOException {
        byte[] out = new byte[n];
        if (channel != null) {
            ByteBuffer b = ByteBuffer.wrap(out);
            long at = offset;
            while (b.hasRemaining()) {
                int r = channel.read(b, at);
                if (r < 0) {
                    throw new IOException("figwig: " + path + " ended at " + at + ", expected " + (offset + n));
                }
                at += r;
            }
        } else if (HTTP) {
            // One range request per read on a stream of its own, so IGV's loader threads read concurrently.
            try (SeekableStream s = IGVSeekableStreamFactory.getInstance().getStreamFor(path)) {
                fill(s, offset, out);
            }
        } else {
            synchronized (stream) {
                fill(stream, offset, out);
            }
        }
        chunkReads.incrementAndGet();
        chunkBytesRead.addAndGet(n);
        return out;
    }

    private void fill(SeekableStream s, long offset, byte[] out) throws IOException {
        s.seek(offset);
        int done = 0;
        while (done < out.length) {
            int r = s.read(out, done, out.length - done);
            if (r < 0) {
                throw new IOException("figwig: " + path + " ended at " + (offset + done) + ", expected " + (offset + out.length));
            }
            done += r;
        }
    }

    public void close() throws IOException {
        closed = true;
        clearCache();
        if (summaryIndex != null) {
            summaryIndex.clear();
        }
        stream.close();
        if (channel != null) {
            channel.close();
        }
    }



    private static int lowerBound(int[] a, long x) {
        int lo = 0, hi = a.length;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if ((a[mid] & 0xffffffffL) < x) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        return lo;
    }

}
