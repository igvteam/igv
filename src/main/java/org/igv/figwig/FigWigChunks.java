package org.igv.figwig;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * The bytes of figwig chunks and blocks, apart from any file: the element types, inflating, the checks
 * of a matrix chunk, and the decoding of one track's values from a chunk into runs.
 */
final class FigWigChunks {

    private FigWigChunks() {
    }

    static final int[] WIDTHS = {1, 2, 4, 1, 2, 4, 4, 0, 8};
    static final VarHandle I16 =
            MethodHandles.byteArrayViewVarHandle(short[].class, ByteOrder.LITTLE_ENDIAN);
    static final VarHandle I32 =
            MethodHandles.byteArrayViewVarHandle(int[].class, ByteOrder.LITTLE_ENDIAN);
    static final VarHandle F32 =
            MethodHandles.byteArrayViewVarHandle(float[].class, ByteOrder.LITTLE_ENDIAN);
    static final VarHandle F64 =
            MethodHandles.byteArrayViewVarHandle(double[].class, ByteOrder.LITTLE_ENDIAN);
    /**
     * A chunk or block as decoded once for every track that reads it: its bytes, and for a sparse or
     * runs chunk where each slot's offsets and values start.
     */
    static final class Loaded {
        final byte[] raw;
        final int[] slotByte;    // sparse and runs: the first byte of each slot's offsets, and the end
        final int[] slotEntry;   // the index of each slot's first value, and the total
        final int valuesStart;

        Loaded(byte[] raw, int[] slotByte, int[] slotEntry, int valuesStart) {
            this.raw = raw;
            this.slotByte = slotByte;
            this.slotEntry = slotEntry;
            this.valuesStart = valuesStart;
        }
    }
    /** A matrix chunk's header, pointers and track numbers, before any of its offsets are read. */
    static void checkMatrix(byte[] raw, int tracks, int code) {
        check(raw.length >= 4, "a matrix chunk is cut short");
        long nPresent = (int) I32.get(raw, 0) & 0xffffffffL;
        check(nPresent <= tracks && raw.length >= 8 + 8 * nPresent, "a matrix chunk is corrupt");
        int n = (int) nPresent;
        int ptr = 4 + 4 * n;
        long total = (int) I32.get(raw, ptr + 4 * n) & 0xffffffffL;
        long valuesAt = (8 + 8L * n + 2 * total + 3) & ~3L;
        check((int) I32.get(raw, ptr) == 0 && valuesAt + total * WIDTHS[code] <= raw.length, "a matrix chunk is corrupt");
        for (int j = 0; j < n; j++) {
            long t = (int) I32.get(raw, 4 + 4 * j) & 0xffffffffL;
            long lo = (int) I32.get(raw, ptr + 4 * j) & 0xffffffffL;
            long hi = (int) I32.get(raw, ptr + 4 * (j + 1)) & 0xffffffffL;
            check(t < tracks && (j == 0 || t > ((int) I32.get(raw, 4 * j) & 0xffffffffL)) && lo <= hi && hi <= total,
                    "a matrix chunk's tracks or pointers are out of order");
        }
    }
    static void decodeChunk(Loaded loaded, int kindByte, int tracks, int slot, int tileBases, int tileStart,
                                    int tileEnd, int from, int to, FigWigFile.Runs runs) {
        byte[] b = loaded.raw;
        int kind = kindByte & 7;
        int code = (kindByte >> 3) & 15;
        if (kind > 4 || code > 8) {
            throw new IllegalStateException("figwig: kind byte " + kindByte);
        }
        switch (kind) {
            case 0:
                return;
            case 1: {  // dense
                check(code != 7, "a dense chunk with the one code");
                int w = WIDTHS[code];
                int at = slot * tileBases * w + (from - tileStart) * w;
                switch (code) {
                    case 0:
                        for (int p = from; p < to; p++, at++) runs.add(p, p + 1, b[at] & 0xff);
                        break;
                    case 1:
                        for (int p = from; p < to; p++, at += 2) runs.add(p, p + 1, (short) I16.get(b, at) & 0xffff);
                        break;
                    case 6:
                        for (int p = from; p < to; p++, at += 4) runs.add(p, p + 1, (float) F32.get(b, at));
                        break;
                    default:
                        for (int p = from; p < to; p++, at += w) runs.add(p, p + 1, value(b, code, at));
                }
                return;
            }
            case 2:
            case 4: {  // sparse, runs
                check(!(kind == 4 && code == 7), "a runs chunk with the one code");
                int total = loaded.slotEntry[tracks];
                int mine = loaded.slotByte[slot];
                int myLength = loaded.slotByte[slot + 1] - mine;
                int before = loaded.slotEntry[slot];
                int w = WIDTHS[code];
                int valuesStart = loaded.valuesStart;
                int idx = before;
                int p = mine;
                int end = mine + myLength;
                int offset = -1;
                int pendingStart = -1;
                double pendingValue = 0;
                while (p < end) {
                    int v = 0;
                    int shift = 0;
                    int c;
                    do {
                        if (p >= end || shift > 14) {
                            throw new IllegalStateException("figwig: a chunk's offsets are corrupt");
                        }
                        c = b[p++] & 0xff;
                        v |= (c & 0x7f) << shift;
                        shift += 7;
                    } while ((c & 0x80) != 0);
                    offset = offset < 0 ? v : offset + v;
                    if (offset >= tileBases) {
                        throw new IllegalStateException("figwig: a chunk's offset passes its tile");
                    }
                    double value;
                    if (code == 6) {
                        value = (float) F32.get(b, valuesStart + 4 * idx);
                    } else if (code == 7) {
                        value = 1.0;
                    } else if (code == 0) {
                        value = b[valuesStart + idx] & 0xff;
                    } else {
                        value = value(b, code, valuesStart + idx * w);
                    }
                    idx++;
                    int pos = tileStart + offset;
                    if (kind == 2) {
                        if (pos >= from && pos < to) {
                            runs.add(pos, pos + 1, value);
                        }
                    } else {
                        if (pendingStart >= 0) {
                            runs.add(Math.max(pendingStart, from), Math.min(pos, to), pendingValue);
                        }
                        pendingStart = pos;
                        pendingValue = value;
                    }
                }
                if (kind == 4 && pendingStart >= 0) {
                    runs.add(Math.max(pendingStart, from), Math.min(tileEnd, to), pendingValue);
                }
                check(tracks == 1 ? idx == total : idx == loaded.slotEntry[slot + 1], "a chunk's offsets do not match its total");
                return;
            }
            case 3: {  // matrix
                int nPresent = (int) I32.get(b, 0);
                int j = -1;
                int lo = 0, hi = nPresent - 1;
                while (lo <= hi) {
                    int mid = (lo + hi) >>> 1;
                    int s = (int) I32.get(b, 4 + 4 * mid);
                    if (s < slot) {
                        lo = mid + 1;
                    } else if (s > slot) {
                        hi = mid - 1;
                    } else {
                        j = mid;
                        break;
                    }
                }
                if (j < 0) {
                    return;
                }
                int ptr = 4 + 4 * nPresent;
                int total = (int) I32.get(b, ptr + 4 * nPresent);
                int a = (int) I32.get(b, ptr + 4 * j);
                int z = (int) I32.get(b, ptr + 4 * (j + 1));
                int offsets = 8 + 8 * nPresent;
                int w = WIDTHS[code];
                int valuesStart = (offsets + 2 * total + 3) & ~3;
                for (int i = a; i < z; i++) {
                    int off = (short) I16.get(b, offsets + 2 * i) & 0xffff;
                    if (off >= tileBases || (i > a && off <= ((short) I16.get(b, offsets + 2 * (i - 1)) & 0xffff))) {
                        throw new IllegalStateException("figwig: a matrix chunk's offsets are out of order");
                    }
                    int pos = tileStart + off;
                    if (pos >= from && pos < to) {
                        runs.add(pos, pos + 1, code == 7 ? 1.0 : value(b, code, valuesStart + i * w));
                    }
                }
                return;
            }
            default:
                throw new IllegalStateException("figwig: kind " + kind);
        }
    }
    static double value(byte[] b, int code, int at) {
        switch (code) {
            case 0:
                return b[at] & 0xff;
            case 1:
                return (short) I16.get(b, at) & 0xffff;
            case 2:
                return (int) I32.get(b, at) & 0xffffffffL;
            case 3:
                return b[at];
            case 4:
                return (short) I16.get(b, at);
            case 5:
                return (int) I32.get(b, at);
            case 6:
                return (float) F32.get(b, at);
            case 7:
                return 1.0;
            case 8:
                return (double) F64.get(b, at);
            default:
                throw new IllegalStateException("figwig: element code " + code);
        }
    }
    /** One Inflater per thread, reset for each chunk: creating and ending one per chunk costs native setup. */
    static final ThreadLocal<Inflater> INFLATERS = ThreadLocal.withInitial(Inflater::new);

    static byte[] inflate(byte[] in, int raw) throws IOException {
        Inflater inflater = INFLATERS.get();
        inflater.reset();
        try {
            inflater.setInput(in);
            byte[] out = new byte[raw];
            int done = 0;
            while (done < raw && !inflater.finished()) {
                int r = inflater.inflate(out, done, raw - done);
                if (r == 0 && (inflater.needsInput() || inflater.needsDictionary())) {
                    break;
                }
                done += r;
            }
            if (done != raw) {
                throw new IllegalStateException("figwig: a chunk inflated to " + done + " bytes, not " + raw);
            }
            if (!inflater.finished()) {
                byte[] extra = new byte[1];
                check(inflater.inflate(extra) == 0 && inflater.finished(), "a chunk inflates past its raw size");
            }
            check(inflater.getRemaining() == 0, "bytes after a chunk's zlib stream");
            return out;
        } catch (DataFormatException e) {
            throw new IOException("figwig: a corrupt chunk: " + e.getMessage(), e);
        }
    }
    static long crc32(byte[] b, int off, int len) {
        CRC32 crc = new CRC32();
        crc.update(b, off, len);
        return crc.getValue();
    }
    static void check(boolean ok, String message) {
        if (!ok) {
            throw new IllegalStateException("figwig: " + message);
        }
    }
}
