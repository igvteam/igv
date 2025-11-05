package org.broad.igv.hic;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;
import org.broad.igv.ucsc.twobit.UnsignedByteBufferImpl;
import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.collections.LRUCache;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.util.*;

/**
 * Partial conversion of the JavaScript HicFile to Java.
 * - Expects a FileChannel positioned on the underlying .hic file.
 * - Relies on existing UnsignedByteBuffer, Matrix, ContactRecord, NormalizationVector, LRUCache classes.
 *
 * Note: Remote/file-provider abstractions, rate limiting, and browser-specific behavior are omitted.
 * This class provides synchronous IO methods (FileChannel reads). Methods retain original names but throw IOException.
 */
public class HicFile {

    private final SeekableStream fileChannel;
    private final Map<String, Object> config;

    private boolean initialized = false;
    private String magic;
    private Integer version;
    private long footerPosition;

    private String genomeId;
    private long normVectorIndexPosition;
    private long normVectorIndexSize;
    private Map<String, IndexEntry> masterIndex;
    private Map<String, Long> expectedValueVectors;
    private Map<String, Object> attributes;
    private List<Chromosome> chromosomes = new ArrayList<>();
    private Map<String, Integer> chromosomeIndexMap = new HashMap<>();
    private List<Integer> bpResolutions = new ArrayList<>();
    private List<Integer> fragResolutions = new ArrayList<>();
    private Map<String, String> chrAliasTable = new HashMap<>();
    private Map<String, IndexEntry> normVectorIndex;
    private LRUCache<String, NormalizationVector> normVectorCache = new LRUCache<>(10);
    private LRUCache<String, Matrix> matrixCache = new LRUCache<>(10);
    private BlockCache blockCache = new BlockCache();
    private List<String> normalizationTypes = new ArrayList<>(Collections.singletonList("NONE"));
    private Long normExpectedValueVectorsPosition;

    public static HicFile create(String path) throws IOException {
        return create(path, Collections.emptyMap());
    }

    public static HicFile create(String path, Map<String, Object> config) throws IOException {
        SeekableStream stream = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        HicFile hicFile = new HicFile(stream, config);
        hicFile.init();
        return hicFile;
    }

    public HicFile(SeekableStream fileChannel, Map<String, Object> config) {
        this.fileChannel = fileChannel;
        this.config = config != null ? config : Collections.emptyMap();

    }

    public synchronized void init() throws IOException {
        if (initialized) return;
        readHeaderAndFooter();
        initialized = true;
    }

    public int getVersion() {
        return this.version;
    }

    private void readHeaderAndFooter() throws IOException {
        // Read initial fields magic, version, and footer position
        byte[] header = readBytes(0, 16);
        if (header == null || header.length == 0) throw new IOException("File content is empty");

        UnsignedByteBuffer bp = UnsignedByteBufferImpl.wrap(header);
        this.magic = bp.getString();
        this.version = bp.getInt();
        if (this.version < 5) throw new IOException("Unsupported hic version: " + this.version);
        this.footerPosition = bp.getLong();

        readFooter();

        // determine body position (min start of master index entries)
        long bodyPosition = Long.MAX_VALUE;
        for (IndexEntry e : masterIndex.values()) {
            bodyPosition = Math.min(bodyPosition, e.start);
        }
        long remainingSize = bodyPosition - 16;
        byte[] body = readBytes(16, (int) remainingSize);
        UnsignedByteBuffer bodyParser = UnsignedByteBufferImpl.wrap(body);

        this.genomeId = bodyParser.getString();

        if (this.version >= 9) {
            this.normVectorIndexPosition = bodyParser.getLong();
            this.normVectorIndexSize = bodyParser.getLong();
        }

        // attributes
        this.attributes = new HashMap<>();
        int nAttributes = bodyParser.getInt();
        while (nAttributes-- > 0) {
            String k = bodyParser.getString();
            String v = bodyParser.getString();
            this.attributes.put(k, v);
        }

        // chromosomes
        this.chromosomes = new ArrayList<>();
        this.chromosomeIndexMap = new HashMap<>();
        int nChrs = bodyParser.getInt();
        for (int i = 0; i < nChrs; i++) {
            String name = bodyParser.getString();
            long size = this.version < 9 ? bodyParser.getInt() : bodyParser.getLong();
            Chromosome chr = new Chromosome(i, name, (int) size);
            if ("all".equalsIgnoreCase(name)) {
                // whole genome handling omitted other fields
            }
            this.chromosomes.add(chr);
            this.chromosomeIndexMap.put(name, i);
        }

        // bp resolutions
        int nBp = bodyParser.getInt();
        for (int i = 0; i < nBp; i++) {
            this.bpResolutions.add(bodyParser.getInt());
        }

        // frag resolutions (optional)
        if (Boolean.TRUE.equals(config.get("loadFragData"))) {
            int nFrag = bodyParser.getInt();
            for (int i = 0; i < nFrag; i++) {
                this.fragResolutions.add(bodyParser.getInt());
            }
        }

        // build alias table
        for (String chrName : chromosomeIndexMap.keySet()) {
            if (chrName.startsWith("chr")) chrAliasTable.put(chrName.substring(3), chrName);
            else if ("MT".equals(chrName)) chrAliasTable.put("chrM", chrName);
            else chrAliasTable.put("chr" + chrName, chrName);
        }

        // meta (not stored in separate object here)
    }

    private void readFooter() throws IOException {
        int skip = this.version < 9 ? 8 : 12;
        byte[] data = readBytes(this.footerPosition, skip);
        if (data == null) return;
        UnsignedByteBuffer bp = UnsignedByteBufferImpl.wrap(data);
        long nBytes = this.version < 9 ? bp.getInt() : bp.getLong();
        int nEntries = bp.getInt();

        int miSize = nEntries * (100 + 64 + 32);
        byte[] miData = readBytes(this.footerPosition + skip, Math.min(miSize, (int) nBytes));
        UnsignedByteBuffer miParser = UnsignedByteBufferImpl.wrap(miData);

        this.masterIndex = new HashMap<>();
        while (nEntries-- > 0) {
            String key = miParser.getString();
            long pos = miParser.getLong();
            int size = miParser.getInt();
            masterIndex.put(key, new IndexEntry(pos, size));
        }

        // expected values start after master index; compute normExpectedValueVectorsPosition
        if (this.version > 5) {
            int skip2 = this.version < 9 ? 4 : 8;
            this.normExpectedValueVectorsPosition = this.footerPosition + skip2 + nBytes;
        }
    }

    public Matrix getMatrix(int chrIdx1, int chrIdx2) throws IOException {
        String key = Matrix.getKey(chrIdx1, chrIdx2);
        if (matrixCache.containsKey(key)) return matrixCache.get(key);
        Matrix m = readMatrix(chrIdx1, chrIdx2);
        if (m != null) matrixCache.put(key, m);
        return m;
    }

    private Matrix readMatrix(int chrIdx1, int chrIdx2) throws IOException {
        init();
        if (chrIdx1 > chrIdx2) {
            int tmp = chrIdx1;
            chrIdx1 = chrIdx2;
            chrIdx2 = tmp;
        }
        String key = Matrix.getKey(chrIdx1, chrIdx2);
        IndexEntry idx = masterIndex.get(key);
        if (idx == null) return null;
        byte[] data = readBytes(idx.start, idx.size);
        if (data == null) return null;
        // Matrix.parseMatrix expects byte[] and chromosome list
        return Matrix.parseMatrix(data, this.chromosomes);
    }

    public List<Integer> getBpResolutions() {
        return bpResolutions;
    }

    public List<ContactRecord> getContactRecords(String normalization,
                                                 Region region1,
                                                 Region region2,
                                                 String units,
                                                 int binSize,
                                                 boolean allRecords) throws IOException {

        init();

        int idx1 = chromosomeIndexMap.getOrDefault(getFileChrName(region1.chr()), -1);
        int idx2 = chromosomeIndexMap.getOrDefault(getFileChrName(region2.chr()), -1);
        boolean transpose = idx1 > idx2 || (idx1 == idx2 && region1.start() >= region2.end());

        if (transpose) {
            Region tmp = region1;
            region1 = region2;
            region2 = tmp;
        }

        List<Block> blocks = getBlocks(region1, region2, units, binSize);
        if (blocks == null || blocks.isEmpty()) return Collections.emptyList();

        List<ContactRecord> contactRecords = new ArrayList<>();
        double x1 = (double) region1.start() / binSize;
        double x2 = (double) region1.end() / binSize;
        double y1 = (double) region2.start() / binSize;
        double y2 = (double) region2.end() / binSize;
        int nvX1 = (int) Math.floor(x1);
        int nvX2 = (int) Math.ceil(x2);
        int nvY1 = (int) Math.floor(y1);
        int nvY2 = (int) Math.ceil(y2);

        boolean isNorm = normalization != null && !"NONE".equals(normalization);
        String chr1 = getFileChrName(region1.chr());
        String chr2 = getFileChrName(region2.chr());

        for (Block block : blocks) {
            if (block == null) continue;

            NormalizationVector nv1 = null;
            NormalizationVector nv2 = null;
            if (isNorm) {
                nv1 = getNormalizationVector(normalization, chr1, units, binSize);
                nv2 = chr1.equals(chr2) ? nv1 : getNormalizationVector(normalization, chr2, units, binSize);
                if (nv1 == null || nv2 == null) {
                    isNorm = false;
                }
            }

            for (ContactRecord rec : block.records) {
                if (allRecords || (rec.bin1() >= x1 && rec.bin1() < x2 && rec.bin2() >= y1 && rec.bin2() < y2)) {
                    if (isNorm) {
                        int x = rec.bin1();
                        int y = rec.bin2();
                        double[] v1 = nv1.getValues(nvX1, nvX2);
                        double[] v2 = nv2.getValues(nvY1, nvY2);
                        double nvnv = v1[x - nvX1] * v2[y - nvY1];
                        if (nvnv != 0 && !Double.isNaN(nvnv)) {
                            double counts = rec.counts() / nvnv;
                            contactRecords.add(new ContactRecord(x, y, counts));
                        }
                    } else {
                        contactRecords.add(rec);
                    }
                }
            }
        }

        return contactRecords;
    }

    public List<Block> getBlocks(Region region1, Region region2, String unit, int binSize) throws IOException {
        init();
        String chr1 = getFileChrName(region1.chr());
        String chr2 = getFileChrName(region2.chr());
        Integer idx1 = chromosomeIndexMap.get(chr1);
        Integer idx2 = chromosomeIndexMap.get(chr2);
        if (idx1 == null || idx2 == null) return Collections.emptyList();

        Matrix matrix = getMatrix(idx1, idx2);
        if (matrix == null) return Collections.emptyList();

        MatrixZoomData zd = matrix.getZoomData(binSize, unit);
        if (zd == null) throw new IOException("No data available for resolution: " + binSize);

        List<Integer> blockNumbers = zd.getBlockNumbers(region1, region2, this.version);
        List<Block> blocks = new ArrayList<>();
        List<Integer> toQuery = new ArrayList<>();
        for (Integer num : blockNumbers) {
            String key = zd.getKey() + "_" + num;
            if (blockCache.has(binSize, key)) {
                blocks.add(blockCache.get(binSize, key));
            } else {
                toQuery.add(num);
            }
        }

        for (Integer bn : toQuery) {
            Block b = readBlock(bn, zd);
            if (b != null) {
                blockCache.set(binSize, zd.getKey() + "_" + b.blockNumber, b);
            }
            blocks.add(b);
        }

        return blocks;
    }

    public Block readBlock(int blockNumber, MatrixZoomData zd) throws IOException {
        StaticBlockIndex.BlockIndexEntry idx = zd.getBlockIndex().getBlockIndexEntry(blockNumber);
        if (idx == null) return null;
        byte[] data = readBytes(idx.filePosition, idx.size);
        if (data == null) return null;

        // decompress
        byte[] plain = (new CompressionUtils()).decompress(data);

        UnsignedByteBuffer parser = UnsignedByteBufferImpl.wrap(plain);
        int nRecords = parser.getInt();
        List<ContactRecord> records = new ArrayList<>();

        if (this.version < 7) {
            for (int i = 0; i < nRecords; i++) {
                int binX = parser.getInt();
                int binY = parser.getInt();
                double counts = parser.getFloat();
                records.add(new ContactRecord(binX, binY, counts));
            }
        } else {
            int binXOffset = parser.getInt();
            int binYOffset = parser.getInt();
            boolean useFloatContact = parser.get() == 1;
            boolean useIntXPos = this.version < 9 ? false : parser.get() == 1;
            boolean useIntYPos = this.version < 9 ? false : parser.get() == 1;
            int type = parser.get();

            final int Short_MIN_VALUE = -32768;

            if (type == 1) {
                int rowCount = useIntYPos ? parser.getInt() : parser.getShort();
                for (int i = 0; i < rowCount; i++) {
                    int dy = useIntYPos ? parser.getInt() : parser.getShort();
                    int binY = binYOffset + dy;
                    int colCount = useIntXPos ? parser.getInt() : parser.getShort();
                    for (int j = 0; j < colCount; j++) {
                        int dx = useIntXPos ? parser.getInt() : parser.getShort();
                        int binX = binXOffset + dx;
                        double counts = useFloatContact ? parser.getFloat() : parser.getShort();
                        records.add(new ContactRecord(binX, binY, counts));
                    }
                }
            } else if (type == 2) {
                int nPts = parser.getInt();
                int w = parser.getShort();
                for (int i = 0; i < nPts; i++) {
                    int row = i / w;
                    int col = i - row * w;
                    int bin1 = binXOffset + col;
                    int bin2 = binYOffset + row;
                    if (useFloatContact) {
                        double counts = parser.getFloat();
                        if (!Double.isNaN(counts)) records.add(new ContactRecord(bin1, bin2, counts));
                    } else {
                        int counts = parser.getShort();
                        if (counts != Short_MIN_VALUE) records.add(new ContactRecord(bin1, bin2, counts));
                    }
                }
            } else {
                throw new IOException("Unknown block type: " + type);
            }
        }

        return new Block(blockNumber, zd, records, idx);
    }

    public boolean hasNormalizationVector(String type, Object chr, String unit, int binSize) throws IOException {
        init();
        int chrIdx;
        if (chr instanceof Integer) chrIdx = (Integer) chr;
        else chrIdx = chromosomeIndexMap.getOrDefault(getFileChrName(chr.toString()), -1);
        String key = getNormalizationVectorKey(type, chrIdx, unit, binSize);
        Map<String, IndexEntry> nvi = getNormVectorIndex();
        return nvi != null && nvi.containsKey(key);
    }

    public NormalizationVector getNormalizationVector(String type, Object chr, String unit, int binSize) throws IOException {
        init();
        int chrIdx;
        if (chr instanceof Integer) chrIdx = (Integer) chr;
        else chrIdx = chromosomeIndexMap.getOrDefault(getFileChrName(chr.toString()), -1);
        String key = getNormalizationVectorKey(type, chrIdx, unit, binSize);
        if (normVectorCache.containsKey(key)) return normVectorCache.get(key);
        Map<String, IndexEntry> nvi = getNormVectorIndex();
        if (nvi == null) return null;
        IndexEntry idx = nvi.get(key);
        if (idx == null) return null;
        byte[] header = readBytes(idx.start(), 8);
        if (header == null) return null;
        UnsignedByteBuffer parser = UnsignedByteBufferImpl.wrap(header);
        long nValues = this.version < 9 ? parser.getInt() : parser.getLong();
        int dataType = this.version < 9 ? BinaryTypes.DOUBLE : BinaryTypes.FLOAT;
        long filePos = this.version < 9 ? idx.start() + 4 : idx.start() + 8;
        NormalizationVector nv = new NormalizationVector(fileChannel, filePos, (int) nValues, dataType);
        normVectorCache.put(key, nv);
        return nv;
    }

    public Map<String, IndexEntry> getNormVectorIndex() throws IOException {
        if (this.version < 6) return null;
        if (this.normVectorIndex == null) {
            // If config contains "nvi" string (start,size) use it
            if (config.containsKey("nvi")) {
                String nviStr = (String) config.get("nvi");
                String[] parts = nviStr.split(",");
                long start = Long.parseLong(parts[0]);
                int size = Integer.parseInt(parts[1]);
                readNormVectorIndex(new IndexEntry(start, size));
            } else {
                try {
                    readNormExpectedValuesAndNormVectorIndex();
                } catch (IOException e) {
                    // File may not contain norm vectors
                }
            }
        }
        return normVectorIndex;
    }

    public void readNormVectorIndex(IndexEntry range) throws IOException {
        init();
        byte[] data = readBytes(range.start, range.size);
        UnsignedByteBuffer bp = UnsignedByteBufferImpl.wrap(data);
        this.normVectorIndex = new HashMap<>();
        int nEntries = bp.getInt();
        while (nEntries-- > 0) {
            parseNormVectorEntry(bp);
        }
    }

    public void readNormExpectedValuesAndNormVectorIndex() throws IOException {
        init();
        if (this.normExpectedValueVectorsPosition == null) return;
        long nviStart = skipExpectedValues(this.normExpectedValueVectorsPosition);
        byte[] header = readBytes(nviStart, 4);
        if (header == null || header.length == 0) return;
        UnsignedByteBuffer bp = UnsignedByteBufferImpl.wrap(header);
        int nEntries = bp.getInt();
        int sizeEstimate = nEntries * 30;
        byte[] data = readBytes((int) (nviStart + 4), sizeEstimate);
        this.normVectorIndex = new HashMap<>();
        // process entries (simplified iterative parsing)
        UnsignedByteBuffer parser = UnsignedByteBufferImpl.wrap(data);
        while (parser.remaining() > 0) {
            parseNormVectorEntry(parser);
        }
        this.config.put("nvi", nviStart + "," + 4);
    }

    public long skipExpectedValues(long start) throws IOException {
        int versionLocal = this.version;
        // Use a Buffered read style - simplified
        byte[] header = readBytes(start, 4);
        UnsignedByteBuffer bp = UnsignedByteBufferImpl.wrap(header);
        int nEntries = bp.getInt();
        if (nEntries == 0) return start + 4;
        long p = start + 4;
        for (int i = 0; i < nEntries; i++) {
            byte[] small = readBytes(p, 500);
            UnsignedByteBuffer sp = UnsignedByteBufferImpl.wrap(small);
            String type = sp.getString();
            String unit = sp.getString();
            int binSize = sp.getInt();
            long nValues = versionLocal < 9 ? sp.getInt() : sp.getLong();
            long chunkSize = sp.position() + nValues * (versionLocal < 9 ? BinaryTypes.DOUBLE : BinaryTypes.FLOAT);
            byte[] after = readBytes((int) (p + chunkSize), 4);
            UnsignedByteBuffer ap = UnsignedByteBufferImpl.wrap(after);
            int nChrScaleFactors = ap.getInt();
            chunkSize += 4 + nChrScaleFactors * (4 + (versionLocal < 9 ? BinaryTypes.DOUBLE : BinaryTypes.FLOAT));
            p += chunkSize;
        }
        return p;
    }

    private void parseNormVectorEntry(UnsignedByteBuffer parser) {
        String type = parser.getString();
        int chrIdx = parser.getInt();
        String unit = parser.getString();
        int binSize = parser.getInt();
        long filePosition = parser.getLong();
        long sizeInBytes = this.version < 9 ? parser.getInt() : parser.getLong();
        String key = getNormalizationVectorKey(type, chrIdx, unit, binSize);
        if (!this.normalizationTypes.contains(type)) this.normalizationTypes.add(type);
        this.normVectorIndex.put(key, new IndexEntry(filePosition, (int) sizeInBytes));
    }

    private String getFileChrName(String chrAlias) {
        return chrAliasTable.getOrDefault(chrAlias, chrAlias);
    }

    public static String getNormalizationVectorKey(String type, int chrIdx, String unit, int resolution) {
        return type + "_" + chrIdx + "_" + unit + "_" + resolution;
    }

    private byte[] readBytes(long position, long size) throws IOException {
        if (size <= 0) return new byte[0];
        byte[] byteArray = new byte[(int) size];
        int read = 0;
        fileChannel.seek(position);
        while (read < size) {
            int r = fileChannel.read(byteArray, read, (int) (size - read));
            if (r < 0) break;
            read += r;
        }
        if (read == 0) return null;
        return byteArray;
    }

    private record IndexEntry(long start, int size) {
    }

    public static class Block {
        public final int blockNumber;
        public final MatrixZoomData zoomData;
        public final List<ContactRecord> records;
        public final StaticBlockIndex.BlockIndexEntry idx;

        public Block(int blockNumber, MatrixZoomData zd, List<ContactRecord> records, StaticBlockIndex.BlockIndexEntry idx) {
            this.blockNumber = blockNumber;
            this.zoomData = zd;
            this.records = records;
            this.idx = idx;
        }
    }

    public static class BlockCache {
        private Integer resolution = null;
        private LRUCache<String, Block> map = new LRUCache<>(6);

        public void set(int resolution, String key, Block value) {
            if (this.resolution == null || this.resolution != resolution) {
                this.map.clear();
                this.resolution = resolution;
            }
            this.map.put(key, value);
        }

        public Block get(int resolution, String key) {
            return (this.resolution != null && this.resolution == resolution) ? this.map.get(key) : null;
        }

        public boolean has(int resolution, String key) {
            return (this.resolution != null && this.resolution == resolution) && this.map.containsKey(key);
        }
    }

    // Placeholder utilities and classes used in this file (to be provided by project)
    private static class BinaryTypes {
        static final int DOUBLE = 8;
        static final int FLOAT = 4;
    }

}