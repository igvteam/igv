package org.igv.hic;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.igv.feature.Chromosome;
import org.igv.feature.genome.Genome;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.util.CompressionUtils;
import org.igv.util.collections.CaseInsensitiveMap;
import org.igv.util.collections.LRUCache;
import org.igv.util.stream.IGVSeekableStreamFactory;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;

/**
 * Representation of a .hic file, with methods to access contact recrods and normalization vectors.
 *
 * @author Jim Robinson
 */
public class HicFile {

    private static final Logger log = LogManager.getLogger(HicFile.class);

    enum DataType {
        FLOAT(4),
        DOUBLE(8);

        private final int byteSize;

        DataType(int byteSize) {
            this.byteSize = byteSize;
        }

        public int getByteSize() {
            return byteSize;
        }
    }

    private final SeekableStream fileChannel;
    private final Map<String, Object> config;
    private final String path;
    private final Genome genome;

    private boolean initialized = false;
    private String magic;
    private Integer version;
    private long footerPosition;

    private String genomeId;
    private long normVectorIndexPosition;
    private int normVectorIndexSize;
    private Map<String, IndexEntry> masterIndex;
    private Map<String, Long> expectedValueVectors;
    private Map<String, Object> attributes;
    private List<Chromosome> chromosomes = new ArrayList<>();
    private Map<String, Integer> chromosomeIndexMap = new CaseInsensitiveMap<Integer>();
    private Integer wgResolution = null;
    private List<Integer> bpResolutions = new ArrayList<>();
    private List<Integer> fragResolutions = new ArrayList<>();
    private Map<String, String> chrAliasTable = new HashMap<>();
    private Map<String, IndexEntry> normVectorIndex;
    private LRUCache<String, NormalizationVector> normVectorCache = new LRUCache<>(10);
    private LRUCache<String, Matrix> matrixCache = new LRUCache<>(10);
    private BlockCache blockCache = new BlockCache();
    private List<String> normalizationTypes = new ArrayList<>(Collections.singletonList("NONE"));
    private Long normExpectedValueVectorsPosition;

    // Map to track pending block read requests to prevent duplicate I/O operations
    // Key is the block key, value is a CompletableFuture that will complete when the block is loaded
    private final Map<String, CompletableFuture<Block>> pendingBlockRequests = new ConcurrentHashMap<>();


    public HicFile(String path, Genome genome) throws IOException {
        this.path = path;
        this.genome = genome;
        this.fileChannel = IGVSeekableStreamFactory.getInstance().getStreamFor(path);

        // TODO -- support customization options here.  This is an igv.js artifact.
        this.config = Collections.emptyMap();

        init();
    }

    private synchronized void init() throws IOException {
        if (initialized) return;
        readHeaderAndFooter();
        initialized = true;
    }

    public int getVersion() {
        return this.version;
    }

    public String getNVIString() {
        if (this.normVectorIndexPosition > 0 && this.normVectorIndexSize > 0) {
            return this.normVectorIndexPosition + "," + this.normVectorIndexSize;
        } else {
            return null;
        }
    }

    private static ByteBuffer wrap(byte[] data) {
        ByteBuffer bb = ByteBuffer.wrap(data);
        bb.order(ByteOrder.LITTLE_ENDIAN);
        return bb;
    }

    private void readHeaderAndFooter() throws IOException {

        // Read initial fields magic, version, and footer position
        byte[] header = readBytes(0, 16);
        if (header == null || header.length == 0) throw new IOException("File content is empty");

        ByteBuffer bp = wrap(header);
        this.magic = getString(bp);
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
        ByteBuffer bodyParser = wrap(body);

        this.genomeId = getString(bodyParser);  // This is not used, but could be useful for validation

        // The normalization vector index position and size.
        // This allows skipping reading expected values, which is expensive.
        if (this.version >= 9) {
            this.normVectorIndexPosition = bodyParser.getLong();
            this.normVectorIndexSize = (int) bodyParser.getLong();
        } else {
            // See if normVectorIndex is in table
            String nviStr = NVI.getNVI(this.path);
            if (nviStr != null) {
                try {
                    String[] parts = nviStr.split(",");
                    this.normVectorIndexPosition = Long.parseLong(parts[0]);
                    this.normVectorIndexSize = Integer.parseInt(parts[1]);
                } catch (NumberFormatException e) {
                    log.error("Error parsing NVI string: " + nviStr, e);
                }
            }
        }

        // attributes
        this.attributes = new HashMap<>();
        int nAttributes = bodyParser.getInt();
        while (nAttributes-- > 0) {
            String k = getString(bodyParser);
            String v = getString(bodyParser);
            this.attributes.put(k, v);
        }

        // chromosomes
        this.chromosomes = new ArrayList<>();
        this.chromosomeIndexMap = new CaseInsensitiveMap<>();
        int nChrs = bodyParser.getInt();
        for (int i = 0; i < nChrs; i++) {
            String name = getString(bodyParser);
            long size = this.version < 9 ? bodyParser.getInt() : bodyParser.getLong();
            Chromosome chr = new Chromosome(i, name, (int) size);
            this.chromosomes.add(chr);

            String canonicalName = genome == null ? name : genome.getCanonicalChrName(name);
            chrAliasTable.put(canonicalName, name);

            this.chromosomeIndexMap.put(name, i);

        }

        // bp resolutions
        int nBp = bodyParser.getInt();
        for (int i = 0; i < nBp; i++) {
            this.bpResolutions.add(bodyParser.getInt());
        }

        // frag resolutions --  For possible future use, this will never evaluate to true in IGV.
        boolean loadFragData = false;
        if (loadFragData) {
            int nFrag = bodyParser.getInt();
            for (int i = 0; i < nFrag; i++) {
                this.fragResolutions.add(bodyParser.getInt());
            }
        }
    }

    private void readFooter() throws IOException {
        int skip = this.version < 9 ? 8 : 12;
        byte[] data = readBytes(this.footerPosition, skip);
        if (data == null) return;
        ByteBuffer bp = wrap(data);
        long nBytes = this.version < 9 ? bp.getInt() : bp.getLong();
        int nEntries = bp.getInt();

        int miSize = nEntries * (100 + 64 + 32);
        byte[] miData = readBytes(this.footerPosition + skip, Math.min(miSize, (int) nBytes));
        ByteBuffer miParser = wrap(miData);

        this.masterIndex = new HashMap<>();
        while (nEntries-- > 0) {
            String key = getString(miParser);
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

    private Matrix getMatrix(int chrIdx1, int chrIdx2) throws IOException {
        String key = Matrix.getKey(chrIdx1, chrIdx2);
        if (matrixCache.containsKey(key)) return matrixCache.get(key);
        Matrix m = readMatrix(chrIdx1, chrIdx2);
        if (m != null) matrixCache.put(key, m);
        return m;
    }

    private Matrix readMatrix(int chrIdx1, int chrIdx2) throws IOException {

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

    public List<ContactRecord> getContactRecords(Region region1,
                                                 Region region2,
                                                 String units,
                                                 int binSize,
                                                 String normalization,
                                                 boolean allRecords) throws IOException {
        return this.getContactRecords(region1, region2, units, binSize, normalization, allRecords, 1);

    }

    public List<ContactRecord> getContactRecords(Region region1,
                                                 Region region2,
                                                 String units,
                                                 int binSize,
                                                 String normalization,
                                                 boolean allRecords,
                                                 int countsTreshold) throws IOException {

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


        boolean useNormalization = normalization != null && !"NONE".equals(normalization);
        NormalizationVector nv = useNormalization ? getNormalizationVector(normalization, region1.chr(), "BP", binSize) : null;
        if (nv == null) {
            useNormalization = false;
        }


        for (Block block : blocks) {
            if (block == null) continue;

            double[] normVector = null;

            // Compute binMin and binMax from block.records
            int binMin = Integer.MAX_VALUE;
            int binMax = Integer.MIN_VALUE;
            if (useNormalization) {
                for (ContactRecord rec : block.records) {
                    binMin = Math.min(binMin, Math.min(rec.bin1(), rec.bin2()));
                    binMax = Math.max(binMax, Math.max(rec.bin1(), rec.bin2()));
                }
                normVector = nv.getValues(binMin, binMax);
            }

            for (ContactRecord rec : block.records) {
                if ((allRecords || (rec.bin1() >= x1 && rec.bin1() < x2 && rec.bin2() >= y1 && rec.bin2() < y2)) && rec.counts() > countsTreshold) {
                    if (normVector == null) {
                        contactRecords.add(rec);
                    } else {
                        double nvnv = normVector[rec.bin1() - binMin] * normVector[rec.bin2() - binMin];
                        if (!Double.isNaN(nvnv)) {
                            float normCounts = (float) (rec.counts() / nvnv);
                            ContactRecord normRec = new ContactRecord(rec.bin1(), rec.bin2(), rec.counts(), normCounts);
                            contactRecords.add(normRec);
                        }
                    }
                }
            }
        }

        return contactRecords;
    }

    public int getWGResolution() {
        if (wgResolution == null) {
            try {
                Integer idx = chromosomeIndexMap.get("all");
                if (idx == null) return -1;
                Matrix matrix = getMatrix(idx, idx);
                if (matrix == null) return -1;
                List<MatrixZoomData> zdArray = matrix.getBpZoomData();
                if (zdArray.isEmpty()) return -1;
                wgResolution = zdArray.get(0).getZoom().binSize();
            } catch (IOException e) {
                log.error(e.getMessage());
                wgResolution = -1;
            }
        }
        return wgResolution;
    }

    public int getBinSize(String chr, double bpPerPixel) {

        if ("all".equalsIgnoreCase(chr)) {
            // Special case, the whole-genome psuedo-chromosome all has a single resolution
            return this.getWGResolution();
        }

        // choose resolution
        List<Integer> resolutions = this.getBpResolutions();
        int index = 0;
        for (int i = resolutions.size() - 1; i >= 0; i--) {
            if (resolutions.get(i) >= bpPerPixel) {
                index = i;
                break;
            }
        }
        int binSize = resolutions.get(index);
        return binSize;

    }

    private List<Block> getBlocks(Region region1, Region region2, String unit, int binSize) throws IOException {
        init();
        String chr1 = getFileChrName(region1.chr());
        String chr2 = getFileChrName(region2.chr());
        Integer idx1 = chromosomeIndexMap.get(chr1);
        Integer idx2 = chromosomeIndexMap.get(chr2);
        if (idx1 == null || idx2 == null) return Collections.emptyList();

        Matrix matrix = getMatrix(idx1, idx2);
        if (matrix == null) return Collections.emptyList();

        MatrixZoomData zd = matrix.getZoomData(binSize, unit);
        if (zd == null) {
            log.info("No data available for resolution: " + binSize + " for chromosome pair: " + chr1 + ", " + chr2);
            return Collections.emptyList();
        }

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

        // Fetch blocks that aren't in cache, with deduplication for concurrent requests
        for (Integer bn : toQuery) {
            String key = zd.getKey() + "_" + bn;
            Block block = getBlockWithDeduplication(key, bn, binSize, zd);
            blocks.add(block);
        }

        return blocks;
    }

    /**
     * Get a block with request deduplication. If another thread is already fetching the same block,
     * this method will wait for that operation to complete and return the same result instead of
     * initiating a redundant read.
     */
    private Block getBlockWithDeduplication(String key, int blockNumber, int binSize, MatrixZoomData zd) throws IOException {
        // Double-check cache (in case another thread just added it)
        if (blockCache.has(binSize, key)) {
            return blockCache.get(binSize, key);
        }

        // Create a new future for this request
        CompletableFuture<Block> newFuture = new CompletableFuture<>();

        // Try to be the first to register this request
        CompletableFuture<Block> existingFuture = pendingBlockRequests.putIfAbsent(key, newFuture);

        if (existingFuture != null) {
            // Another thread is already fetching this block - wait for its result
            try {
                return existingFuture.get();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new IOException("Interrupted while waiting for block data", e);
            } catch (ExecutionException e) {
                Throwable cause = e.getCause();
                if (cause instanceof IOException) {
                    throw (IOException) cause;
                }
                throw new IOException("Error reading block", cause);
            }
        }

        // This thread won the race - perform the actual read
        try {
            Block block = readBlock(blockNumber, zd);
            if (block != null) {
                blockCache.set(binSize, key, block);
            }
            newFuture.complete(block);
            return block;
        } catch (IOException e) {
            newFuture.completeExceptionally(e);
            throw e;
        } finally {
            // Clean up the pending request
            pendingBlockRequests.remove(key);
        }
    }

    private Block readBlock(int blockNumber, MatrixZoomData zd) throws IOException {
        StaticBlockIndex.BlockIndexEntry idx = zd.getBlockIndex().getBlockIndexEntry(blockNumber);
        if (idx == null) return null;
        byte[] data = readBytes(idx.filePosition, idx.size);
        if (data == null) return null;

        // decompress
        byte[] plain = (new CompressionUtils()).decompress(data);

        ByteBuffer parser = wrap(plain);
        int nRecords = parser.getInt();
        List<ContactRecord> records = new ArrayList<>();

        if (this.version < 7) {
            for (int i = 0; i < nRecords; i++) {
                int binX = parser.getInt();
                int binY = parser.getInt();
                float counts = parser.getFloat();
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
                        float counts = useFloatContact ? parser.getFloat() : parser.getShort();
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
                        float counts = parser.getFloat();
                        if (!Float.isNaN(counts)) records.add(new ContactRecord(bin1, bin2, counts));
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

    public boolean hasNormalizationVector(String type, String chr, String unit, int binSize) {

        int chrIdx = chromosomeIndexMap.getOrDefault(getFileChrName(chr), -1);

        String key = getNormalizationVectorKey(type, chrIdx, unit, binSize);
        if (normVectorCache.containsKey(key)) {
            return true;
        }

        try {
            Map<String, IndexEntry> nvi = getNormVectorIndex();
            return nvi == null ? false : nvi.containsKey(key);
        } catch (IOException e) {
            log.error("Error reading norm vector index", e);
            return false;
        }

    }

    /**
     * Look up and return the normalization vector for the given parameters.
     * @param type
     * @param chr
     * @param unit
     * @param binSize
     * @return
     * @throws IOException
     */
    public NormalizationVector getNormalizationVector(String type, String chr, String unit, int binSize) throws IOException {
        init();

        int chrIdx = chromosomeIndexMap.getOrDefault(getFileChrName(chr), -1);


        String key = getNormalizationVectorKey(type, chrIdx, unit, binSize);
        if (normVectorCache.containsKey(key)) {
            return normVectorCache.get(key);
        }

        Map<String, IndexEntry> nvi = getNormVectorIndex();
        if (nvi == null) return null;

        IndexEntry idx = nvi.get(key);
        if (idx == null) return null;

        byte[] header = readBytes(idx.start(), 8);
        if (header == null) return null;

        ByteBuffer parser = wrap(header);
        long nValues = this.version < 9 ? parser.getInt() : parser.getLong();
        DataType dataType = this.version < 9 ? DataType.DOUBLE : DataType.FLOAT;
        long filePos = this.version < 9 ? idx.start() + 4 : idx.start() + 8;
        NormalizationVector nv = new NormalizationVector(fileChannel, filePos, (int) nValues, dataType);
        normVectorCache.put(key, nv);
        return nv;
    }

    private Map<String, IndexEntry> getNormVectorIndex() throws IOException {
        if (this.version < 6) return null;
        if (this.normVectorIndex == null) {

            // If we know the norm vector index position and size, read it directly.
            if (this.normVectorIndexPosition > 0 && this.normVectorIndexSize > 0) {
                readNormVectorIndex(new IndexEntry(this.normVectorIndexPosition, this.normVectorIndexSize));
            } else {
                try {
                    // The norm vector index is located after the expected value vectors.  We need to skip over those first.
                    // This is in general slow, but works for all file versions.
                    readNormExpectedValuesAndNormVectorIndex();
                } catch (IOException e) {
                    // File may not contain norm vectors.  This rare but not neccessarily an error.
                    log.warn("Error reading norm vector index.  This could indicate normalization vectors are not present", e);
                }
            }
        }
        return normVectorIndex;
    }

    private void readNormVectorIndex(IndexEntry range) throws IOException {
        init();
        byte[] data = readBytes(range.start, range.size);
        ByteBuffer bp = wrap(data);
        this.normVectorIndex = new HashMap<>();
        int nEntries = bp.getInt();
        while (nEntries-- > 0) {
            parseNormVectorEntry(bp);
        }
    }

    private void readNormExpectedValuesAndNormVectorIndex() throws IOException {

        init();

        if (this.normExpectedValueVectorsPosition == null) {
            return;
        }
        long nviStart = skipExpectedValues(this.normExpectedValueVectorsPosition);

        byte[] data = readBytes(nviStart, 4);
        if (data == null || data.length == 0) {
            // This is possible if the are no norm verctors.  Not common.
            return;
        }
        ByteBuffer bp = wrap(data);
        int nEntries = bp.getInt();
        int sizeEstimate = nEntries * 30;
        data = readBytes(nviStart + 4, sizeEstimate);


        this.normVectorIndex = new HashMap<>();

        // process entries
        int byteCount = 4;
        ByteBuffer parser = wrap(data);
        while (nEntries-- > 0) {
            if (parser.remaining() < 100) {
                nEntries++;   // Reset counter as entry is not processed

                byteCount += parser.position();
                sizeEstimate = Math.max(1000, nEntries * 30);
                data = readBytes(nviStart + byteCount, sizeEstimate);
                parser = wrap(data);
            }
            parseNormVectorEntry(parser);

        }


        //   this.config.put("nvi", nviStart + "," + 4);
    }

    private long skipExpectedValues(long start) throws IOException {

        final int INT_SIZE = 4;
        final int DOUBLE_SIZE = 8;
        final int FLOAT_SIZE = 4;

        byte[] data = readBytes(start, INT_SIZE);
        if (data == null) {
            return start;
        }
        ByteBuffer buffer = wrap(data);
        int nEntries = buffer.getInt();

        if (nEntries == 0) {
            return start + INT_SIZE;
        }

        long currentPosition = start + INT_SIZE;
        for (int i = 0; i < nEntries; i++) {
            // Read the header of the expected value record to determine its size
            // The size is variable, so we read a reasonable chunk first
            byte[] chunkHeader = readBytes(currentPosition, 500);
            if (chunkHeader == null) {
                throw new IOException("Unexpected end of file while reading expected values.");
            }
            ByteBuffer chunkBuffer = wrap(chunkHeader);

            // Skip type and unit strings
            getString(chunkBuffer); // type
            getString(chunkBuffer); // unit
            chunkBuffer.getInt();   // binSize

            // Get # of values and compute size
            long nValues = (version < 9) ? chunkBuffer.getInt() : chunkBuffer.getLong();
            long valuesSize = nValues * ((version < 9) ? DOUBLE_SIZE : FLOAT_SIZE);

            // Position after the values array
            long posAfterValues = currentPosition + chunkBuffer.position() + valuesSize;

            // Read nChrScaleFactors
            byte[] scaleFactorsHeader = readBytes(posAfterValues, INT_SIZE);
            if (scaleFactorsHeader == null) {
                throw new IOException("Unexpected end of file while reading scale factors.");
            }
            ByteBuffer scaleFactorsBuffer = wrap(scaleFactorsHeader);
            int nChrScaleFactors = scaleFactorsBuffer.getInt();

            long scaleFactorsSize = (long) nChrScaleFactors * (INT_SIZE + ((version < 9) ? DOUBLE_SIZE : FLOAT_SIZE));

            // The start of the next entry is the position after the current one.
            currentPosition = posAfterValues + INT_SIZE + scaleFactorsSize;
        }

        return currentPosition;

    }

    private void parseNormVectorEntry(ByteBuffer parser) {
        String type = getString(parser);
        int chrIdx = parser.getInt();
        String unit = getString(parser);
        int binSize = parser.getInt();
        long filePosition = parser.getLong();
        long sizeInBytes = this.version < 9 ? parser.getInt() : parser.getLong();
        String key = getNormalizationVectorKey(type, chrIdx, unit, binSize);
        if (!this.normalizationTypes.contains(type)) {
            this.normalizationTypes.add(type);
        }
        this.normVectorIndex.put(key, new IndexEntry(filePosition, (int) sizeInBytes));
    }

    private String getFileChrName(String chrName) {
        return chrAliasTable.getOrDefault(chrName, chrName);
    }

    public List<String> getNormalizationTypes() {
        try {
            getNormVectorIndex();   // Populate normalization types as a sid
        } catch (IOException e) {
            log.error("Error reading norm vector index", e);
        }
        return normalizationTypes;
    }

    private static String getNormalizationVectorKey(String type, int chrIdx, String unit, int resolution) {
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

    /**
     * Read a null-terminated string from the buffer
     * @param buffer
     * @return
     */
    private static String getString(ByteBuffer buffer) {
        ByteArrayOutputStream bis = new ByteArrayOutputStream(1000);
        int b;
        while ((b = buffer.get()) != 0) {
            bis.write((byte) b);
        }
        return new String(bis.toByteArray());
    }

    public void setNVIString(String nviString) {
        String[] parts = nviString.split(",");
        if (parts.length < 2) {
            log.error("Invalid NVI string: " + nviString);
            this.normVectorIndexPosition = -1;
            this.normVectorIndexSize = -1;
            return;
        }
        try {
            this.normVectorIndexPosition = Long.parseLong(parts[0]);
            this.normVectorIndexSize = Integer.parseInt(parts[1]);
        } catch (NumberFormatException e) {
            log.error("Error parsing NVI string: " + nviString, e);
            this.normVectorIndexPosition = -1;
            this.normVectorIndexSize = -1;
        }
    }

    record IndexEntry(long start, int size) {
    }

    private static class Block {
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

    private static class BlockCache {
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

}