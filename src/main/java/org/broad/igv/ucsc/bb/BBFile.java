package org.broad.igv.ucsc.bb;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.NamedFeature;
import org.broad.igv.data.BasicScore;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureType;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.ChromAlias;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ucsc.BPTree;
import org.broad.igv.ucsc.Trix;
import org.broad.igv.ucsc.bb.codecs.*;
import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;
import org.broad.igv.ucsc.twobit.UnsignedByteBufferImpl;
import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.net.URL;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.stream.Collectors;

/* bigWig/bigBed file structure:
 *     fixedWidthHeader
 *         magic# 		4 bytes
 *         version              2 bytes
 *	   zoomLevels		2 bytes
 *         chromosomeTreeOffset	8 bytes
 *         fullDataOffset	8 bytes
 *	   fullIndexOffset	8 bytes
 *         fieldCount           2 bytes (for bigWig 0)
 *         definedFieldCount    2 bytes (for bigWig 0)
 *         autoSqlOffset        8 bytes (for bigWig 0) (0 if no autoSql information)
 *         totalSummaryOffset   8 bytes (0 in earlier versions of file lacking totalSummary)
 *         uncompressBufSize    4 bytes (Size of uncompression buffer.  0 if uncompressed.)
 *         extensionOffset      8 bytes (Offset to header extension 0 if no such extension)
 *     zoomHeaders		there are zoomLevels number of these
 *         reductionLevel	4 bytes
 *	   reserved		4 bytes
 *	   dataOffset		8 bytes
 *         indexOffset          8 bytes
 *     autoSql string (zero terminated - only present if autoSqlOffset non-zero)
 *     totalSummary - summary of all data in file - only present if totalSummaryOffset non-zero
 *         basesCovered        8 bytes
 *         minVal              8 bytes float (for bigBed minimum depth of coverage)
 *         maxVal              8 bytes float (for bigBed maximum depth of coverage)
 *         sumData             8 bytes float (for bigBed sum of coverage)
 *         sumSquared          8 bytes float (for bigBed sum of coverage squared)
 *     extendedHeader
 *         extensionSize       2 size of extended header in bytes - currently 64
 *         extraIndexCount     2 number of extra fields we will be indexing
 *         extraIndexListOffset 8 Offset to list of non-chrom/start/end indexes
 *         reserved            48 All zeroes for now
 *     extraIndexList - one of these for each extraIndex
 *         type                2 Type of index.  Always 0 for bPlusTree now
 *         fieldCount          2 Number of fields used in this index.  Always 1 for now
 *         indexOffset         8 offset for this index in file
 *         reserved            4 All zeroes for now
 *         fieldList - one of these for each field being used in _this_ index
 *            fieldId          2 index of field within record
 *            reserved         2 All zeroes for now
 *     chromosome b+ tree       bPlusTree index
 *     full data
 *         sectionCount		8 bytes (item count for bigBeds)
 *         section data		section count sections, of three types (bed data for bigBeds)
 *     full index               cirTree index
 *     zoom info             one of these for each zoom level
 *         zoom data
 *             zoomCount	4 bytes
 *             zoom data	there are zoomCount of these items
 *                 chromId	4 bytes
 *	           chromStart	4 bytes
 *                 chromEnd     4 bytes
 *                 validCount	4 bytes
 *                 minVal       4 bytes float
 *                 maxVal       4 bytes float
 *                 sumData      4 bytes float
 *                 sumSquares   4 bytes float
 *         zoom index        	cirTree index
 *     extraIndexes [optional]  bPlusTreeIndex for each extra field that is indexed
 *     magic# 		4 bytes - same as magic number at start of header
 */
public class BBFile {


    public enum Type {BIGWIG, BIGBED}

    static Logger log = LogManager.getLogger(BBFile.class);


    static public final int BBFILE_HEADER_SIZE = 64;
    static public final long BIGWIG_MAGIC = 2291137574l; // BigWig Magic
    static public final long BIGBED_MAGIC = 2273964779l; // BigBed Magic
    static public final int BBFILE_EXTENDED_HEADER_HEADER_SIZE = 64;
    private Trix trix;
    private String autosql;
    private ChromTree chromTree;
    private double featureDensity;
    private Map<String, String> chrAliasTable;
    private BBTotalSummary totalSummary;
    private BPTree[] _searchTrees;
    private Map<Long, RPTree> rTreeCache;
    private BBCodec bedCodec;
    private String path;
    private Type type;
    private FeatureType featureType;
    private BBHeader header = null;
    private BBZoomHeader[] zoomHeaders;
    private ByteOrder byteOrder;
    private Genome genome;
    private byte[] _preloadBytes;

    public BBFile(String path, Genome genome) throws IOException {
        this.path = path;
        this.genome = genome;
        this.chrAliasTable = new HashMap<>();
        this.rTreeCache = new HashMap<>();
        init();
    }

    public void preload() throws IOException {
        if (FileUtils.isRemote(path)) {
            this._preloadBytes = HttpUtils.getInstance().getContentsAsBytes(new URL(this.path), null);
        }
    }

    public BBFile(String path, Genome genome, String trixPath) throws IOException {
        this(path, genome);
        this.trix = new Trix(trixPath + "x", trixPath);
    }

    void init() throws IOException {
        this.header = readHeader();
    }

    public int getChromosomeCount() {
        return (int) chromTree.getItemCount();
    }

    public Type getType() {
        return type;
    }

    public FeatureType getFeatureType() {
        if (featureType == null) {
            if (isBigWigFile()) {
                featureType = FeatureType.WIG;
            } else if (isBigBedFile()) {
                String as = getAutoSQL();
                if (as != null) {
                    BBUtils.ASTable astable = BBUtils.parseAutosql(autosql);
                    if (astable != null && "bigRmskBed".equals(astable.name)) {
                        featureType = FeatureType.RMSK;
                    } else if (astable != null && astable.name.toLowerCase().contains("interact")) {
                        featureType = FeatureType.INTERACT;
                    }
                }
                if (featureType == null) {
                    featureType = FeatureType.BED;
                }
            }
        }
        return featureType;
    }

    public BBHeader getHeader() {
        return header;
    }

    public Genome getGenome() {
        return genome;
    }

    public boolean isBigWigFile() {
        return type == Type.BIGWIG;
    }

    public boolean isBigBedFile() {
        return type == Type.BIGBED;
    }

    public double getFeatureDensity() {
        return featureDensity;
    }

    public String getAutoSQL() {
        return autosql;
    }

    BBHeader readHeader() throws IOException {

        // The common header
        ByteOrder order = ByteOrder.LITTLE_ENDIAN;
        UnsignedByteBuffer buffer = loadBinaryBuffer(this.path, order, 0, BBFILE_HEADER_SIZE);
        long magic = buffer.getUInt();
        if (magic == BIGWIG_MAGIC) {
            this.type = Type.BIGWIG;
        } else if (magic == BIGBED_MAGIC) {
            this.type = Type.BIGBED;
        } else {
            //Try big endian order
            order = ByteOrder.BIG_ENDIAN;
            buffer = loadBinaryBuffer(this.path, order, 0, BBFILE_HEADER_SIZE);
            magic = buffer.getUInt();
            if (magic == BIGWIG_MAGIC) {
                this.type = Type.BIGWIG;
            } else if (magic == BIGBED_MAGIC) {
                this.type = Type.BIGBED;
            } else {
                throw new RuntimeException("Bad magic number " + magic);
            }
        }
        this.byteOrder = order;

        BBHeader header = new BBHeader();
        header.version = buffer.getUShort();
        header.nZoomLevels = buffer.getUShort();
        header.chromTreeOffset = buffer.getLong();
        header.fullDataOffset = buffer.getLong();
        header.fullIndexOffset = buffer.getLong();
        header.fieldCount = buffer.getUShort();
        header.definedFieldCount = buffer.getUShort();
        header.autoSqlOffset = buffer.getLong();
        header.totalSummaryOffset = buffer.getLong();
        header.uncompressBuffSize = buffer.getInt();
        header.extensionOffset = buffer.getLong();
        this.header = header;

        // Zoom headers, autosql, and total summary if present
        int size = (int) (header.totalSummaryOffset > 0 ?
                header.totalSummaryOffset - BBFILE_HEADER_SIZE + 40 :
                Math.min(header.fullDataOffset, header.chromTreeOffset) - BBFILE_HEADER_SIZE);

        buffer = loadBinaryBuffer(this.path, order, BBFILE_HEADER_SIZE, size);

        // Zoom headers -- immediately follows the common header
        this.zoomHeaders = new BBZoomHeader[header.nZoomLevels];
        for (int i = 0; i < header.nZoomLevels; ++i) {
            BBZoomHeader zlh = new BBZoomHeader();
            zlh.reductionLevel = buffer.getUInt();
            zlh.reserved = buffer.getInt();
            zlh.dataOffset = buffer.getLong();
            zlh.indexOffset = buffer.getLong();
            this.zoomHeaders[i] = zlh;
        }
        // Sort in order of decreasing reduction level (increasing resolution
        Arrays.sort(zoomHeaders, (o1, o2) -> (int) (o2.reductionLevel - o1.reductionLevel));

        // Autosql -- spec implies this follows the zoom headers
        final int startOffset = BBFILE_HEADER_SIZE;
        if (header.autoSqlOffset > 0) {
            buffer.position((int) (header.autoSqlOffset - startOffset));
            this.autosql = buffer.getString();
        }

        // Total summary -- present in versions >= 2.  Follows the zoom headers and autosql
        if (header.version > 1 && header.totalSummaryOffset > 0) {
            buffer.position((int) (header.totalSummaryOffset - startOffset));
            this.totalSummary = BBTotalSummary.parseSummary(buffer);
        }

        // Chromosome tree
        this.chromTree = new ChromTree(this.path, this.header.chromTreeOffset);

        if (type == Type.BIGBED) {
            //Total data count -- for bigbed this is the number of features, for bigwig it is number of sections
            buffer = loadBinaryBuffer(this.path, order, header.fullDataOffset, 4);
            header.dataCount = buffer.getInt();
            this.featureDensity = ((double) header.dataCount) / this.chromTree.estimateGenomeSize();

            bedCodec = BBCodecFactory.getCodec(autosql, header.definedFieldCount);
        }

        //extension
        if (header.extensionOffset > 0) {
            this.loadExtendedHeader(header.extensionOffset);
        }


        return header;

    }

    void loadExtendedHeader(long offset) throws IOException {

        UnsignedByteBuffer binaryParser = loadBinaryBuffer(this.path, byteOrder, offset, BBFILE_EXTENDED_HEADER_HEADER_SIZE);

        int extensionSize = binaryParser.getUShort();
        int extraIndexCount = binaryParser.getUShort();
        long extraIndexListOffset = binaryParser.getLong();
        if (extraIndexCount == 0) return;

        int sz = extraIndexCount * (2 + 2 + 8 + 4 + 10 * (2 + 2));
        binaryParser = loadBinaryBuffer(this.path, byteOrder, extraIndexListOffset, sz);

        // const type = []
        // const fieldCount = []
        // const reserved = []
        // const indexOffset = []
        long[] indexOffset = new long[extraIndexCount];
        for (int i = 0; i < extraIndexCount; i++) {

            //type.push(binaryParser.getUShort())
            int type = binaryParser.getUShort();

            int fc = binaryParser.getUShort();
            //fieldCount.push(fc)


            indexOffset[i] = binaryParser.getLong();
            //reserved.push(binaryParser.getInt())
            binaryParser.getInt();

            for (int j = 0; j < fc; j++) {
                int fieldId = binaryParser.getUShort();
                //const field = this.autoSql.fields[fieldId]
                //console.log(field)

                //reserved.push(binaryParser.getUShort())
                binaryParser.getUShort();
            }
        }
        this.header.extraIndexCount = extraIndexCount;
        this.header.extraIndexOffsets = indexOffset;
    }

    List<byte[]> getLeafChunks(int chrIdx1, int bpStart, int chrIdx2, int bpEnd, long treeOffset) throws IOException {

        if (this.header == null) {
            this.header = this.readHeader();
        }

        // Load the R Tree and fine leaf items
        RPTree rpTree = rTreeCache.get(treeOffset);
        if (rpTree == null) {
            rpTree = RPTree.loadTree(this.path, treeOffset);
            rTreeCache.put(treeOffset, rpTree);
        }

        List<byte[]> leafChunks = new ArrayList<>();
        List<RPTree.Item> leafItems = rpTree.findLeafItemsOverlapping(chrIdx1, bpStart, chrIdx2, bpEnd);
        if (leafItems != null && leafItems.size() > 0) {

            // Consolidate leaf items and get all data at once
            long start = Long.MAX_VALUE;
            long end = 0;
            for (RPTree.Item item : leafItems) {
                start = Math.min(start, item.dataOffset);
                end = Math.max(end, item.dataOffset + item.dataSize);
            }
            int size = (int) (end - start);
            int uncompressBufSize = this.header.uncompressBuffSize;

            byte[] buffer = getBytes(start, size);

            for (RPTree.Item item : leafItems) {
                int offset = (int) (item.dataOffset - start);
                int end_ = (int) (offset + item.dataSize);
                byte[] itemBuffer = leafItems.size() == 1 ? buffer : Arrays.copyOfRange(buffer, offset, end_);
                leafChunks.add(itemBuffer);
            }

        }

        return leafChunks;

    }


    /**
     * Return the ID for the given chromosome name.  If there is no direct match, search for a chromosome alias.
     *
     * @param chr -- the canonical chromsome name for the current reference
     */
    Integer getIdForChr(String chr) throws IOException {

        if (chrAliasTable.get(chr) != null) {
            chr = chrAliasTable.get(chr);
        }

        Integer chrIdx = chromTree.getIdForName(chr);

//        // Try alias
        if (chrIdx == null && genome != null) {
            String alias = null;
            ChromAlias aliasRecord = genome.getAliasRecord(chr);
            if (aliasRecord != null) {
                for (String v : aliasRecord.values()) {
                    chrIdx = chromTree.getIdForName(v);
                    if (chrIdx != null) {
                        alias = v;
                        break;
                    }
                }
            }
            this.chrAliasTable.put(chr, alias);  // alias may be undefined => no alias exists. Setting prevents repeated attempts
        }

        return chrIdx;
    }

    /**
     * Return the zoom header that most closely matches the given resolution.  Resolution is in BP / Pixel.
     *
     * @param bpPerPixel -- the resolution in bp per pixel.
     * @return A zoom header, or null if no appropriate zoom data is available for the resolution.
     */

    BBZoomHeader zoomLevelForScale(double bpPerPixel) {
        return zoomLevelForScale(bpPerPixel, 2);
    }

    BBZoomHeader zoomLevelForScale(double bpPerPixel, int tolerance) {
        if(this.zoomHeaders.length == 0) {
            return null;
        }
        for (BBZoomHeader zl : this.zoomHeaders) {
            if (zl.reductionLevel < bpPerPixel) {
                return zl;
            }
        }
        // For the lowest resolution, consider a match if within a factor "tolerance" of the requested resolution
        BBZoomHeader lastLevel = this.zoomHeaders[this.zoomHeaders.length - 1];
        return lastLevel.reductionLevel / tolerance < bpPerPixel ? lastLevel : null;
    }

    public boolean isSearchable() {
        return header.extraIndexCount > 0;
    }

    /**
     * Search the extended BP tree for the search term, and return any matching features.  This only works
     * for BB sources with an "extended" BP tree for searching.
     * <p>
     * Currently we don't support multiple hits for a search.  Keep the largest.
     *
     * @param term
     * @returns {Promise<void>}
     */


    public IGVFeature search(String term) throws IOException {

        if (this.header == null) {
            this.readHeader();
        }
        if (this.header.extraIndexCount == 0) {
            return null;
        }

        if (this.trix != null) {
            String termLower = term.toLowerCase();
            Map<String, String[]> results = trix.search(termLower);
            if (results != null && results.containsKey(termLower)) {
                String[] exactMatches = results.get(termLower);
                if (exactMatches.length > 0) term = exactMatches[0];
            }
        }


        long[] region = this.searchForRegions(term);  // Either 1 or no (undefined) reginos returned for now
        if (region != null) {
            long start = region[0];
            int size = (int) region[1];
            byte[] buffer = this.getBytes(start, size);
            List<IGVFeature> features = decodeFeatures(null, buffer, -1, -1, -1);


            // Filter features to those matching term
            final String searchTerm = term;
            return features.stream().filter(f -> {
                return f.getName().equalsIgnoreCase(searchTerm) ||
                        f.getAttributes().values().stream().anyMatch(v -> v.equalsIgnoreCase(searchTerm));
            }).collect(Collectors.toList());

        }
        return null;
    }


    List<IGVFeature> decodeFeatures(String chr, byte[] buffer, int chrIdx, int start, int end) {

        List<IGVFeature> features = new ArrayList<>();

        byte[] uncompressed;
        if (this.header.uncompressBuffSize > 0) {
            uncompressed = (new CompressionUtils()).decompress(buffer, this.header.uncompressBuffSize);
        } else {
            uncompressed = buffer;    // use uncompressed read buffer directly
        }


        UnsignedByteBufferImpl bb = UnsignedByteBufferImpl.wrap(uncompressed, byteOrder);
        while (bb.remaining() > 0) {

            int chromId = bb.getInt();
            int chromStart = bb.getInt();
            int chromEnd = bb.getInt();
            String restOfFields = bb.getString();

            if (chrIdx > 0) {
                if (chromId < chrIdx || (chromId == chrIdx && chromEnd < start)) continue;
                else if (chromId > chrIdx || (chromId == chrIdx && chromStart >= end)) break;
            }


            // When decoding features from a search the chromosome name is not known until the feature is decoded.
            // Try to get the chromosome name from the id.
            if (chr == null) {
                chr = this.chromTree.getNameForId(chromId);
            }
            if (chr != null) {
                final BedData bedData = new BedData(chr, chromStart, chromEnd, restOfFields);
                final IGVFeature feature = bedCodec.decode(bedData);
                features.add(feature);
            }

        }
        return features;
    }

    List<LocusScore> decodeZoomData(String chr, byte[] buffer, int chrIdx, int start, int end, WindowFunction windowFunction, List<LocusScore> features) {

        byte[] uncompressed;
        if (header.uncompressBuffSize > 0) {
            uncompressed = (new CompressionUtils()).decompress(buffer, header.uncompressBuffSize);
        } else {
            uncompressed = buffer;    // use uncompressed read buffer directly
        }

        UnsignedByteBufferImpl bb = UnsignedByteBufferImpl.wrap(uncompressed, byteOrder);
        while (bb.remaining() > 0) {

            int chromId = bb.getInt();
            int chromStart = bb.getInt();
            int chromEnd = bb.getInt();
            int validCount = bb.getInt();
            float minVal = bb.getFloat();
            float maxVal = bb.getFloat();
            float sumData = bb.getFloat();
            float sumSquares = bb.getFloat();

            if (chrIdx > 0) {
                if (chromId < chrIdx || (chromId == chrIdx && chromEnd < start)) continue;
                else if (chromId > chrIdx || (chromId == chrIdx && chromStart >= end)) break;
            }

            if (chr == null) {
                chr = this.chromTree.getNameForId(chromId);
            }

            if (chr != null) {
                float value;
                switch (windowFunction) {
                    case min:
                        value = minVal;
                        break;
                    case max:
                        value = maxVal;
                        break;
                    case mean:
                        value = sumData / validCount;
                        break;
                    default:
                        throw new RuntimeException("Unsupported window function: " + windowFunction);

                }

                final LocusScore feature = new WigDatum(chr, chromStart, chromEnd, value);
                features.add(feature);
            }
        }
        return features;
    }

    List<LocusScore> decodeWigData(byte[] buffer, int chrIdx, int start, int end, List<LocusScore> features) {

        byte[] uncompressed;
        if (header.uncompressBuffSize > 0) {
            uncompressed = (new CompressionUtils()).decompress(buffer, header.uncompressBuffSize);
        } else {
            uncompressed = buffer;    // use uncompressed read buffer directly
        }

        UnsignedByteBuffer binaryParser = UnsignedByteBufferImpl.wrap(uncompressed, byteOrder);
        int chromId = binaryParser.getInt();
        int blockStart = binaryParser.getInt();
        int chromStart = blockStart;
        int chromEnd = binaryParser.getInt();
        int itemStep = binaryParser.getInt();
        int itemSpan = binaryParser.getInt();
        byte type = binaryParser.get();
        byte reserved = binaryParser.get();
        int itemCount = binaryParser.getUShort();

        if (chromId >= chrIdx && chromId <= chrIdx) {

            int idx = 0;
            while (itemCount-- > 0) {
                float value = Float.POSITIVE_INFINITY;
                switch (type) {
                    case 1:
                        chromStart = binaryParser.getInt();
                        chromEnd = binaryParser.getInt();
                        value = binaryParser.getFloat();
                        break;
                    case 2:
                        chromStart = binaryParser.getInt();
                        value = binaryParser.getFloat();
                        chromEnd = chromStart + itemSpan;
                        break;
                    case 3:  // Fixed step
                        value = binaryParser.getFloat();
                        chromStart = blockStart + idx * itemStep;
                        chromEnd = chromStart + itemSpan;
                        ++idx;
                        break;
                }

                if (chromId < chrIdx || (chromId == chrIdx && chromEnd < start)) continue;
                else if (chromId > chrIdx || (chromId == chrIdx && chromStart >= end)) break;

                if (Float.isFinite(value)) {
                    //  const chr = chrDict[chromId]
                    features.add(new BasicScore(chromStart, chromEnd, value));
                }
            }
        }

        return features;
    }


    private long[] searchForRegions(String term) throws IOException {

        BPTree[] searchTrees = this.getSearchTrees();
        if (searchTrees != null) {

            // Use a trix index if we have one to map entered term to indexed value in bb file
            if (this.trix != null) {
                String termLower = term.toLowerCase();
                Map<String, String[]> trixResults = this.trix.search(termLower);
                if (trixResults != null && trixResults.containsKey(termLower)) {   // <= exact matches only for now
                    term = trixResults.get(termLower)[0];
                }
            }

            // For now take the first match, we don't support multiple results
            for (BPTree bpTree : searchTrees) {
                if (bpTree != null) {
                    long[] result = bpTree.searchLongLong(term);
                    if (result != null) {
                        return result;
                    }
                }
            }
        }
        return null;
    }

    BPTree[] getSearchTrees() throws IOException {
        if (this._searchTrees == null &&
                this.header.extraIndexOffsets != null &&
                this.header.extraIndexOffsets.length > 0) {
            this._searchTrees = new BPTree[this.header.extraIndexOffsets.length];
            int idx = 0;
            for (long offset : this.header.extraIndexOffsets) {
                BPTree bpTree = BPTree.loadBPTree(this.path, offset);
                this._searchTrees[idx++] = bpTree;
            }
        }
        return this._searchTrees;
    }

    UnsignedByteBuffer loadBinaryBuffer(String path, ByteOrder order, long start, int size) throws IOException {
        byte[] bytes = getBytes(start, size);
        return UnsignedByteBufferImpl.wrap(bytes, byteOrder);
    }

    private byte[] getBytes(long start, int size) throws IOException {

        if (_preloadBytes != null) {
            return Arrays.copyOfRange(_preloadBytes, (int) start, (int) start + size);
        } else {
            SeekableStream is = null;
            try {
                is = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
                byte[] bytes = new byte[size];
                is.seek(start);
                is.readFully(bytes);
                return bytes;
            } finally {
                if (is != null) {
                    try {
                        is.close();
                    } catch (IOException e) {
                        log.error("Error closing stream for " + path, e);
                    }
                }
            }
        }
    }

    /**
     * Method computes bigwig "totalSummary" statistics for comparison with values in file header.  The purpose is to
     * validate interpretation of file header values.   "sumOfSquares" is literally the sum of squares of the data,
     * not sum of the differences wrt the mean.  The bigwig specification is not clear on this point.   It states that
     * "the standard deviation can be easily computed", but its not clear how.
     *
     * From Jim Kent https://github.com/ucscGenomeBrowser/kent/blob/7d0bd2b7d089f94c9a8260c417246d4dbfce2523/src/lib/bbiWrite.c#L515
     * void bbiAddRangeToSummary(bits32 chromId, bits32 chromSize, bits32 start, bits32 end,
     * 	                         double val, int reduction, struct bbiSummary **pOutList)
     *       int size = end - start;
     *       double sum = size * val;
     *       double sumSquares = sum * val;
     *
     *
     * @throws IOException
     */
    public void computeStats() throws IOException {

        long rTreeOffset = getHeader().fullIndexOffset;
        double sum = 0;
        long basesCovered = 0;
        double sumSquares = 0;
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        for (int idx = 0; idx < getChromosomeCount(); idx++) {
            List<byte[]> chunks = this.getLeafChunks(idx, 0, idx, Integer.MAX_VALUE, rTreeOffset);
            for (byte[] c : chunks) {
                List<LocusScore> features = new ArrayList<>();
                decodeWigData(c, idx, 0, Integer.MAX_VALUE, features);
                for (LocusScore score : features) {

                    min = Math.min(min, score.getScore());
                    max = Math.max(max, score.getScore());
                    final float val = score.getScore() * (score.getEnd() - score.getStart());
                    sum += val;
                    sumSquares += (score.getScore() * val);
                    basesCovered += score.getEnd() - score.getStart();
                }
            }
        }
        System.out.println(totalSummary.printString());
        double mean = sum / basesCovered;
        System.out.println("basesCovered: " + basesCovered);
        System.out.println("min: " + min);
        System.out.println("max: " + max);
        System.out.println("sum: " + sum);
        System.out.println("sumSquares: " + sumSquares);
        System.out.println("mean: " + mean);
    }

}
