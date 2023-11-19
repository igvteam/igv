package org.broad.igv.ucsc.bb;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ucsc.UnsignedByteBuffer;
import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

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

    static public final int BBFILE_HEADER_SIZE = 64;
    static public final long BIGWIG_MAGIC = 2291137574l; // BigWig Magic Low to High
    static public final long BIGWIG_MAGIC_HTL = 654085990l; // BigWig Magic High to Low
    static public final long BIGBED_MAGIC = 2273964779l; // BigBed Magic Low to High
    static public final long BIGBED_MAGIC_HTL = 3958540679l; // BigBed Magic High to Low
    static public final int BBFILE_EXTENDED_HEADER_HEADER_SIZE = 64;
    String autosql;
    private BBTotalSummary totalSummary;
    private ChromTree chromTree;
    private HashSet<String> chrNames;

    enum Type {BIGWIG, BIGBED}

    String path;
    Type type;
    BBHeader header = null;

    BBZoomHeader[] zoomHeaders;
    ByteOrder byteOrder;
    Genome genome;

    public static BBFile openFile(String path) throws IOException {
        BBFile bbf = new BBFile(path);
        bbf.init();
        return bbf;
    }

    private BBFile(String path) throws IOException {
        this.path = path;
    }

    void init() throws IOException {
        this.header = readHeader();
    }

    BBHeader readHeader() throws IOException {

        ByteOrder order = ByteOrder.LITTLE_ENDIAN;
        UnsignedByteBuffer buffer = UnsignedByteBuffer.loadBinaryBuffer(this.path, order, 0, BBFILE_HEADER_SIZE);
        long magic = buffer.getUInt();
        if (magic == BIGWIG_MAGIC) {
            this.type = Type.BIGWIG;
        } else if (magic == BIGBED_MAGIC) {
            this.type = Type.BIGBED;
        } else {
            //Try big endian order
            order = ByteOrder.BIG_ENDIAN;
            buffer = UnsignedByteBuffer.loadBinaryBuffer(this.path, order, 0, BBFILE_HEADER_SIZE);
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

        // Read rest of fields up to full data offset
        buffer = UnsignedByteBuffer.loadBinaryBuffer(this.path, order, BBFILE_HEADER_SIZE, (int) (header.fullDataOffset - BBFILE_HEADER_SIZE));

        // Zoom headers
        this.zoomHeaders = new BBZoomHeader[header.nZoomLevels];
        for (int i = 1; i <= header.nZoomLevels; i++) {
            int zoomNumber = header.nZoomLevels - i;
            BBZoomHeader zlh = new BBZoomHeader();
            zlh.reductionLevel = buffer.getInt();
            zlh.reserved = buffer.getInt();
            zlh.dataOffset = buffer.getLong();
            zlh.indexOffset = buffer.getLong();
            this.zoomHeaders[zoomNumber] = zlh;
        }

        // Autosql
        final int startOffset = BBFILE_HEADER_SIZE;
        if (header.autoSqlOffset > 0) {
            buffer.position((int) (header.autoSqlOffset - startOffset));
            this.autosql = buffer.getString();
        }

        // Total summary -- present in versions >= 2
        if (header.totalSummaryOffset > 0) {
            buffer.position((int) (header.totalSummaryOffset - startOffset));
            this.totalSummary = BBTotalSummary.parseSummary(buffer);
        }

        // Chromosome tree
        buffer.position((int) (header.chromTreeOffset - startOffset));
        this.chromTree = ChromTree.parseTree(buffer, startOffset, this.genome);
        this.chrNames = new HashSet<>(Arrays.asList(this.chromTree.idToName));

        //Finally total data count
        //buffer.position((int)(header.fullDataOffset - startOffset));
        //header.dataCount = buffer.getInt();

        this.header = header;

        //extension
        if (header.extensionOffset > 0) {
            this.loadExtendedHeader(header.extensionOffset);
        }

        //this.setDefaultVisibilityWindow(header);

        return header;

    }

    void loadExtendedHeader(long offset) throws IOException {

        UnsignedByteBuffer binaryParser = UnsignedByteBuffer.loadBinaryBuffer(this.path, byteOrder, offset, BBFILE_EXTENDED_HEADER_HEADER_SIZE);

        int extensionSize = binaryParser.getUShort();
        int extraIndexCount = binaryParser.getUShort();
        long extraIndexListOffset = binaryParser.getLong();
        if (extraIndexCount == 0) return;

        int sz = extraIndexCount * (2 + 2 + 8 + 4 + 10 * (2 + 2));
        binaryParser = UnsignedByteBuffer.loadBinaryBuffer(this.path, byteOrder, extraIndexListOffset, sz);

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

    List<byte[]> getLeafChunks(String chr1, int bpStart, String chr2, int bpEnd, double bpPerPixel) throws IOException {

        if (this.header == null) {
            this.header = this.readHeader();
        }

        int chrIdx1 = this.getIdForChr(chr1);
        int chrIdx2 = this.getIdForChr(chr2);

//        if (chrIdx1 === undefined || chrIdx2 === undefined) {
//            return []
//        }
//
        long treeOffset;
        if (this.type == Type.BIGWIG) {
            // Select a biwig "zoom level" appropriate for the current resolution.
            BBZoomHeader zoomLevelHeader = zoomLevelForScale(bpPerPixel);
            if (zoomLevelHeader != null) {
                treeOffset = zoomLevelHeader.indexOffset;
            } else {
                treeOffset = this.header.fullIndexOffset;
            }
        } else {
            // bigbed, zoom data is not currently used in igv for bed type features
            treeOffset = this.header.fullIndexOffset;
        }


        // Load the R Tree and fine leaf items
        RPTree rpTree = RPTree.loadTree(this.path, treeOffset);
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

            try (SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(this.path)) {
                byte[] buffer = new byte[size];
                is.seek(start);
                is.readFully(buffer);

                for (RPTree.Item item : leafItems) {
                    int offset = (int) (item.dataOffset - start);
                    int end_ = (int) (offset + item.dataSize);
                    byte[] itemBuffer = leafItems.size() == 1 ? buffer : Arrays.copyOfRange(buffer, offset, end_);
                    if (uncompressBufSize > 0) {
                        byte[] uncompressed = (new CompressionUtils()).decompress(itemBuffer, uncompressBufSize);
                        leafChunks.add(uncompressed);
                    } else {
                        leafChunks.add(itemBuffer);    // use uncompressed read buffer directly
                    }
                }
            }
        }

        return leafChunks;


//
//            // Parse data and return features
//            const features = []
//            for (let item of leafItems) {
//                const uint8Array = new Uint8Array(arrayBuffer, item.dataOffset - start, item.dataSize)
//                let plain
//                const isCompressed = this.header.uncompressBuffSize > 0
//                if (isCompressed) {
//                    plain = BGZip.inflate(uint8Array)
//                } else {
//                    plain = uint8Array
//                }
//                decodeFunction.call(this, new DataView(plain.buffer), chrIdx1, bpStart, chrIdx2, bpEnd, features, this.chromTree.idToName, windowFunction)
//            }
//
//            features.sort(function (a, b) {
//                return a.start - b.start
//            })
//
//            return features
//        }
    }


    /**
     * Return the ID for the given chromosome name.  If there is no direct match, search for a chromosome alias.
     *
     * @param chr
     * @returns {Promise<*>}
     */
    int getIdForChr(String chr) {

//        if (this.chrAliasTable.has(chr)) {
//            chr = this.chrAliasTable.get(chr)
//            if (chr === undefined) {
//                return undefined
//            }
//        }

        int chrIdx = this.chromTree.nameToId.get(chr);

//        // Try alias
//        if (chrIdx === undefined) {
//            const aliasRecord = await this.genome.getAliasRecord(chr)
//            let alias
//            if (aliasRecord) {
//                const aliases = Object.keys(aliasRecord)
//                        .filter(k => k !== "start" && k !== "end")
//                    .map(k => aliasRecord[k])
//                    .filter(a => this.chromTree.nameToId.has(a))
//                if (aliases.length > 0) {
//                    alias = aliases[0]
//                    chrIdx = this.chromTree.nameToId.get(aliases[0])
//                }
//            }
//            this.chrAliasTable.set(chr, alias)  // alias may be undefined => no alias exists. Setting prevents repeated attempts
//        }
        return chrIdx;
    }

    String getChrForId(int chrIdx) {
        return chromTree.idToName[chrIdx];
    }

    BBZoomHeader zoomLevelForScale(double bpPerPixel) {
        BBZoomHeader level = null;
        for (BBZoomHeader zl : this.zoomHeaders) {
            if (zl.reductionLevel < bpPerPixel) {
                level = zl;
                break;
            }
        }
        return level;
    }

}
