package org.broad.igv.ucsc.bb;

public class BBHeader {

    public int magic;                // 4 byte mMagic Number
    public int version;            // 2 byte version ID; currently 3
    public int nZoomLevels;         // 2 byte count of zoom sumary levels
    public long chromTreeOffset;     // 8 byte offset to mChromosome B+ Tree index
    public long fullDataOffset;      // 8 byte offset to unzoomed data dataCount
    public long fullIndexOffset;     // 8 byte offset to R+ Tree index of items
    public int fieldCount;         // 2 byte number of fields in bed. (0 for bigWig)
    public int definedFieldCount;  // 2 byte number of fields that are bed fields
    public long autoSqlOffset;       // 8 byte offset to 0 terminated string with .as spec
    public long totalSummaryOffset;  // 8 byte offset to file summary data block
    public int uncompressBuffSize;  // 4 byte maximum size for decompressed buffer
    public long extensionOffset;
    public int extraIndexCount = 0;
    public long[] extraIndexOffsets;
    public int dataCount;
}
