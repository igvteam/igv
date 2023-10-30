package org.broad.igv.ucsc;

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
public class BBIfile {
}
