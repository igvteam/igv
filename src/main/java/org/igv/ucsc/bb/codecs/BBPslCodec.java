package org.igv.ucsc.bb.codecs;

import org.igv.Globals;
import org.igv.feature.Exon;
import org.igv.feature.PSLRecord;
import org.igv.feature.Strand;
import org.igv.ucsc.bb.BBUtils;
import org.igv.ucsc.bb.BedData;

/**
 * Codec for UCSC "bigPsl" files (bed12+13).  Records are decoded directly into PSLRecords, populated
 * the same way PSLCodec populates them from ascii PSL.
 * <p>
 * Notes on the bigPsl format (see kent/src/hg/utils/pslToBigPsl):
 * - target blocks are always on the + reference strand, as required by BED.  No reversal is needed.
 * - "strand" is the PSL query strand.  oChromStart/oChromEnd/oChromSize/oChromStarts are PSL
 * qStart/qEnd/qSize/qStarts, chromSize is tSize, match/misMatch/repMatch/nCount are copied verbatim.
 * - the insert counts (qNumInsert, qBaseInsert, tNumInsert, tBaseInsert) are not stored and are
 * computed from the blocks (kent pslComputeInsertCounts).
 * - for protein alignments (seqType == 2) blockSizes are in bases (3x residues) while query
 * coordinates are in residues.
 * <p>
 * oStrand, oSequence, oCDS, seqType and thickStart/thickEnd are not carried onto the record.
 */
public class BBPslCodec extends BBRecordCodec {

    // Indices into "restOfFields" (which starts at the 4th bed column) for the standard bed12 fields
    private static final int NAME = 0, STRAND = 2, BLOCK_COUNT = 6, BLOCK_SIZES = 7, CHROM_STARTS = 8;

    private final int oChromStart, oChromEnd, oChromSize, oChromStarts, chromSize, match, misMatch, repMatch,
            nCount, seqType;

    public BBPslCodec(int standardFieldCount, BBUtils.ASTable autosql) {
        super("bigPsl", standardFieldCount, autosql);
        oChromStart = fieldIndex(autosql, "oChromStart");
        oChromEnd = fieldIndex(autosql, "oChromEnd");
        oChromSize = fieldIndex(autosql, "oChromSize");
        oChromStarts = fieldIndex(autosql, "oChromStarts");
        chromSize = fieldIndex(autosql, "chromSize");
        match = fieldIndex(autosql, "match");
        misMatch = fieldIndex(autosql, "misMatch");
        repMatch = fieldIndex(autosql, "repMatch");
        nCount = fieldIndex(autosql, "nCount");
        seqType = fieldIndex(autosql, "seqType");
    }

    @Override
    PSLRecord decodeRecord(BedData data) {

        String[] tokens = Globals.tabPattern.split(data.restOfFields, -1);
        if (tokens.length <= CHROM_STARTS) {
            return null;
        }

        String chr = data.chr;
        int start = data.chromStart;
        int end = data.chromEnd;
        Strand strand = tokens[STRAND].startsWith("+") ? Strand.POSITIVE : Strand.NEGATIVE;

        int blockCount = Integer.parseInt(tokens[BLOCK_COUNT]);
        String[] blockSizes = Globals.commaPattern.split(tokens[BLOCK_SIZES]);
        String[] chromStarts = Globals.commaPattern.split(tokens[CHROM_STARTS]);
        String qStartsString = stringValue(tokens, oChromStarts, "");
        String[] qStarts = Globals.commaPattern.split(qStartsString);
        if (blockSizes.length != blockCount || chromStarts.length != blockCount || qStarts.length != blockCount) {
            return null;
        }

        PSLRecord f = new PSLRecord();
        f.setName(tokens[NAME]);
        f.setChr(chr);
        f.setStart(start);
        f.setEnd(end);
        f.setStrand(strand);

        // Blocks, and insert counts (kent pslComputeInsertCounts).  bigPsl block sizes are in bases, query
        // coordinates are in residues for protein alignments.
        int qDiv = "2".equals(stringValue(tokens, seqType, "0")) ? 3 : 1;
        int qNumInsert = 0, qBaseInsert = 0, tNumInsert = 0, tBaseInsert = 0;
        int prevSize = 0, prevTStart = 0, prevQStart = 0;
        for (int i = 0; i < blockCount; i++) {
            int size = Integer.parseInt(blockSizes[i].trim());
            int tStart = start + Integer.parseInt(chromStarts[i].trim());
            int qStart = Integer.parseInt(qStarts[i].trim());
            f.addExon(new Exon(chr, tStart, tStart + size, strand));
            if (i > 0) {
                int qGap = qStart - (prevQStart + prevSize / qDiv);
                int tGap = tStart - (prevTStart + prevSize);
                if (qGap > 0) {
                    qNumInsert++;
                    qBaseInsert += qGap;
                }
                if (tGap > 0) {
                    tNumInsert++;
                    tBaseInsert += tGap;
                }
            }
            prevSize = size;
            prevTStart = tStart;
            prevQStart = qStart;
        }

        int match = intValue(tokens, this.match);
        int misMatch = intValue(tokens, this.misMatch);
        int repMatch = intValue(tokens, this.repMatch);
        int ns = intValue(tokens, this.nCount);
        int qSize = intValue(tokens, this.oChromSize);

        // Score and description, as in PSLCodec
        float score = (1000.0f * (match + repMatch - misMatch - qNumInsert - tNumInsert)) / qSize;

        f.setMatch(match);
        f.setMisMatch(misMatch);
        f.setRepMatch(repMatch);
        f.setNs(ns);
        f.setQGapCount(qNumInsert);
        f.setQGapBases(qBaseInsert);
        f.setTGapCount(tNumInsert);
        f.setTGapBases(tBaseInsert);
        f.setQSize(qSize);
        f.setTSize(intValue(tokens, this.chromSize));
        f.setqStart(intValue(tokens, this.oChromStart));
        f.setqEnd(intValue(tokens, this.oChromEnd));
        f.setBlockQueryStarts(qStartsString);
        f.setScore(score);

        StringBuffer desc = new StringBuffer();
        desc.append("matches = " + match);
        desc.append("<br>");
        desc.append("mismatches = " + misMatch);
        desc.append("<br>");
        desc.append("repeat matches = " + repMatch);
        desc.append("<br>");
        desc.append("# inserts in query = " + qNumInsert);
        desc.append("<br>");
        desc.append("# inserts in target = " + tNumInsert);
        f.setDescription(desc.toString());

        return f;
    }
}
