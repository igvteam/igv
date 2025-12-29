package org.broad.igv.feature;

import org.broad.igv.feature.*;


/**
 * @author jrobinso
 * Date: 11/29/12
 * Time: 6:59 PM
 * <p>
 * <p>
 * * 1. matches - Number of bases that match that aren't repeats
 * * 2. misMatches - Number of bases that don't match
 * * 3. repMatches - Number of bases that match but are part of repeats
 * * 4. nCount - Number of 'N' bases
 * * 5. qNumInsert - Number of inserts in query
 * * 6. qBaseInsert - Number of bases inserted in query
 * * 7. tNumInsert - Number of inserts in target
 * * 8. tBaseInsert - Number of bases inserted in target
 * * 9. strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
 * * 10. qName - Query sequence name
 * * 11. qSize - Query sequence size
 * * 12. qStart - Alignment start position in query
 * * 13. qEnd - Alignment end position in query
 * * 14. tName - Target sequence name
 * * 15. tSize - Target sequence size
 * * 16. tStart - Alignment start position in target
 * * 17. tEnd - Alignment end position in target
 * * 18. blockCount - Number of blocks in the alignment (a block contains no gaps)
 * * 19. blockSizes - Comma-separated list of sizes of each block
 * * 20. qStarts - Comma-separated list of starting positions of each block in query
 * * 21. tStarts - Comma-separated list of starting positions of each block in target
 * *
 * * 59 9 0 0 1 823 1 96 +- FS_CONTIG_48080_1 1955 171 1062 chr22    47748585 13073589 13073753 2 48,20,  171,1042,  34674832,34674976,* *
 */

public class PSLRecord extends BasicFeature {


    private int tSize;
    private int match;
    private int misMatch;
    private int repMatch;
    private int qNumInsert;
    private int tNumInsert;
    private int qGapCount;
    private int tGapCount;
    private int qSize;
    private int ns;
    private int qGapBases;
    private int tGapBases;

    private String blockQueryStarts;

    private int qStart;
    private int qEnd;

    public void setMatch(int match) {
        this.match = match;
    }

    public void setMisMatch(int misMatch) {
        this.misMatch = misMatch;
    }

    public void setRepMatch(int repMatch) {
        this.repMatch = repMatch;
    }

    public void setQGapCount(int QNumInsert) {
        this.qGapCount = QNumInsert;
    }

    public void setTGapCount(int TNumInsert) {
        this.tGapCount = TNumInsert;
    }

    public void setqStart(int qStart) {
        this.qStart = qStart;
    }

    public void setqEnd(int qEnd) {
        this.qEnd = qEnd;
    }

    public void setQSize(int qSize) {
        this.qSize = qSize;
    }

    public void setNs(int ns) {
        this.ns = ns;
    }

    public void setQGapBases(int qGapBases) {
        this.qGapBases = qGapBases;
    }

    public void setTGapBases(int tGapBases) {
        this.tGapBases = tGapBases;
    }

    public int getTSize() {
        return tSize;
    }

    public int getMatch() {
        return match;
    }

    public int getMisMatch() {
        return misMatch;
    }

    public int getRepMatch() {
        return repMatch;
    }

    public int getQNumInsert() {
        return qNumInsert;
    }

    public int getTNumInsert() {
        return tNumInsert;
    }

    public int getQGapCount() {
        return qGapCount;
    }

    public int getTGapCount() {
        return tGapCount;
    }

    public int getqSize() {
        return qSize;
    }

    public int getNs() {
        return ns;
    }

    public int getQGapBases() {
        return qGapBases;
    }

    public int getTGapBases() {
        return tGapBases;
    }

    public void setBlockQueryStarts(String blockQueryStarts) {
        this.blockQueryStarts = blockQueryStarts;
    }

// 128999999952YourSeq	1289128chr3	9364096243640974310119,	0	36409615,
    public String getText() {
        StringBuffer buffer = new StringBuffer();
        buffer.append(match + "\t");   // 0
        buffer.append(misMatch + "\t");  // 1
        buffer.append(repMatch + "\t");  // 2
        buffer.append(ns + "\t");        // 3
        buffer.append(qGapCount + "\t"); // 4
        buffer.append(qGapBases + "\t"); // 5
        buffer.append(tGapCount + "\t"); // 6
        buffer.append(tGapBases + "\t"); // 7
        buffer.append((strand == Strand.POSITIVE ? '+' : '-') + "\t"); // 8
        buffer.append(name + "\t"); // 9
        buffer.append(qSize + "\t");  // 10
        buffer.append(qStart + "\t"); // 11
        buffer.append(qEnd + "\t");  // 12
        buffer.append(chr + "\t"); // 13
        buffer.append(tSize + "\t"); // 14
        buffer.append(start + "\t");  // 15
        buffer.append(end + "\t");  // 16
        buffer.append(exons.size() + "\t"); // 17

        // block sizes  // 18
        for (Exon exon : exons) {
            buffer.append(exon.getLength() + ",");
        }
        buffer.append("\t");

        buffer.append(blockQueryStarts + "\t");  // 19

        // block target starts  // 20
        for (Exon exon : exons) {
            buffer.append(exon.getStart() + ",");
        }
        buffer.append("\t");


        /*

         */
        return buffer.toString();
    }
}
