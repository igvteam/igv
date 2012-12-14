package org.broad.igv.feature;

import org.broad.igv.feature.*;


/**
 * @author jrobinso
 *         Date: 11/29/12
 *         Time: 6:59 PM
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
    private String text;


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

    public void setText(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }
}
