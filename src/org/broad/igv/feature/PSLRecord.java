/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
