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

package org.broad.igv.data.rnai;

/**
 * A simple value class representing a hairpin score for a specific batch
 *
 * @author jrobinso
 */
public class RNAIHairpinValue {

    private String name;
    private float scoreMean;
    private float scoreSTD;
    private boolean hasScoreSTD;

    public RNAIHairpinValue(String name, float scoreMean) {
        this.name = name;
        this.scoreMean = scoreMean;
        hasScoreSTD = false;
    }

    public RNAIHairpinValue(String name, float scoreMean, float scoreSTD) {
        this.name = name;
        this.scoreMean = scoreMean;
        this.scoreSTD = scoreSTD;
        hasScoreSTD = true;
    }

    public String getName() {
        return name;
    }

    public float getScoreMean() {
        return scoreMean;
    }

    public float getScoreSTD() {
        return scoreSTD;
    }

    public boolean hasScoreSTD() {
        return hasScoreSTD;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final RNAIHairpinValue other = (RNAIHairpinValue) obj;
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 41 * hash + (this.name != null ? this.name.hashCode() : 0);
        hash = 41 * hash + Float.floatToIntBits(this.scoreMean);
        hash = 41 * hash + Float.floatToIntBits(this.scoreSTD);
        return hash;
    }


}

