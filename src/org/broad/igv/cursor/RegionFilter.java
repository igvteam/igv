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

package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 2/9/14
 *         Time: 9:30 PM
 */
public class RegionFilter {

    public List<Clause> getClauses() {
        return clauses;
    }

    public enum Pred {AND, OR}

    public enum Condition {exists, GT, LT}

    Pred pred;
    List<Clause> clauses;

    public RegionFilter() {
        this(Pred.AND);
    }

    public RegionFilter(Pred pred) {
        this.pred = pred;
        this.clauses = new ArrayList<Clause>();
    }

    public boolean pass(CursorRegion region) {

        if(clauses.isEmpty()) return true;

        if (pred == Pred.AND) {
            for (Clause c : clauses) {
                if (!c.pass(region)) return false;
            }
            return true;
        } else {
            for (Clause c : clauses) {
                if (c.pass(region)) return true;
            }
            return false;
        }

    }

    public static class Clause {

        Condition condition;
        double threshold;
        CursorTrack track;

        boolean pass(CursorRegion region) {

            List<BasicFeature> features = track.getFeatures(region.getChr());
            int lf = track.getLongestFeatureLength(region.getChr());
            double score = region.getScore(track, CursorModel.frameBPWidth);

            if (score < 0) return false;
            else {
                if (condition == Condition.GT) return score > threshold;
                else if (condition == Condition.LT) return score < threshold;
                else return true;
            }

        }
    }

}
