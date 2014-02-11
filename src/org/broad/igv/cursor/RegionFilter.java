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
            double score = region.getScore(features, lf, CursorModel.frameBPWidth);

            if (score < 0) return false;
            else {
                if (condition == Condition.GT) return score > threshold;
                else if (condition == Condition.LT) return score < threshold;
                else return true;
            }

        }
    }

}
