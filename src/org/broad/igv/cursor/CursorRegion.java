package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.EncodePeakFeature;
import org.broad.igv.feature.FeatureUtils;
import org.broad.tribble.Feature;

import java.util.List;

/**
 * @author jrobinso
 *         Date: 1/14/14
 *         Time: 2:04 PM
 */
public class CursorRegion implements Feature {

    String chr;  //
    int location;  // Genomic location in BP.  Zero based coords
    private int MAGIC_NUMBER_CHANGE_ME = 500;

    public CursorRegion(String chr, int location) {
        this.chr = chr;
        this.location = location;
    }

    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public int getStart() {
        return location;
    }

    @Override
    public int getEnd() {
        return location + 1;
    }

    public int getLocation() {
        return location;
    }

    public double getScore(List<BasicFeature> features) {

        double score = 0;

        int bpStart = location - MAGIC_NUMBER_CHANGE_ME;
        int bpEnd = location + MAGIC_NUMBER_CHANGE_ME;

        int i0 = FeatureUtils.getIndexBefore(location - MAGIC_NUMBER_CHANGE_ME, features);
        if (i0 > 0) {

            // Just use max
            for (int fIdx = i0; fIdx < features.size(); fIdx++) {
                BasicFeature f = features.get(fIdx);
                if (f.getStart() > bpEnd) break;
                else if (f.getEnd() < bpStart) continue;
                else {
                    // score = Math.max(f instanceof EncodePeakFeature ? ((EncodePeakFeature) f).getSignal() : f.getScore(), score);
                    score = Math.max(f.getScore(), score);
                }
            }

            // Do a weighted average
//            int count = 0;
//            for (int fIdx = i0; fIdx < features.size(); fIdx++) {
//
//                EncodePeakFeature f = features.get(fIdx);
//                if (f.getStart() > bpEnd) break;
//                else if (f.getEnd() < bpStart) continue;
//                else {
//                    count++;
//                    double weight = Math.min(1, ((double) f.getLength()) / (2*MAGIC_NUMBER_CHANGE_ME));
//                    score +=  weight * f.getScore();
//                }
//            }
//            if(count > 0) score /= count;
        }
        return score;
    }
}
