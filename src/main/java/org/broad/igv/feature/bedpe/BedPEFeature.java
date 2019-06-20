package org.broad.igv.feature.bedpe;

import java.awt.*;
import java.util.Map;

/**
 * Created by jrobinso on 6/29/18.
 */
public class BedPEFeature {

    String chr1;
    int start1;
    int end1;
    String chr2;
    int start2;
    int end2;
    String name;
    double score;
    Color color;
    int thickness = 1;
    String type;
    Map<String, String> attributes;


    public int getStart() {
        return Math.min(start1, start2);
    }

    public int getEnd() {
        return Math.max(end1, end2);
    }

    public boolean isSameChr() {
        return chr1.equals(chr2);
    }

}
