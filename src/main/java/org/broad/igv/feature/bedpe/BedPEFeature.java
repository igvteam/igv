package org.broad.igv.feature.bedpe;

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
    String score;


    public int getStart() {
        return Math.min(start1, start2);
    }

    public int getEnd() {
        return Math.max(end1, end2);
    }

}
