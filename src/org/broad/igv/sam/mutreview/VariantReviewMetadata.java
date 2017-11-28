package org.broad.igv.sam.mutreview;

/**
 * Created by jrobinso on 11/13/17.
 */
public class VariantReviewMetadata {

    public VariantReviewMetadata() {
        id = System.currentTimeMillis();
    }

    public long id;
    public int score;
    public String scoreString;
    public String userId;
    public String userEmail;
    public String chrom;
    public int pos;
    public int windowSize;
    public int readDepth;
    public char ref;
    public int refCount;
    public int refCountPos;
    public int refCountNeg;
    public String alt;
    public String altCount;
    public String altCountPos;
    public String altCountNeg;
    public int delCount;

}
