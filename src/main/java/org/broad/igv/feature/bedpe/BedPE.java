package org.broad.igv.feature.bedpe;

public interface BedPE {

    public BedPEFeature get();
    public String getChr();
    public int getStart();
    public int getEnd();
    public double getScore();
    public boolean isSameChr();

}
