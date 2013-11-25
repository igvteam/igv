package org.broad.igv.feature;

/**
 * Representation of a feature from an Encode "peak" file
 *
 * TODO Extending BasicFeature is overkill, no exons for example.
 *
 * @author jrobinso
 *         Date: 11/5/13
 *         Time: 1:11 PM
 */
public class EncodePeakFeature extends BasicFeature {

    private int peakPosition = -1;

    public EncodePeakFeature(String chr, int start, int end) {
        super(chr, start, end);
    }

    public int getPeakPosition() {
        return peakPosition;
    }

    public void setPeakPosition(int peakPosition) {
        this.peakPosition = peakPosition;
    }
}
