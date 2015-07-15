package org.broad.igv.feature;

import org.broad.igv.track.WindowFunction;

/**
 * Representation of a feature from an Encode "peak" file
 * <p/>
 * TODO Extending BasicFeature is overkill, no exons for example.
 *
 * @author jrobinso
 *         Date: 11/5/13
 *         Time: 1:11 PM
 */
public class EncodePeakFeature extends BasicFeature implements SignalFeature {

    private int peakPosition = -1;
    private float signal;
    private float PValue;
    private float QValue;


    public EncodePeakFeature(String chr, int start, int end) {
        super(chr, start, end);
    }

    public int getPeakPosition() {
        return peakPosition;
    }

    public void setPeakPosition(int peakPosition) {
        this.peakPosition = peakPosition;
    }

    public void setSignal(float signal) {
        this.signal = signal;
    }

    public float getSignal() {
        return signal;
    }

    public void setPValue(float PValue) {
        this.PValue = PValue;
    }

    public float getPValue() {
        return PValue;
    }

    public void setQValue(float QValue) {
        this.QValue = QValue;
    }

    public float getQValue() {
        return QValue;
    }

    @Override
    public String getValueString(double position, WindowFunction ignored) {

        StringBuffer desc = new StringBuffer();
        desc.append(super.getValueString(position, ignored));

        desc.append("<br>Signal Value: " + signal);
        desc.append("<br>pValue (-log10): " + PValue);
        desc.append("<br>qValue (-log10): " + QValue);
        if (peakPosition > 0) {
            desc.append("<br>Peak: " + (peakPosition + 1));
        }
        return desc.toString();
    }
}
