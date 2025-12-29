package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;

/**
 * @author jrobinso
 */
public class CompositeScore implements LocusScore {

    float[] data;
    String[] probes;
    float value;
    private int start;
    private int end;
    private WindowFunction windowFunction;


    public CompositeScore(int start, int end, float value, float[] data, String[] probes, WindowFunction windowFunction) {
        this.start = start;
        this.end = end;
        this.value = value;
        this.probes = probes;
        this.data = data;

        // Only keep 5 representative values
        if (data.length > 5) {
            float[] temp = new float[5];
            System.arraycopy(data, 0, temp, 0, 5);
            this.data = temp;

            if (probes != null && probes.length > 5) {
                String[] temp2 = new String[5];
                System.arraycopy(probes, 0, temp2, 0, 5);
                this.probes = temp2;
            }
        }

    }

    public CompositeScore(CompositeScore sc) {
        this.start = sc.start;
        this.end = sc.end;
        this.value = sc.value;
        this.probes = sc.probes;
        this.data = sc.data;
    }


    public CompositeScore copy() {
        return new CompositeScore(this);
    }


    public float getScore() {
        return value;
    }

    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        if (data == null) {
            return "";
        }
        StringBuffer buf = new StringBuffer();
        buf.append("Composite value = " + value + " (" + windowFunction + ")<br>");
        buf.append("-------------------------------<br>");
        for (int j = 0; j < data.length; j++) {
            buf.append(String.valueOf(data[j]));
            String probe = (probes != null && j < probes.length) ? probes[j] : null;
            if (probe != null && probe.length() > 0) {
                buf.append("&nbsp;(");
                buf.append(probes[j]);
                buf.append(")");
            }
            buf.append("<br>");
        }
        if (data.length >= 5) {
            buf.append("...<br>");
        }
        return buf.toString();
    }

    public String getChr() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getContig() {
        return null;
    }

    public int getStart() {
        return start;
    }


    public void setStart(int start) {
        this.start = start;
    }


    public int getEnd() {
        return end;
    }

    /**
     * Method description
     *
     * @param end
     */
    public void setEnd(int end) {
        this.end = end;
    }

    public int getExtendedStart() {
        return getStart();
    }

    public int getExtendedEnd() {
        return getEnd();
    }
}
