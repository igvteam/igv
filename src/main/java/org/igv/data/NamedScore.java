package org.igv.data;

import org.igv.track.WindowFunction;

/**
 * @author jrobinso
 * @date May 19, 2011
 */
public class NamedScore extends BasicScore {

    private String probe;

    public NamedScore(int start, int end, float score, String probe) {
        super(start, end, score);
        this.probe = probe;
    }

    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        StringBuffer buf = new StringBuffer();
        buf.append("Value: " + score);
        if(probe != null && probe.length() > 0) {
            buf.append("&nbsp;(");
            buf.append(probe);
            buf.append(")");
        }
        buf.append("<br>");

        return buf.toString();
    }

}
