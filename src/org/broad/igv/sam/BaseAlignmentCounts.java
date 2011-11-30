package org.broad.igv.sam;

import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.panel.ReferenceFrame;

/**
 * @author Jim Robinson
 * @date 11/29/11
 */
abstract public class BaseAlignmentCounts implements AlignmentCounts {

    static char[] nucleotides = {'a', 'c', 'g', 't', 'n'};


    public String getValueStringAt(int pos) {

        StringBuffer buf = new StringBuffer();
        int totalCount = getTotalCount(pos);
        buf.append("Total count: " + totalCount + "<br>");
        for (char c : nucleotides) {
            int negCount = getNegCount(pos, (byte) c);
            int posCount = getPosCount(pos, (byte) c);
            int count = negCount + posCount;
            int percent = (int) Math.round(((float) count) * 100 / totalCount);
            char cU = Character.toUpperCase(c);
            buf.append(cU + "      : " + count);
            if (count == 0) {
                buf.append("<br>");
            } else {
                buf.append("  (" + percent + "%,     " + posCount + "+,   " + negCount + "- )<br>");
            }
        }

        int delCount = getDelCount(pos);
        if (delCount > 0) {
            buf.append("DEL: " + delCount);
        }
        int insCount = getInsCount(pos);
        if (insCount > 0) {
            buf.append("INS: " + insCount);
        }

        return buf.toString();

    }

}
