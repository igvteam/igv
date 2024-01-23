package org.broad.igv.ultima;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SAMAlignment;

/*
 * misc utilities to handle flow base (Ultima) data
 */
public class FlowUtil {

    public static final String TAG_T0 = "t0";
    public static final double MIN_PROB_DEFAULT = 0.01;
    public static String ATTR_TP = "tp";
    private static String RG_ATTR_PL = "PL";
    private static String RG_ATTR_MC = "mc";
    private static String RG_ATTR_PL_ULTIMA = "ULTIMA";

    static public boolean isFlow(SAMRecord record) {

        return isUltimaFlowReadGroup(record.getReadGroup()) && record.hasAttribute(ATTR_TP);
    }

    static public boolean isFlow(Alignment alignment) {

        return (alignment instanceof SAMAlignment) && isFlow(((SAMAlignment)alignment).getRecord());
    }

    public static boolean isUltimaFlowReadGroup(SAMReadGroupRecord readGroup) {

        // must have a read group to begin with
        if ( readGroup == null )
            return false;

        // modern files have an ultima platform
        if ( RG_ATTR_PL_ULTIMA.equals(readGroup.getAttribute(RG_ATTR_PL)) )
            return true;

        // fall back on the presence of an mc (max-class)
        if ( readGroup.getAttribute(RG_ATTR_MC) != null )
            return true;

        // if here, probably not an ultima flow read group
        return false;
    }
}
