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
    public static String ATTR_TI = "ti";
    private static String RG_ATTR_PL = "PL";
    private static String RG_ATTR_MC = "mc";
    private static String RG_ATTR_PL_ULTIMA = "ULTIMA";

    public enum UltimaFileFormat {
        NON_FLOW,
        BASE_TI,
        BASE_TP
    }

    static public boolean isFlow(SAMRecord record) {

        return isUltimaFlowReadGroup(record.getReadGroup())
                && (record.hasAttribute(ATTR_TP) || record.hasAttribute(ATTR_TI));
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

    static public UltimaFileFormat getUltimaFileVersion(Alignment alignment)
    {
        // check for preconditions for a flow file.  note that we extend only SAMAlignment instances
        if ( !(alignment instanceof SAMAlignment) ) {
            return FlowUtil.UltimaFileFormat.NON_FLOW;
        }
        final SAMAlignment samAlignment = (SAMAlignment)alignment;
        if ( !FlowUtil.isUltimaFlowReadGroup(samAlignment.getRecord().getReadGroup()) ) {
            return FlowUtil.UltimaFileFormat.NON_FLOW;
        }

        if ( alignment.getAttribute(FlowUtil.ATTR_TI) != null )
            return FlowUtil.UltimaFileFormat.BASE_TI;
        else if ( alignment.getAttribute(FlowUtil.ATTR_TP) != null )
            return FlowUtil.UltimaFileFormat.BASE_TP;
        else
            return FlowUtil.UltimaFileFormat.NON_FLOW;
    }
}
