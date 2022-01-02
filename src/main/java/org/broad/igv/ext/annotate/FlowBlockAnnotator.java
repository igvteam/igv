package org.broad.igv.ext.annotate;

import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang3.StringUtils;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.ByteSubarray;
import org.broad.igv.sam.SAMAlignment;

public class FlowBlockAnnotator implements IAlignmentBlockAnnotationExtension {

    private static final String KEY_ATTR = "ti,tp";

    @Override
    public boolean extendsContext(final Object context) {

        // must be a block with qualities
        if ( !(context instanceof AlignmentBlock) ) {
            return false;
        }
        final AlignmentBlock block = (AlignmentBlock)context;

        // sanity: bases and qualities
        if ( (block.getBases() == null) || (block.getQualities() == null) ) {
            return false;
        }


        return true;
    }

    @Override
    public void appendBlockQualityAnnotation(SAMAlignment samAlignment, AlignmentBlock block, StringBuffer buf) {

        if ( isFlow(samAlignment.getRecord()) ) {
            buf.append(" @ QV " + qualsAsString(block.getQualities()) + attrAsString(samAlignment, block, "ti,tp", -1));
        }
    }

    @Override
    public void appendBlockAttrAnnotation(SAMAlignment samAlignment, AlignmentBlock block, int offset, StringBuffer buf) {

        if ( isFlow(samAlignment.getRecord()) ) {
            buf.append(attrAsString(samAlignment, block, KEY_ATTR, offset));
        }
    }

    private boolean isFlow(SAMRecord record) {

        // must have one of the key attributes
        for ( String name : KEY_ATTR.split(",") ) {
            if ( record.hasAttribute(name) ) {
                return true;
            }
        }
        return false;
    }

    private int[] attrAsIntegers(SAMAlignment samAlignment, AlignmentBlock block, String name, int offset, StringBuilder sbName) {

        Object  value = null;
        for ( String name1 : name.split(",") ) {
            value = samAlignment.getRecord().getAttribute(name1);
            if ( value != null ) {
                name = name1;
                if ( sbName != null )
                    sbName.append(name);
                break;
            }
        }
        if ( value == null )
            return new int[0];
        byte[]  arr = (byte[])value;

        if ( offset >= 0 ) {
            int[]   integers = new int[1];
            int     start = block.getIndexOnRead();
            integers[0] = (int)(arr[start + offset]);
            return integers;

        } else {
            int     start = block.getIndexOnRead();
            int     length = block.getLength();
            int[]   integers = new int[length];
            for ( int ofs = start ; ofs < (start + length) ; ofs++ ) {
                integers[ofs - start] = arr[ofs];
            }
            return integers;
        }
    }

    private String attrAsString(SAMAlignment samAlignment, AlignmentBlock block, String name, int offset) {

        StringBuilder   sb = new StringBuilder(" ");
        int[]       integers = attrAsIntegers(samAlignment, block, name, offset, sb);
        if ( integers == null )
            return "";

        sb.append(" ");
        sb.append(StringUtils.join(integers, ','));

        String  r = sb.toString();
        if ( r.length() > 40 )
            r = r.substring(0, 40) + "...";

        return r;
    }

    static public String qualsAsString(ByteSubarray quals) {

        StringBuilder       sb = new StringBuilder();
        final int           digestLength = 25;
        boolean             isDigest = quals.length >= (digestLength * 2);

        // translate quals back into ASCII
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final byte q = quals.getByte(i);
            if ( !isDigest || (i < digestLength) || (i > quals.length - digestLength) ) {
                if ( i > 0 )
                    sb.append(",");
                if (q == 255)
                    sb.append('?');
                else {
                    sb.append(Integer.toString(q));
                }
            }
            if ( isDigest && (i == digestLength) )
                sb.append("...");
            i++;
        }

        return sb.toString();
    }
}
