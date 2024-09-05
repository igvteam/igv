package org.broad.igv.ultima.render;

import org.broad.igv.sam.Alignment;

public class ColorByTagValueList {

    public boolean handlesTag(final String tag) {


        // must contain list prefix
        return tag.contains(":");
    }

    public Object getValueForColorByTag(Alignment alignment, String tag) {

        // parse tag as a tag followed by a series of values, seperated by a column :
        String  toks[] = tag.split(":");
        Object  value = alignment.getAttribute(toks[0]);
        if ( value == null )
            return value;

        // scan for contained values
        String  strValue = value.toString();
        for ( int i = 1 ; i < toks.length ; i++ ) {
            if (strValue.contains(toks[i]))
                return Integer.toString(i);
        }

        // if here, not found
        return "0";
    }
}
