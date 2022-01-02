package org.broad.igv.ext.render;

import org.broad.igv.sam.Alignment;

public class ColorByTagValueList implements IColorByTagExtension {

    @Override
    public boolean extendsContext(Object context) {

        // context is tag
        if ( !(context instanceof String) ) {
            return false;
        }
        final String tag = (String)context;

        // must contain list prefix
        return tag.contains(":");
    }

    @Override
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
