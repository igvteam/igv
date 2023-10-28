package org.broad.igv.ext;

import org.broad.igv.ext.annotate.FlowBlockAnnotator;
import org.broad.igv.ext.load.LoadMultipleTracks;
import org.broad.igv.ext.render.ColorByTagValueList;
import org.broad.igv.ext.render.FlowIndelRendering;

public class ExtensionManager {

    private static final IExtension EXTENSIONS[] = {
            new LoadMultipleTracks(),
            new FlowBlockAnnotator(),
            new ColorByTagValueList(),
            new FlowIndelRendering()
    };

    public static IExtension getExtentionFor(final Class<? extends IExtension> type, final Object context) {
        for ( IExtension ext : EXTENSIONS ) {
            if ( type.isInstance(ext) ) {
                if ( ext.extendsContext(context) )
                    return ext;
            }
        }
        return null;
    }
}
