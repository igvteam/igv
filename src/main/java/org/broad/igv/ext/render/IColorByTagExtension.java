package org.broad.igv.ext.render;

import org.broad.igv.ext.IExtension;
import org.broad.igv.sam.Alignment;

public interface IColorByTagExtension extends IExtension {

    public Object getValueForColorByTag(Alignment alignment, String tag);
}
