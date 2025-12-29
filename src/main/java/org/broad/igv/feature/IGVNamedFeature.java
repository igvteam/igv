package org.broad.igv.feature;

/**
 * @author jrobinso
 * @date Sep 16, 2010
 */
public interface IGVNamedFeature extends htsjdk.tribble.NamedFeature {

    default String getDisplayName(String property) {
        return getName();
    }
}
