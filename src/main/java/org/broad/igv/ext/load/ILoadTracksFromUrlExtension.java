package org.broad.igv.ext.load;

import org.broad.igv.ext.IExtension;
import org.broad.igv.util.ResourceLocator;

import java.util.Collection;

public interface ILoadTracksFromUrlExtension extends IExtension {

    public Collection<ResourceLocator> locatorsForUrl(final String url, final String indexUrl);
}
