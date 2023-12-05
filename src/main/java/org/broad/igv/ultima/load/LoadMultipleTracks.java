package org.broad.igv.ultima.load;

import org.broad.igv.util.GoogleUtils;
import org.broad.igv.ui.action.LoadFromURLMenuAction;
import org.broad.igv.util.ResourceLocator;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

public class LoadMultipleTracks {

    final static String SPLIT_REGEX = "\\s+";

    public boolean handlesUrl(final String url) {

        // url must not be a session
        if (url.endsWith(".xml") || url.endsWith(".session")) {
            return false;
        }

        // url should not be an s3:// url
        if ( url.startsWith("s3://") ) {
            return false;
        }

        // must contain multiple
        String[]        toks = url.split(SPLIT_REGEX);
        if ( toks.length <= 1 )
            return false;

        // if here, url can be split into individual urls
        return true;
    }

    public Collection<ResourceLocator> locatorsForUrl(final String url, final String indexUrl) {

        List<ResourceLocator>   locators = new LinkedList<>();
        for ( String tok : url.split(SPLIT_REGEX) ) {

            ResourceLocator rl = new ResourceLocator(LoadFromURLMenuAction.mapURL(tok));

            if (indexUrl != null) {

                final String index = indexUrl.trim();
                if (GoogleUtils.isGoogleDrive(index) || GoogleUtils.isGoogleDrive(index)) {
                    LoadFromURLMenuAction.enableGoogleMenu();
                }

                rl.setIndexPath(index);
            }

            locators.add(rl);
        }

        return locators;
    }
}
