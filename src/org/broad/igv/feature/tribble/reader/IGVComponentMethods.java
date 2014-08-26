/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature.tribble.reader;

import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.AbstractFeatureReader;

import java.io.IOException;
import java.net.URL;

public class IGVComponentMethods extends AbstractFeatureReader.ComponentMethods {

    public boolean isTabix(String resourcePath, String tabxIndex) throws IOException {
        ResourceLocator locator = new ResourceLocator(resourcePath);
        if (tabxIndex == null) {
            if (HttpUtils.isRemoteURL(locator.getPath())) {
                final URL url = new URL(locator.getPath());
                String path = url.getPath();
                String indexPath = path + ".tbi";   // Strip off parameters
                tabxIndex = locator.getPath().replace(path, indexPath);
            } else {
                tabxIndex = locator.getPath() + ".tbi";
            }
        }
        boolean isTabix =  (locator.getPath().endsWith(".gz") || locator.getPath().endsWith(".bgz")) && FileUtils.resourceExists(tabxIndex);
        if(isTabix) {
            locator.setIndexPath(tabxIndex);
        }
        return isTabix;
    }

}
