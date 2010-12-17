/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.track.loader;

import org.broad.igv.util.ResourceLocator;

/**
 * An Abstract base class for track loaders.  Provides shared utility functions.
 *
 *
 * @author jrobinso
 * @date Dec 16, 2010
 */
public class AbstractLoader {

    /**
     * Returns a string indicating the resource type.  This is generally the file extension, after removing optional
     * .txt, and .gz extensions, but the type string can also be explicitly set in the ResourceLocator.
     *
     * @param locator
     * @return
     */
    protected String getTypeString(ResourceLocator locator) {
        String typeString = locator.getType();
        if (typeString == null) {
            typeString = locator.getPath().toLowerCase();
            if (!typeString.endsWith("_sorted.txt") &&
                    (typeString.endsWith(".txt") || typeString.endsWith(".gz"))) {
                typeString = typeString.substring(0, typeString.lastIndexOf("."));
            }
        }
        return typeString.toLowerCase();
    }
}
