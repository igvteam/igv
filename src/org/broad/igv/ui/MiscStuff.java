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

/**
 *  Misc globally visible parameters.  These need a home.
 *
 */
package org.broad.igv.ui;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;

/**
 * Placeholder for miscelleneous stuff that hasn't found a home yet.
 */
public class MiscStuff {

    private static Logger log = Logger.getLogger(MiscStuff.class);


    /**
     * Object representing the complete set of genes for the selected genomeId
     */
    //private static GeneManager geneData;
    /**
     * Special feature track.  Always present and always first in feature list.
     * Created when genomeId is loaded.
     */
    private static boolean showMissingData = false;

    static {
        String showMissingDataPref =
                PreferenceManager.getInstance().get(
                        PreferenceManager.SHOW_MISSING_DATA_KEY, "false");

        if (showMissingDataPref != null && !showMissingDataPref.trim().equals("")) {
            showMissingData = Boolean.parseBoolean(showMissingDataPref);
        }
    }

    public static void setShowMissingDataEnabled(boolean value) {
        showMissingData = value;
    }

    public static boolean isShowMissingDataEnabled() {
        return showMissingData;
    }

    public Class getRendererByType(String type) {

        Class renderer = null;

        if (type.equals(type)) {

        }
        return renderer;
    }
}
