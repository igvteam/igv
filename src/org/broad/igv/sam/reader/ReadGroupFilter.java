/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.sam.reader;

import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Oct 1, 2009
 * Time: 7:09:37 PM
 * To change this template use File | Settings | File Templates.
 */

public class ReadGroupFilter {


    private Set<String> filteredReadGroups;

    private ReadGroupFilter(Set<String> filteredReadGroups) {
        this.filteredReadGroups = filteredReadGroups;
    }

    public boolean filterAlignment(Alignment alignment) {
        return filteredReadGroups.contains(alignment.getReadGroup());
    }


    static Map<String, ReadGroupFilter> filterCache = new HashMap();

    public static synchronized ReadGroupFilter getFilter() {

        PreferenceManager samPrefs = PreferenceManager.getInstance();

        if (samPrefs.getAsBoolean(PreferenceManager.SAM_FILTER_ALIGNMENTS)) {

            String filterURL = samPrefs.get(PreferenceManager.SAM_FILTER_URL);

            ReadGroupFilter filter = filterURL == null ? null : filterCache.get(filterURL);

            if (filter == null && filterURL != null && filterURL.trim().length() > 0) {
                Set<String> readGroups = new HashSet();
                AsciiLineReader reader = null;
                try {
                    reader = ParsingUtils.openAsciiReader(new ResourceLocator(filterURL));
                    String nextLine;
                    while ((nextLine = reader.readLine()) != null) {
                        readGroups.add(nextLine.trim());
                    }
                    filter = new ReadGroupFilter(readGroups);
                    filterCache.put(filterURL, filter);
                }
                catch (Exception e) {
                    MessageUtils.showMessage("Error reading read filter list: " + e.getMessage());
                }
            }
            return filter;
        }
        return null;

    }

}
