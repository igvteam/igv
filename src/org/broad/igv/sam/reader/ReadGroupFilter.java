/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.sam.reader;

import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

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
