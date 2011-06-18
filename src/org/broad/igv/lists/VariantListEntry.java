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

package org.broad.igv.lists;

import java.util.List;

/**
 * An experimental class to represent a list of variants, and associated samples.
 *
 * @author jrobinso
 * @date Apr 27, 2011
 */
public class VariantListEntry {

    String name;
    String listLabel;
    String chr;
    int position;   // <= 1-based
    String[] samples;

    public VariantListEntry(String chr, int position, String[] samples) {
        this.chr = chr;
        this.position = position;
        this.samples = samples;
        this.name = chr + ":" + position;

        StringBuffer buf = new StringBuffer();
        buf.append(name);
        buf.append(" ");
        boolean first = true;
        for (String s : samples) {
            if (!first) buf.append(",");
            buf.append(s);
            first = false;
        }
        listLabel = buf.toString();
    }

    public String toString() {
        return listLabel;
    }
}
