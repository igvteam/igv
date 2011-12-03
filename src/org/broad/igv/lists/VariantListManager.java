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

import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Apr 27, 2011
 */
public class VariantListManager {

    static List<VariantListEntry> variantList = new ArrayList();

    static Map<String, String> samplePathMap = new HashMap();

    public static void openNavigator(Frame frame) {
        (new VariantListNavigator(frame)).setVisible(true);
    }

    public static void loadVariants(ResourceLocator locator) throws IOException {

        BufferedReader reader = null;

        try {
            reader = ParsingUtils.openBufferedReader(locator.getPath());
            String nextLine;

            // Assumed format chr position  sampleNames
            while ((nextLine = reader.readLine()) != null) {

                if (nextLine.startsWith("#")) {
                    continue;
                }
                String[] tokens = nextLine.split("\\s+");
                if (tokens.length != 3) {
                    System.err.println("Unexpected number of tokens: " + tokens.length);
                }

                String chr = tokens[0];
                int position = Integer.parseInt(tokens[1]);
                String[] samples = tokens[2].split(",");
                variantList.add(new VariantListEntry(chr, position, samples));
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    public static void loadSamplePathMap(ResourceLocator locator) throws IOException {

        BufferedReader reader = null;

        try {
            reader = ParsingUtils.openBufferedReader(locator.getPath());
            String nextLine;

            // Assumed format chr position  sampleNames
            while ((nextLine = reader.readLine()) != null) {

                if (nextLine.startsWith("#")) {
                    continue;
                }
                String[] tokens = nextLine.split("\\s+");
                if (tokens.length != 2) {
                    System.err.println("Unexpected number of tokens: " + tokens.length);
                }

                String sample = tokens[0];
                String path = tokens[1];
                samplePathMap.put(sample, path);
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

}
