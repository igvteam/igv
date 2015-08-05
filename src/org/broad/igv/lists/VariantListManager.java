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
