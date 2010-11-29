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


package org.broad.igv.synteny;

//~--- JDK imports ------------------------------------------------------------

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

/**
 * @author jrobinso
 */
public class SyntenyUtils {


    /**
     * region R:chr1:chr2:D1 chr1 47243 59894 + chr2 111274749 111285190 +
     * anchor A:chr1:3928:chr4:34422 chr1 17429161 17429302 + chr4 140098907 140099048 - 96.7
     *
     * @param file
     */
    public static Map<String, List<SyntenyMapping>> loadMappings(String file, boolean reverse) {

        BufferedReader reader = null;
        Map<String, List<SyntenyMapping>> mappings = new HashMap();

        try {
            reader = new BufferedReader(new FileReader(file));

            String nextLine;
            Pattern.compile("\t");
            while ((nextLine = reader.readLine()) != null) {
                if (!(nextLine.startsWith("region") || nextLine.startsWith("anchor"))) {
                    continue;
                }
                String[] tokens = nextLine.split(" ");
                String type = tokens[0];
                String name = tokens[1];
                String fromChr = tokens[2];
                int fromStart = Integer.parseInt(tokens[3]);
                int fromEnd = Integer.parseInt(tokens[4]);
                String fromStrand = tokens[5];
                String toChr = tokens[6];
                int toStart = Integer.parseInt(tokens[7]);
                int toEnd = Integer.parseInt(tokens[8]);
                String toStrand = tokens[9];

                SyntenyMapping syntenyMapping = null;
                if (reverse) {
                    syntenyMapping = new SyntenyMapping(name, toChr, toStart, toEnd, toStrand,
                            fromChr, fromStart, fromEnd, fromStrand);
                } else {
                    syntenyMapping = new SyntenyMapping(name, fromChr, fromStart, fromEnd, fromStrand,
                            toChr, toStart, toEnd, toStrand);
                }

                List<SyntenyMapping> syntenyMappingList = mappings.get(syntenyMapping.getFromChr());
                if (syntenyMappingList == null) {
                    syntenyMappingList = new ArrayList(1000);
                    mappings.put(syntenyMapping.getFromChr(), syntenyMappingList);
                }
                syntenyMappingList.add(syntenyMapping);

            }


            for (List<SyntenyMapping> syntenyMappingList : mappings.values()) {
                sortMappingList(syntenyMappingList);
            }

            return mappings;

        }
        catch (IOException exception) {
            exception.printStackTrace();
            return null;
        }
        finally {
            if (reader != null) {
                try {

                    reader.close();

                }
                catch (IOException iOException) {
                }
            }
        }
    }


    public static SyntenyMapping getMappingContaining(List<SyntenyMapping> mappings, int fromPosition) {

        int idx = getIndexBefore(mappings, fromPosition);
        for (int i = idx; i < mappings.size(); i++) {
            SyntenyMapping mapping = mappings.get(i);
            if (mapping.containsFromPosition(fromPosition)) {
                return mapping;
            }
        }
        return null;
    }

    public static List<SyntenyMapping> getMappingsOverlapping(List<SyntenyMapping> mappings, int fromStart, int fromEnd) {

        int i1 = -1;
        int i2 = -1;
        int idx = getIndexBefore(mappings, fromStart);
        for (int i = idx; i < mappings.size(); i++) {
            SyntenyMapping mapping = mappings.get(i);
            if (mapping.containsFromPosition(fromStart)) {
                i1 = idx;
                break;
            }
            else if(mapping.getFromStart() > fromStart) {
                i1 = idx;
                break;
            }
        }
        idx = getIndexBefore(mappings, fromEnd);
        for (int i = idx; i < mappings.size(); i++) {
            SyntenyMapping mapping = mappings.get(i);
            if (mapping.containsFromPosition(fromEnd)) {
                i2 = idx;
                break;
            }
            else if(mapping.getFromStart() > fromStart) {
                i2 = Math.max(i1, idx - 1);
                break;
            }
        }

        if(i1 < 0 || i2  < 0) {
            return null;
        }
        return mappings.subList(i1, i2+1);
    }

    /**
     * Sort the feature list by ascending start value
     *
     * @param features
     */
    private static void sortMappingList(List<SyntenyMapping> features) {

        Collections.sort(features, new Comparator<SyntenyMapping>() {

            public int compare(SyntenyMapping o1, SyntenyMapping o2) {

                return (int) (o1.getFromStart() - o2.getFromStart());
            }
        });
    }

    private static int getIndexBefore(List<SyntenyMapping> values, int x) {
        return getIndexBefore(values, x, 0, values.size());
    }

    /**
     * Method description
     *
     * @param values
     * @param x
     * @param leftBound
     * @param rightBound
     * @return
     */
    private static int getIndexBefore(List<SyntenyMapping> values, int x, int leftBound, int rightBound) {

        int idx = (leftBound + rightBound) / 2;

        if ((idx == 0) || (idx == values.size() - 1)) {
            return idx;
        }
        if (values.get(idx).getFromStart() == x) {
            return idx;
        }

        if (values.get(idx).getFromStart() < x) {
            if (values.get(idx + 1).getFromStart() >= x) {
                return idx;
            } else {
                leftBound = idx;
                return getIndexBefore(values, x, leftBound, rightBound);
            }
        } else {    // values[idx] > x
            if (values.get(idx - 1).getFromStart() <= x) {
                return idx - 1;
            } else {
                rightBound = idx;
                return getIndexBefore(values, x, leftBound, rightBound);
            }
        }
    }
}
