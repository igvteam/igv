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


package org.broad.igv.feature;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Class description
 *
 * @author Enter your name here...
 * @version Enter version here..., 08/10/16
 */
public class CytoBandFileParser {

    private static Logger logger = Logger.getLogger(CytoBandFileParser.class);


    /**
     * Validate cytoband file
     *
     * @param reader
     * @param filename
     * @return
     */
    public static boolean isValid(AsciiLineReader reader, String filename) {

        if (reader == null) {
            return false;
        }

        try {
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String chr = tokens[0].trim();
                Cytoband cytoData = new Cytoband(chr);
                parseData(tokens, cytoData);
            }
            return true;

        } catch (Exception e) {
            logger.error("Invalid Cytoband file data : file=" + filename, e);
            return false;
        }
    }

    /**
     * Method description
     *
     * @param reader
     * @return
     */
    public static LinkedHashMap<String, List<Cytoband>> loadData(BufferedReader reader) {

        LinkedHashMap<String, List<Cytoband>> dataMap = new LinkedHashMap<String, List<Cytoband>>();
        try {

            String nextLine;
            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                String[] data = nextLine.split("\t");
                String chr = data[0].trim();
                List<Cytoband> cytobands = dataMap.get(chr);
                if (cytobands == null) {
                    cytobands = new ArrayList<Cytoband>();
                    dataMap.put(chr, cytobands);
                }
                Cytoband cytoData = new Cytoband(chr);
                parseData(data, cytoData);
                cytobands.add(cytoData);
            }

            reader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
        return dataMap;

    }

    private static void parseData(String[] tokens, Cytoband cytoData) {

        cytoData.setStart(Integer.parseInt(tokens[1].trim()));
        cytoData.setEnd(Integer.parseInt(tokens[2].trim()));
        if (tokens.length > 3) {
            cytoData.setName(tokens[3]);
        }
        if (tokens.length > 4) {
            if (tokens[4].equals("acen")) {
                cytoData.setType('c');
            } else {
                cytoData.setType(tokens[4].charAt(1));
                if (cytoData.getType() == 'p') {
                    cytoData.setStain(Short.parseShort(tokens[4].substring(4).trim()));
                }
            }
        }

    }

    private static String parseChromosome(String[] data) {
        String chr = data[0].substring(3);
        int underscore = chr.indexOf('_');
        if (underscore > 0) {
            chr = chr.substring(0, underscore);
        }
        return chr;
    }
}
