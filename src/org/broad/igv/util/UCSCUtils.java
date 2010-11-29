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

package org.broad.igv.util;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: May 30, 2010
 * Time: 10:42:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class UCSCUtils {

    static Set<String> types = new HashSet(Arrays.asList("SINE", "LINE", "LTR", "DNA", "Simple_repeat",
            "Low_complexity", "Satellite", "RNA", "Other", "Unknown", "Uncategorized"));

    public static void main(String[] args) throws IOException {
        String ifile = "/Users/jrobinso/IGV/mm9/mm9_repmsk.txt";
        String output = "/Users/jrobinso/IGV/mm9//repeat_masker";
        String prefix = "mm9_repmask_";
        splitRepeatMasker(ifile, output, prefix);
    }

    public static void splitRepeatMasker(String iFile, String outputDirectory, String prefix) throws IOException {


        BufferedReader br = null;
        Map<String, PrintWriter> pws = new HashMap();

        try {
            br = new BufferedReader(new FileReader(iFile));
            File dir = new File(outputDirectory);
            for (String type : types) {
                File f = new File(dir, prefix + type + ".bed");
                pws.put(type, new PrintWriter(new BufferedWriter(new FileWriter(f))));
            }

            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("#")) continue;
                String[] tokens = nextLine.split("\t");
                String type = getType(tokens[5]);

                // Rerrange columns for legal bed
                pws.get(type).println(tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] + "\t" +
                        tokens[4] + "\t" + tokens[3]);
            }


        }
        finally {
            if (br != null) {
                br.close();
            }
            for (PrintWriter pw : pws.values()) {
                pw.close();
            }
        }

    }


    public static String getType(String s) {

        s = s.replace("?", "");

        if (s.contains("RNA")) {
            return "RNA";
        } else if (s.equals("RC")) {
            return "Other";
        } else if (s.equals("repClass")) {
            return "Other";
        } else if (types.contains(s)) {
            return s;
        } else {
            return "Uncategorized";
        }

    }


}
