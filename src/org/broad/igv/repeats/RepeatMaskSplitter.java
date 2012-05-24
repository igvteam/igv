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
package org.broad.igv.repeats;

import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Properties;

/**
 * Splits a repeat mask file downloaded from UCSC into multiple files,  one per repeat class.
 * Assumes downloaded columns as follows (use table browser, and "select columns" option
 * <p/>
 * genoName  genoStart  genoEnd  strand  repName repClass repFamily
 * <p/>
 * Assumes file is sorted by chromosome
 *
 * @author jrobinso
 */
public class RepeatMaskSplitter {

    public static void main(String[] args) {
        File file = new File("/Users/jrobinso/Downloads/Repeats/RepMask_3.2.7_hg18.tab");
        split(file);
    }

    public static void split(File file) {

        int binCol = 0;
        int millDivCol = 2;
        int millDelCol = 3;
        int millInsCol = 4;
        int chrCol = 5;
        int startCol = 6;
        int endCol = 7;
        int strandCol = 9;
        int namCol = 10;
        int classCol = 11;
        int famCol = 12;

        Map<String, LinkedHashMap<String, String>> fileMappings = new HashMap();

        AsciiLineReader reader = null;
        HashMap<String, PrintWriter> writers = new HashMap();
        try {
            String lastChr = "";
            reader = new AsciiLineReader(new FileInputStream(file));
            // Skip header
            reader.readLine();
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                String chr = tokens[chrCol];
                if (!chr.equals(lastChr)) {
                    closeWriters(writers);
                }
                lastChr = chr;

                String repClass = tokens[classCol];
                if (repClass.contains("?")) {
                    continue;
                }
                String fileKey = chr + "." + repClass;

                // Get or create file writer for the class + chr combination
                PrintWriter pw = writers.get(fileKey);
                if (pw == null) {

                    File dir = new File(file.getParent(), repClass);
                    if (!dir.exists()) {
                        dir.mkdir();
                        fileMappings.put(repClass, new LinkedHashMap<String, String>());
                    }

                    Map<String, String> fMap = fileMappings.get(repClass);
                    String fn = fileKey + ".bed";
                    fMap.put(chr, fn);

                    File outputFile = new File(dir, fn);
                    pw = new PrintWriter(new FileWriter(outputFile));
                    writers.put(fileKey, pw);
                }

                String name = "Repeat " + tokens[4] + ", family " + tokens[6];

                pw.print(chr);
                pw.print("\t");
                pw.print(Integer.parseInt(tokens[startCol]));
                pw.print("\t");
                pw.print(Integer.parseInt(tokens[endCol]));
                pw.print("\t");
                pw.print(name);
                pw.print("\t");
                pw.print(tokens[strandCol]);
                pw.println();

            }

            // Ouput filemapping files
            for (Map.Entry<String, LinkedHashMap<String, String>> entry : fileMappings.entrySet()) {
                String repClass = entry.getKey();
                File dir = new File(file.getParent(), repClass);
                File listFile = new File(dir, repClass + "_files.list.txt");
                Properties props = new Properties();
                for (Map.Entry<String, String> entry2 : entry.getValue().entrySet()) {
                    props.put(entry2.getKey(), entry2.getValue());
                }
                FileOutputStream os = new FileOutputStream(listFile);
                props.store(os, "");
                os.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            reader.close();
            closeWriters(writers);
        }

    }

    private static void closeWriters(HashMap<String, PrintWriter> writers) {
        for (PrintWriter pw : writers.values()) {
            pw.close();
        }
        writers.clear();
    }
}

