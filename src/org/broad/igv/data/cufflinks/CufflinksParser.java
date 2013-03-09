/*
 * Copyright (c) 2007-2013 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.data.cufflinks;

import org.apache.log4j.Logger;
import org.broad.igv.feature.Locus;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Parses various cufflinks (and cuffdiff) output files as described here:
 * <p/>
 * http://cufflinks.cbcb.umd.edu/manual.html
 *
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 2:30 PM
 */
public class CufflinksParser {

    private static Logger log = Logger.getLogger(CufflinksParser.class);

    public static List<CufflinksValue> parse(String path) throws IOException {

        final String s = path.toLowerCase();
        if (s.endsWith("fpkm_tracking")) {
            return parseFpkmTracking(path);
        } else if (s.endsWith("gene_exp.diff") || s.endsWith("cds_exp.diff")) {
            return parseExpDiff(path);
        } else {
            throw new RuntimeException("Unsupported file type: " + path);
        }


    }

    private static List<CufflinksValue> parseFpkmTracking(String path) throws IOException {

        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(path);

            List<CufflinksValue> values = new ArrayList<CufflinksValue>();

            int geneColumn = 4;
            int locusColumn  = 6;
            int fpkmColumn = 6;
            int confLoColumn = 7;
            int confHiColumn = 8;
            int sigColumn = 12;

            String nextLine = br.readLine(); // Header
            String [] tokens = ParsingUtils.TAB_PATTERN.split(nextLine);
            for(int i=0; i<tokens.length; i++) {
                final String tk = tokens[i];
                if(tk.equals("locus")) locusColumn = i;
                else if(tk.equals("gene_short_name")) geneColumn = i;
                else if(tk.equals("FPKM")) fpkmColumn = i;
                else if(tk.equals("FPKM_conf_lo")) confLoColumn = i;
                else if(tk.startsWith("FPKM_conf_hi")) confHiColumn = i;
            }

            while ((nextLine = br.readLine()) != null) {
                try {
                    tokens = ParsingUtils.TAB_PATTERN.split(nextLine);
                    if (tokens.length >= 12) {
                        String locusString = tokens[locusColumn];
                        Locus locus = new Locus(locusString);
                        if (locus != null) {
                            if (locus.getChr() == null) {
                                continue;
                            }
                            float fpkm = Float.parseFloat(tokens[fpkmColumn]);
                            float confLo = Float.parseFloat(tokens[confLoColumn]);
                            float confHi = Float.parseFloat(tokens[confHiColumn]);
                            String gene = tokens[geneColumn];
                            values.add(new FPKMValue(gene, locus.getChr(), locus.getStart() - 1, locus.getEnd(), fpkm, confLo, confHi));
                        }
                    } else {
                        log.info("Unexpected # of columns.  Expected at least 12,  found " + tokens.length);
                    }
                } catch (NumberFormatException e) {
                    log.info("Skipping line: " + nextLine);
                }
            }
            return values;
        } finally {
            if (br != null) br.close();
        }
    }


    private static List<CufflinksValue> parseExpDiff(String path) throws IOException {

        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(path);

            List<CufflinksValue> values = new ArrayList<CufflinksValue>();

            int geneColumn = 1;
            int locusColumn  = 2;
            int xColumn = 6;
            int yColumn = 7;
            int logRatioColumn = 8;
            int sigColumn = 12;

            String nextLine = br.readLine(); // Header
            String [] tokens = ParsingUtils.TAB_PATTERN.split(nextLine);
            for(int i=0; i<tokens.length; i++) {
                final String tk = tokens[i];
                if(tk.equals("locus")) locusColumn = i;
                else if(tk.equals("gene")) geneColumn = i;
                else if(tk.equals("value_1")) xColumn = i;
                else if(tk.equals("value_2")) yColumn = i;
                else if(tk.startsWith("log2(")) logRatioColumn = i;
                else if(tk.equals("significant")) sigColumn = i;
            }

            while ((nextLine = br.readLine()) != null) {
                try {
                    tokens = ParsingUtils.TAB_PATTERN.split(nextLine);
                    if (tokens.length >= sigColumn) {
                        String locusString = tokens[locusColumn];
                        Locus locus = new Locus(locusString);
                        if (locus != null) {
                            if (locus.getChr() == null) {
                                continue;
                            }
                            float logRatio = Float.parseFloat(tokens[logRatioColumn]);

                            if(Float.isInfinite(logRatio) || Float.isNaN(logRatio)) continue;

                            float fpkmX = Float.parseFloat(tokens[xColumn]);
                            float fpkmY = Float.parseFloat(tokens[yColumn]);
                            String gene = tokens[geneColumn] ;
                            String significant = tokens[sigColumn];
                            values.add(new ExpDiffValue(gene, locus.getChr(), locus.getStart() - 1, locus.getEnd(),
                                    logRatio, fpkmX, fpkmY, significant) {
                            });
                        }
                    } else {
                        log.info("Unexpected # of columns.  Expected at least 12,  found " + tokens.length);
                    }
                } catch (NumberFormatException e) {
                    log.info("Skipping line: " + nextLine);
                }
            }
            return values;
        } finally {
            if (br != null) br.close();
        }

    }
}
