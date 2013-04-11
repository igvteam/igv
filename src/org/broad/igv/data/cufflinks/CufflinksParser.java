/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
        CufflinksCodec codec;
        if (s.endsWith("fpkm_tracking")) {
            codec = new FpkmTrackingCodec();
        } else if (s.endsWith("gene_exp.diff") || s.endsWith("cds_exp.diff")) {
            codec = new ExpDiffCodec();
        } else {
            throw new RuntimeException("Unsupported file type: " + path);
        }
        return parse(codec, path);
    }

    private static List<CufflinksValue> parse(CufflinksCodec codec, String path) throws IOException {
        BufferedReader br = null;
        List<CufflinksValue> values = new ArrayList<CufflinksValue>();

        try {
            br = ParsingUtils.openBufferedReader(path);

            String nextLine = br.readLine(); // Header
            codec.readHeader(nextLine);

            while ((nextLine = br.readLine()) != null) {
                CufflinksValue value = codec.decode(ParsingUtils.TAB_PATTERN.split(nextLine));
                if(value == null){
                    log.info("Skipping line " + nextLine);
                    continue;
                }
                values.add(value);
            }
        }finally {
            if (br != null) br.close();
        }
        return values;
    }

    private static class FpkmTrackingCodec implements CufflinksCodec{
        int geneColumn = 4;
        int locusColumn  = 6;
        int fpkmColumn = 6;
        int confLoColumn = 7;
        int confHiColumn = 8;
        int sigColumn = 12;

        public void readHeader(String headerLine){
            String[] tokens = ParsingUtils.TAB_PATTERN.split(headerLine);
            for(int i=0; i<tokens.length; i++) {
                final String tk = tokens[i];
                if(tk.equals("locus")) locusColumn = i;
                else if(tk.equals("gene_short_name")) geneColumn = i;
                else if(tk.equals("FPKM")) fpkmColumn = i;
                else if(tk.equals("FPKM_conf_lo")) confLoColumn = i;
                else if(tk.startsWith("FPKM_conf_hi")) confHiColumn = i;
            }
        }

        @Override
        public CufflinksValue decode(String[] tokens) {
            if (tokens.length >= 12) {
                String locusString = tokens[locusColumn];
                if (locusString == null) return null;

                Locus locus = new Locus(locusString);
                if(locus.getChr() == null) return null;

                float fpkm = Float.parseFloat(tokens[fpkmColumn]);
                float confLo = Float.parseFloat(tokens[confLoColumn]);
                float confHi = Float.parseFloat(tokens[confHiColumn]);
                String gene = tokens[geneColumn];
                return new FPKMValue(gene, locus.getChr(), locus.getStart() - 1, locus.getEnd(), fpkm, confLo, confHi);
            } else {
                log.info("Unexpected # of columns.  Expected at least 12,  found " + tokens.length);
                return null;
            }
        }
    }

    private static class ExpDiffCodec implements CufflinksCodec{
        int geneColumn = 1;
        int locusColumn  = 2;
        int xColumn = 6;
        int yColumn = 7;
        int logRatioColumn = 8;
        int sigColumn = 12;

        @Override
        public void readHeader(String headerLine){
            String[] tokens = ParsingUtils.TAB_PATTERN.split(headerLine);
            for(int i=0; i<tokens.length; i++) {
                final String tk = tokens[i];
                if(tk.equals("locus")) locusColumn = i;
                else if(tk.equals("gene")) geneColumn = i;
                else if(tk.equals("value_1")) xColumn = i;
                else if(tk.equals("value_2")) yColumn = i;
                else if(tk.startsWith("log2(")) logRatioColumn = i;
                else if(tk.equals("significant")) sigColumn = i;
            }
        }

        @Override
        public CufflinksValue decode(String[] tokens) {
            if (tokens.length >= sigColumn) {
                String locusString = tokens[locusColumn];
                if (locusString == null) return null;

                Locus locus = new Locus(locusString);
                if(locus.getChr() == null) return null;

                String logRatioStr = tokens[logRatioColumn];
                float logRatio = Float.parseFloat(logRatioStr);

                if (Float.isInfinite(logRatio) || Float.isNaN(logRatio)){
                    log.info("LogRatio " + logRatioStr + " cannot be parsed as a float");
                    logRatio = Float.NaN;
                }

                float fpkmX = Float.parseFloat(tokens[xColumn]);
                float fpkmY = Float.parseFloat(tokens[yColumn]);
                String gene = tokens[geneColumn];
                String significant = tokens[sigColumn];
                return new ExpDiffValue(gene, locus.getChr(), locus.getStart() - 1, locus.getEnd(),
                        logRatio, fpkmX, fpkmY, significant);
            } else {
                log.info("Unexpected # of columns.  Expected at least 12,  found " + tokens.length);
                return null;
            }
        }

    }

    private interface CufflinksCodec {
        /**
         * Read the header information, setting column values as necessary
         * @param headerLine
         */
        void readHeader(String headerLine);

        /**
         * Parse a single line of data from a cufflinks file
         * @param tokens
         * @return
         */
        CufflinksValue decode(String[] tokens);

    }
}
