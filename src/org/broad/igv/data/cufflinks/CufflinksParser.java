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
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.Locus;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
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

    public static List<? extends CufflinksValue> parse(String path) throws IOException {

        final String s = path.toLowerCase();
        if (s.endsWith("fpkm_tracking")) {
            AsciiFeatureCodec<FPKMValue> codec = new FpkmTrackingCodec(path);
            return parse(codec, path);
        } else if (s.endsWith("gene_exp.diff") || s.endsWith("cds_exp.diff")) {
            AsciiFeatureCodec<ExpDiffValue> codec = new ExpDiffCodec(path);
            return parse(codec, path);
        } else {
            throw new RuntimeException("Unsupported file type: " + path);
        }

    }

    private static <T extends CufflinksValue> List<T> parse(AsciiFeatureCodec<T> codec, String path) throws IOException {
        List<T> values = new ArrayList<T>();
        AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(path, codec, false);
        Iterator<T> iter = reader.iterator();
        while(iter.hasNext()){
            values.add(iter.next());
        }
        return values;
    }

    private static abstract class CufflinksCodec<T extends CufflinksValue> extends AsciiFeatureCodec<T>{

        String path;

        CufflinksCodec(Class<T> clazz, String path){
            super(clazz);
            this.path = path;
        }

        protected abstract Object readHeader(String[] tokens);

        @Override
        public Object readHeader(LineReader reader){
            String headerLine = null;
            try {
                headerLine = reader.readLine();
                String[] tokens = ParsingUtils.TAB_PATTERN.split(headerLine);
                return readHeader(tokens);
            } catch (IOException e) {
                log.error(e.getMessage(), e);
                throw new DataLoadException("Error reading header: " + e.getMessage(), this.path);
            }
        }
    }

    private static class FpkmTrackingCodec extends CufflinksCodec<FPKMValue>{
        int geneColumn = 4;
        int locusColumn  = 6;
        int fpkmColumn = 6;
        int confLoColumn = 7;
        int confHiColumn = 8;
        int sigColumn = 12;

        FpkmTrackingCodec(String path){
            super(FPKMValue.class, path);
        }

        @Override
        public Object readHeader(String[] tokens) {
            for(int i=0; i<tokens.length; i++) {
                final String tk = tokens[i];
                if(tk.equals("locus")) locusColumn = i;
                else if(tk.equals("gene_short_name")) geneColumn = i;
                else if(tk.equals("FPKM")) fpkmColumn = i;
                else if(tk.equals("FPKM_conf_lo")) confLoColumn = i;
                else if(tk.startsWith("FPKM_conf_hi")) confHiColumn = i;
            }
            return tokens;
        }

        @Override
        public FPKMValue decode(String line) {
            return decode(ParsingUtils.TAB_PATTERN.split(line));
        }

        @Override
        public FPKMValue decode(String[] tokens) {
            //Skip header line
            if (tokens[0].equalsIgnoreCase("tracking_id") || tokens[geneColumn].equalsIgnoreCase("gene_short_name")) {
                return null;
            }
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

    private static class ExpDiffCodec extends CufflinksCodec<ExpDiffValue>{

        private static Logger log = Logger.getLogger(ExpDiffCodec.class);

        int geneColumn = 1;
        int locusColumn = 2;
        int xColumn = 6;
        int yColumn = 7;
        int logRatioColumn = 8;
        int sigColumn = 12;

        protected ExpDiffCodec(String path) {
            super(ExpDiffValue.class, path);
        }

        @Override
        public Object readHeader(String[] tokens) {
            for (int i = 0; i < tokens.length; i++) {
                final String tk = tokens[i];
                if (tk.equals("locus")) locusColumn = i;
                else if (tk.equals("gene")) geneColumn = i;
                else if (tk.equals("value_1")) xColumn = i;
                else if (tk.equals("value_2")) yColumn = i;
                else if (tk.startsWith("log2(")) logRatioColumn = i;
                else if (tk.equals("significant")) sigColumn = i;
            }
            return tokens;
        }

        @Override
        public ExpDiffValue decode(String[] tokens) {
            //Skip header line
            if (tokens[0].equalsIgnoreCase("test_id") || tokens[geneColumn].equalsIgnoreCase("gene_id")) {
                return null;
            }
            if (tokens.length >= sigColumn) {
                String locusString = tokens[locusColumn];
                if (locusString == null) return null;

                Locus locus = new Locus(locusString);
                if (locus.getChr() == null) return null;

                String logRatioStr = tokens[logRatioColumn];
                float logRatio = Float.parseFloat(logRatioStr);

                if (Float.isInfinite(logRatio) || Float.isNaN(logRatio)) {
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

        @Override
        public ExpDiffValue decode(String line) {
            return decode(ParsingUtils.TAB_PATTERN.split(line));
        }


    }
}
