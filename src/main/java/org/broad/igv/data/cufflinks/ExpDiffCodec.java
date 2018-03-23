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

package org.broad.igv.data.cufflinks;

import org.apache.log4j.Logger;
import org.broad.igv.feature.Locus;
import org.broad.igv.util.ParsingUtils;

/**
 * Codec for cufflinks .exp_diff files
* @author jacob
* @date 2013-Apr-18
*/
public class ExpDiffCodec extends CufflinksCodec<ExpDiffValue>{

    private static Logger log = Logger.getLogger(ExpDiffCodec.class);

    int geneColumn = 1;
    int locusColumn = 2;
    int xColumn = 6;
    int yColumn = 7;
    int logRatioColumn = 8;
    int sigColumn = 12;

    public ExpDiffCodec(String path) {
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

    //@Override
    public ExpDiffValue decode(String[] tokens) {
        //Skip header line
        if (tokens[0].equalsIgnoreCase("test_id") || tokens[geneColumn].equalsIgnoreCase("gene_id")) {
            return null;
        }
        if (tokens.length >= sigColumn) {
            String locusString = tokens[locusColumn];
            if (locusString == null) return null;

            Locus locus = Locus.fromString(locusString);
            if (locus == null || locus.getChr() == null) return null;

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
            return new ExpDiffValue(locus.getChr(), locus.getStart() - 1, locus.getEnd(), gene,
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


    @Override
    public boolean canDecode(String path) {
        String s = path.toLowerCase();
        if(s.endsWith(".gz")) s = s.substring(0, s.length()-3);
        return s.endsWith("gene_exp.diff") || s.endsWith("cds_exp.diff");
    }
}
