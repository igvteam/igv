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

    @Override
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
