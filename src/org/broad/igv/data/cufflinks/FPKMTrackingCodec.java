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
 * Codec for Cufflinks FPKM files, extension fpkm_tracking
 *
* @author jacob
* @date 2013-Apr-18
*/
public class FPKMTrackingCodec extends CufflinksCodec<FPKMValue>{

    private static Logger log = Logger.getLogger(FPKMTrackingCodec.class);

    int geneColumn = 4;
    int locusColumn  = 6;

    static final int startfpkmCol = 9;
    static final int colsPerSample = 4;

    private int numSamples = 1;

    public FPKMTrackingCodec(String path){
        super(FPKMValue.class, path);
    }

    @Override
    public Object readHeader(String[] tokens) {
        for(int i=0; i<tokens.length; i++) {
            final String tk = tokens[i];
            if(tk.equals("locus")) locusColumn = i;
            else if(tk.equals("gene_short_name")) geneColumn = i;
        }
        numSamples = (tokens.length - startfpkmCol) / colsPerSample;
        return tokens;
    }

    @Override
    public FPKMValue decode(String line) {
        return decode(ParsingUtils.TAB_PATTERN.split(line));
    }

    //@Override
    public FPKMValue decode(String[] tokens) {
        //Skip header line
        if (tokens[0].equalsIgnoreCase("tracking_id") || tokens[geneColumn].equalsIgnoreCase("gene_short_name")) {
            return null;
        }
        if (tokens.length >= (startfpkmCol + numSamples*colsPerSample)) {
            String locusString = tokens[locusColumn];
            if (locusString == null) return null;

            Locus locus = Locus.fromString(locusString);
            if(locus == null || locus.getChr() == null) return null;

            String gene = tokens[geneColumn];
            float[] fpkm = new float[numSamples];
            float[] confLo = new float[numSamples];
            float[] confHi = new float[numSamples];

            for(int sampNum = 0; sampNum < numSamples; sampNum++){
                int startCol = startfpkmCol + sampNum*colsPerSample;
                fpkm[sampNum] = Float.parseFloat(tokens[startCol]);
                confLo[sampNum] = Float.parseFloat(tokens[startCol+1]);
                confHi[sampNum] = Float.parseFloat(tokens[startCol+2]);
            }

            return new FPKMValue(locus.getChr(), locus.getStart() - 1, locus.getEnd(), gene,
                    fpkm, confLo, confHi);
        } else {
            log.info("Unexpected # of columns.  Expected at least 12,  found " + tokens.length);
            return null;
        }
    }

    public int getNumSamples() {
        return numSamples;
    }

    @Override
    public boolean canDecode(String path) {

        String fn = path.toLowerCase();
        if(fn.endsWith(".gz")) fn = fn.substring(0, fn.length()-3);
        return fn.endsWith("fpkm_tracking");
    }
}
