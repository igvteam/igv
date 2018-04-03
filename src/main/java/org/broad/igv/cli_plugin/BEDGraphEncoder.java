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

package org.broad.igv.cli_plugin;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.util.StringUtils;
import htsjdk.tribble.Feature;

import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Encodes a feature into the bedgraph format.
 * We separate columns with a single tab.
 * <p/>
 * We only encode LocusScore objects, anything else is skipped
 * User: jacob
 * Date: 2012-Aug-23
 */
public class BEDGraphEncoder implements LineFeatureEncoder {

    protected final Pattern splitter = Pattern.compile("\\s+");
    protected String delimiter = "\t";

    @Override
    public String encode(Feature feature) {
        if (feature instanceof LocusScore) {
            return encode((LocusScore) feature);
        } else {
            return null;
        }
    }

    public String encode(LocusScore score) {
        String[] tokens = new String[4];
        tokens[0] = score.getChr();
        tokens[1] = "" + score.getStart();
        tokens[2] = "" + score.getEnd();
        tokens[3] = "" + score.getScore();
        String out = StringUtils.join(tokens, delimiter);
        return out;
    }

    @Override
    public int getNumCols(String line) {
        return splitter.split(line).length;
    }

    @Override
    public String getHeader() {
        return "track type=bedGraph";
    }

    /**
     *
     *
     */
    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap, Argument argument) {
        //pass
    }
}
