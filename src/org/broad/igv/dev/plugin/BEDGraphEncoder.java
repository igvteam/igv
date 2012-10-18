/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.dev.plugin;

import org.apache.commons.lang.StringUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.tribble.Feature;

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
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
        //pass
    }
}
