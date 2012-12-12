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

import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.tribble.IGVBEDCodec;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Codec for reading and writing data with BEDTools
 * User: jacob
 * Date: 2012-Aug-02
 */
public final class BEDToolsDecoder extends AsciiDecoder<BasicFeature> implements LineFeatureDecoder<BasicFeature> {

    private IGVBEDCodec BEDCodec = new IGVBEDCodec();
    private boolean hasSplit = false;

    private int numTracks = 0;

    private int[] numCols = null;

    private boolean multiinter;

    /**
     * The output of the closest and window
     * commands are somewhat different from other commands
     */
    private boolean closestOrSim;

    public BEDToolsDecoder() {
        super.lineFeatureDecoder = this;
    }

    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
        super.setInputs(commands, argumentMap);
        for (String cmdTok : commands) {
            hasSplit |= cmdTok.contains("-split");
        }
        multiinter = commands.contains("multiinter");
        closestOrSim = commands.contains("window") || commands.contains("closest");

        for (Map.Entry<Argument, Object> entry : argumentMap.entrySet()) {
            Argument argument = entry.getKey();
            hasSplit |= argument.getCmdArg().contains("-split");

            switch (argument.getType()) {
                case TEXT:
                    String sVal = (String) entry.getValue();
                    if (sVal == null) continue;
                    hasSplit |= (sVal).contains("-split");
                    break;
                case FEATURE_TRACK:
                    numTracks++;
                    break;
                case MULTI_FEATURE_TRACK:
                    numTracks += ((List) entry.getValue()).size();
                    break;
            }

        }
    }

    @Override
    public BasicFeature decode(String line) {
        return this.decode(Globals.singleTabMultiSpacePattern.split(line));
    }

    public BasicFeature decode(String[] tokens) {
        BasicFeature feat;

        if (hasSplit) {
            //When we split, the returned feature still has the exons
            //We don't want to plot them all a zillion times
            tokens = Arrays.copyOfRange(tokens, 0, Math.min(6, tokens.length));
        }

        if (closestOrSim) {
            String[] closest = Arrays.copyOfRange(tokens, numCols[0], numCols[0] + numCols[1]);
            //If not found, bedtools returns -1 for positions
            if (closest[1].trim().equalsIgnoreCase("-1")) {
                return null;
            }
            feat = BEDCodec.decode(closest);
        } else if (multiinter) {
            //We only look at regions common to ALL inputs
            //Columns: chr \t start \t \end \t # of files which contained this feature \t comma-separated list files +many more
            int numRegions = Integer.parseInt(tokens[3]);
            if (numRegions < numTracks) {
                return null;
            }
            String[] intersection = Arrays.copyOf(tokens, 3);
            feat = BEDCodec.decode(intersection);
        } else {
            feat = BEDCodec.decode(tokens);
        }
        return feat;
    }

    @Override
    public void setAttributes(List<Map<String, Object>> attributes) {
        super.setAttributes(attributes);

        numCols = new int[attributes.size()];
        int ind = 0;
        for (Map<String, Object> attributeMap : attributes) {
            int curOutCols = (Integer) attributeMap.get(AsciiEncoder.NUM_COLS_ATTR);
            numCols[ind++] = curOutCols;
        }
    }
}
