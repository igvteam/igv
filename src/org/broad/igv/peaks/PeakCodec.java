/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

package org.broad.igv.peaks;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.tribble.UCSCCodec;
import org.broad.igv.util.ParsingUtils;

import java.awt.*;


/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 20, 2009
 * Time: 10:15:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class PeakCodec extends UCSCCodec {

    // Declare a static array once, to be reused.

    public Peak decode(String nextLine) {

        if (nextLine.trim().length() == 0 || nextLine.startsWith("#") || nextLine.startsWith("track") ||
                nextLine.startsWith("browser")) {
            return null;
        }

        int tokenCount = ParsingUtils.splitWhitespace(nextLine, tokens);

        String[] tokens = nextLine.split("\t");
        String chr = tokens[0];
        int start = Integer.parseInt(tokens[1]);
        int end = Integer.parseInt(tokens[2]);
        String name = tokens[3];
        float combinedScore = Float.parseFloat(tokens[4]);

        int nTimePoints = tokenCount - 5;

        float[] timePointScores = new float[nTimePoints];
        for (int i = 0; i < nTimePoints; i++) {
            timePointScores[i] = Float.parseFloat(tokens[5 + i]);
        }
        return new Peak(chr, start, end, name, combinedScore, timePointScores);

    }


}