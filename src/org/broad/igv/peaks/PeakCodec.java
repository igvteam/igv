/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.UCSCCodec;


/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 20, 2009
 * Time: 10:15:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class PeakCodec extends UCSCCodec<Peak> {

    Genome genome;

    public PeakCodec(Genome genome) {
        super(Peak.class);
        this.genome = genome;
    }

    // Declare a static array once, to be reused.

    public Peak decode(String nextLine) {

        if (nextLine.trim().length() == 0 || nextLine.startsWith("#") || nextLine.startsWith("track") ||
                nextLine.startsWith("browser")) {
            return null;
        }

        String [] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);
        int tokenCount = tokens.length;

        String chrToken = tokens[0];
        String chr = genome == null ? chrToken : genome.getChromosomeAlias(chrToken);
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

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     * <p/>
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     * <p/>
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(String path) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }


}