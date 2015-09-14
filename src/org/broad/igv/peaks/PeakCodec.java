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