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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.maf.conservation;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ParsingUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

/**
 * @author jrobinso
 */
public class OmegaFileParser {

    public void getDataAsScores(File f, List<LocusScore> scores) {


        AsciiLineReader reader = null;
        try {
            String[] buffer = new String[5];
            reader = new AsciiLineReader(new FileInputStream(f));
            String nextLine = null;
            while ((nextLine = reader.readLine()) != null) {
                ParsingUtils.split(nextLine, buffer, '\t');
                int position = Integer.parseInt(buffer[0]) + 6;
                float score = Float.parseFloat(buffer[4]);
                scores.add(new OmegaScore(position, score));
            }

        } catch (IOException ex) {
            throw new RuntimeException("Error parsing file: " + f.getAbsolutePath(), ex);

        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    public void getDataAsArrays(File f, int[] locations, float[] data, int offset, int start, int end) {

        AsciiLineReader reader = null;
        try {
            String[] buffer = new String[5];
            reader = new AsciiLineReader(new FileInputStream(f));
            String nextLine = null;
            int idx = offset;
            while ((nextLine = reader.readLine()) != null) {
                ParsingUtils.split(nextLine, buffer, '\t');
                int s = Integer.parseInt(buffer[0]) + 6;
                if (s >= start) {
                    locations[idx] = s;
                    data[idx] = Float.parseFloat(buffer[4]);
                } else if (s > end) {
                    return;
                }
            }

        } catch (IOException ex) {
            throw new RuntimeException("Error parsing file: " + f.getAbsolutePath(), ex);

        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    public static class OmegaScore implements LocusScore {

        int position;
        float score;

        public OmegaScore(int position, float score) {
            this.position = position;
            this.score = score;
        }

        public String getChr() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getStart() {
            return position;
        }

        public void setStart(int start) {
            throw new UnsupportedOperationException("setStart is not supported.");
        }

        public int getEnd() {
            return position + 1;
        }

        public void setEnd(int end) {
            throw new UnsupportedOperationException("setEnd is not supported");
        }

        public float getScore() {
            return score;
        }

        public void setConfidence(float confidence) {
            throw new UnsupportedOperationException("setConfidence is not supported yet.");
        }

        public float getConfidence() {
            return 1.0f;
        }

        public LocusScore copy() {
            return new OmegaScore(position, score);
        }

        public String getValueString(double position, WindowFunction windowFunction) {
            return String.valueOf("Omega score = " + score);
        }
    }
}
