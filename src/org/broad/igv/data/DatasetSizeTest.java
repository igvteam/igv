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

package org.broad.igv.data;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.RuntimeUtils;

import java.util.ArrayList;

/**
 * @author jrobinso
 * @date Mar 20, 2011
 */
public class DatasetSizeTest {

    static int nPoints = 100000;
    static int nSample = 1000;

    public static void main(String[] args) {
        for (int i = 0; i < 5; i++) {
            testArrays();
        }
        for (int i = 0; i < 5; i++) {
            testLocusScores();
        }
    }


    static void testArrays() {


        System.gc();

        int[] starts = new int[nPoints];
        int[] ends = new int[nPoints];
        float[][] values = new float[nPoints][nSample];
        String[] probes = new String[nPoints];

        for (int i = 0; i < nPoints; i++) {
            starts[i] = i;
            ends[i] = i;
            for (int j = 0; j < nSample; j++) {
                values[i][j] = i * j;
            }

            probes[i] = "probe" + i;
        }
        long mem = RuntimeUtils.getAvailableMemory(); //Runtime.getRuntime().freeMemory();
        starts = ends = null;
        values = null;
        probes = null;
        System.gc();
        System.gc();
        long mem2 = RuntimeUtils.getAvailableMemory(); //Runtime.getRuntime().freeMemory();
        System.out.println("Memory used = " + (mem2 - mem) / 1000000);


    }


    static void testLocusScores() {


        System.gc();

        ArrayList<TestScore> scores = new ArrayList(nPoints);

        for (int i = 0; i < nPoints; i++) {
            for (int j = 0; j < nSample; j++)
                scores.add(new TestScore(i, "probe" + i, i, i));
        }
        long mem = RuntimeUtils.getAvailableMemory(); //Runtime.getRuntime().freeMemory();
        scores = null;
        System.gc();
        System.gc();
        long mem2 = RuntimeUtils.getAvailableMemory(); //Runtime.getRuntime().freeMemory();
        System.out.println("Memory used = " + (mem2 - mem) / 1000000);


    }

    static class TestScore implements LocusScore {

        int start;
        int end;
        float score;
        String name;

        TestScore(int end, String name, float score, int start) {
            this.end = end;
            this.name = name;
            this.score = score;
            this.start = start;
        }

        public void setStart(int start) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public void setEnd(int end) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public float getScore() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public LocusScore copy() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public String getValueString(double position, WindowFunction windowFunction) {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public String getChr() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getStart() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getEnd() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}
