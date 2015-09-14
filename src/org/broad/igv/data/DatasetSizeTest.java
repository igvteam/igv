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
        @Override
        public String getContig() {
            return null;
        }

        public int getStart() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getEnd() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}
