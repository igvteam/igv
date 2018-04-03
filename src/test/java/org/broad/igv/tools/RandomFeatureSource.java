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

package org.broad.igv.tools;

/**
 * User: jacob
 * Date: 2013-Mar-13
 */

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.tribble.Locus;
import org.broad.igv.track.FeatureSource;
import org.junit.Ignore;

import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 * Generates a feature in the search interval, or null, with some
 * probability of a feature
 */
@Ignore
public class RandomFeatureSource implements FeatureSource<Locus> {

    private float probSuccess = 0.05f;
    private Random generator = new Random();

    public void setGenerator(Random generator) {
        this.generator = generator;
    }

    public void setProbSuccess(float probSuccess) {
        this.probSuccess = probSuccess;
    }

    @Override
    public Iterator<Locus> getFeatures(String chr, int start, int end) throws IOException {
        //System.out.println(String.format("%s:%d-%d", chr, start, end));
        if(generator.nextFloat() <= probSuccess){
            //System.out.println("success");
            Locus feature = new Locus(chr, start, Math.min(start + 50, end));
            return Arrays.asList(feature).iterator();
        }else{
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        return null;
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    @Override
    public int getFeatureWindowSize() {
        return 0;
    }

    @Override
    public void setFeatureWindowSize(int size) {

    }
}
