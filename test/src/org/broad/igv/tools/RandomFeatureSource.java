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
