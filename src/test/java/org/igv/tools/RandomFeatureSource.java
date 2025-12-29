package org.igv.tools;

/**
 * User: jacob
 * Date: 2013-Mar-13
 */

import org.igv.feature.LocusScore;
import org.igv.feature.tribble.Locus;
import org.igv.track.FeatureSource;
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

}
