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

package org.broad.igv.feature.tribble;

import com.google.common.base.Function;
import com.google.common.base.Supplier;
import org.apache.commons.lang.StringUtils;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Dec-11
 */
@Ignore("No Runnable Tests")
public class IGVBEDCodecTest extends AbstractHeadlessTest {

    @Ignore("We always ignore this test anyway, performance not consistent")
    @Test
    public void decodeSpeedTest() throws Exception {
        long benchTime = TestUtils.getBenchmarkTime();

        int nTrials = 5000000;
        BEDStringSupplier supplier = new BEDStringSupplier(nTrials);
        final IGVBEDCodec codec = new IGVBEDCodec();

        Function<String, Void> decodeFunc = new Function<String, Void>(){
            @Override
            public Void apply(String input) {
                BasicFeature feat = codec.decode(input);
                return null;
            }
        };

        supplier.reset();
        String jVersion = System.getProperty(Globals.JAVA_VERSION_STRING);
        System.out.println("\nIGVBEDCodec.decode. java version " + jVersion);
        long[] times = TestUtils.timeMethod(supplier, decodeFunc, nTrials);
        //Calculate average (in nanoseconds)
        double average = TestUtils.average(times);
        long median = times[times.length / 2];


        //we are somewhat forgiving, for the sake of portability. Less than 2 uSec okay, even if it
        //breaks benchmark
        int maxMultiplier = 200000;
        int maxMedian = 3000;
        if (jVersion.contains("1.7")) {
            maxMultiplier = 10000;
            maxMedian = 2000;
        }
        assertTrue("Decoding median speed too slow", median < benchTime / maxMultiplier || median < maxMedian);
    }


    public void timeLoadBigFile() throws Exception {
        //File not in repo
        final String path = "GSM288345_Nanog.bed";
        Supplier<String> supplier = new Supplier<String>() {
            @Override
            public String get() {
                return path;
            }
        };

        final TrackLoader loader = new TrackLoader();

        Function<String, Void> loadFileFunc = new Function<String, Void>() {
            @Override
            public Void apply(String input) {
                try {
                    List<Track> newTrack = loader.load(new ResourceLocator(path), genome);
                } catch (Exception e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
                return null;
            }
        };

        TestUtils.timeMethod(supplier, loadFileFunc, 1);
    }

    private class BEDStringSupplier implements Supplier<String> {

        //final String[] testStrings;
        private int counter = 0;


        private final String chr = "chr1";
        private final int maxLength = 500;

        BEDStringSupplier(int nTrials) {
//            testStrings = new String[nTrials];
//            for(int ft=0; ft < nTrials; ft++){
//                testStrings[ft] = gen(ft);
//            }
        }

        private String gen(int start) {
            String end = "" + (start + (int) Math.random() * maxLength);
            String[] dat = {chr, "" + start, end, "0", "0", "+"};
            return StringUtils.join(dat, '\t');
        }


        @Override
        public String get() {
            return gen(counter++);
            //return testStrings[counter++];
        }

        public void reset() {
            counter = 0;
        }

    }
}
