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

package org.broad.igv.feature.tribble;

import com.google.common.base.Predicate;
import com.google.common.base.Supplier;
import org.apache.commons.lang.StringUtils;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Dec-11
 */
public class IGVBEDCodecTest extends AbstractHeadlessTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 30e5);


    @Test
    public void decodeSpeedTest() throws Exception {
        long benchTime = TestUtils.getBenchmarkTime();

        int nTrials = 5000000;
        BEDStringSupplier supplier = new BEDStringSupplier(nTrials);
        final IGVBEDCodec codec = new IGVBEDCodec();

        Predicate<String> decodePredicate = new Predicate<String>() {
            @Override
            public boolean apply(String input) {
                BasicFeature feat = codec.decode(input);
                return true;
            }
        };

        supplier.reset();
        String jVersion = System.getProperty(Globals.JAVA_VERSION_STRING);
        System.out.println("\nIGVBEDCodec.decode. java version " + jVersion);
        long[] times = TestUtils.timeMethod(supplier, decodePredicate, nTrials);
        //Calculate average (in nanoseconds)
        double average = TestUtils.average(times);

        int maxMultiplier = 200000;
        if (jVersion.contains("1.7")) {
            maxMultiplier = 10000;
        }
        //we are somewhat forgiving, for the sake of portability. Less than 2 uSec okay, even if it
        //breaks benchmark
        assertTrue("Decoding median speed too slow", average < benchTime / maxMultiplier || average < 2000);
    }


    public void timeLoadBigFile() throws Exception {
        final String path = "/Users/jacob/Data/GSM288345_Nanog.bed";
        Supplier<String> supplier = new Supplier<String>() {
            @Override
            public String get() {
                return path;
            }
        };

        final TrackLoader loader = new TrackLoader();

        Predicate<String> loadFile = new Predicate<String>() {
            @Override
            public boolean apply(String input) {
                List<Track> newTrack = loader.load(new ResourceLocator(path), genome);
                return true;
            }
        };

        TestUtils.timeMethod(supplier, loadFile, 1);
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
